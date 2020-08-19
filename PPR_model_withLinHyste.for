c ====================================================================
c Two-Dimensional Linear Cohesive Zone Element for the Intrinsic Model
c Developed by Kyoungsoo Park and Glaucio H. Paulino
c
c Kyoungsoo Park
c School of Civil and Environmental Engineering
c Yonsei University, Seoul, S. Korea
c
c Glaucio H. Paulino
c Department of Civil and Environmental Engineering
c University of Illinois at Urbana-Champaign, U.S.A.
c
c References
c
c K. Park, and G.H. Paulino, 2012, Implementation of the PPR potential-
c    based model in ABAQUS: Educational perspective, Engineering 
c    Fracture Mechanics XX(X), XXX-XXX.
c    DOI: http://dx.doi.org/
c K. Park, G.H. Paulino, and J.R. Roesler, 2009, A unified potential-
c    based cohesive model of mixed-mode fracture, Journal of the 
c    Mechanics and Physics of Solids 57 (6), 891-908.
c    DOI: http://dx.doi.org/10.1016/j.jmps.2008.10.003
c K. Park, 2009, Potential-based fracture mechanics using cohesive zone
c    and virtual internal bond modeling, PhD Thesis, University of
c    Illinois at Urbana-Champaign.
c     reloading degraded by a simple factor 
C
c SUBROUTINE UEL
c  : Main function for the computational implementation of
c    a two-dimensional linear intrinsic cohesive element
c SUBROUTINE k_Cohesive_PPR
c  : Function to compute the traction-separation relationship
c    of the PPR potential-based model
c SUBROUTINE k_Coords_Transform
c  : Function to compute the coordinate transformation matrix
c    between the global and local coordinates
c SUBROUTINE k_Matrix_Zero       : Matrix operation (A = 0)
c SUBROUTINE k_Matrix_Transpose  : Matrix operation (B = A_t)
c SUBROUTINE k_Matrix_PlusScalar : Matrix operation (A = A + c * B)
c SUBROUTINE k_Matrix_Multiply   : Matrix operation (C = A * B)
c
c =======================================================================
       SUBROUTINE UEL (RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS,
     & PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME,
     & DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG,
     & PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, JPROPS,
     & NJPRO, PERIOD)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION RHS(MLVARX,*), AMATRX(NDOFEL,NDOFEL), PROPS(*),
     & SVARS(*), ENERGY(8), COORDS(MCRD, NNODE), U(NDOFEL),
     & DU(MLVARX,*), V(NDOFEL), A(NDOFEL), TIME(2), PARAMS(*),
     & JDLTYP(MDLOAD,*), ADLMAG(MDLOAD,*), DDLMAG(MDLOAD,*),
     & PREDEF(2, NPREDF, NNODE), LFLAGS(*), JPROPS(*)
c Variables defined in the UEL subroutine
c   RHS   : Right-Hand-Side vector
c   AMATRX: Stiffness (Jacobian) matrix
c
c Variables updated in the UEL subroutine
c   SVARS : Maximum separation at each integration point
c
c Variables available in the UEL subroutine
c   U     : Nodal displacement
c   COORDS: Nodal coordinates of an element
c   PROPS : Parameters from an input file
c      PROPS(1): Normal fracture energy (Gn)
c      PROPS(2): Tangential fracture energy (Gt)
c      PROPS(3): Normal cohesive strength (Tn_m)
c      PROPS(4): Tangential cohesive strength (Tt_m)
c      PROPS(5): Normal shape parameter (alph)
c      PROPS(6): Tangential shape parameter (beta)
c      PROPS(7): Normal initial slope indicator (ln)
c      PROPS(8): Tangential initial slope indicator (lt)
c      PROPS(9): Thickness of a cohesive element
c   MCRD  : Largest active degrees of freedom
c   NNODE : Number of nodes
c
c Variables used in the UEL subroutine
       DIMENSION Sc(ndofel,ndofel), Fc(ndofel,nrhs),
     & T(mcrd,nrhs), T_d(mcrd,mcrd), U_l(ndofel), R(mcrd, mcrd),
     & Bc(mcrd,ndofel), Bct(ndofel,mcrd), ShapeN(nnode),
     & del(mcrd), GP(2), GP_w(2), tmp(ndofel,mcrd)
c   Sc  : Element stiffness matrix of a cohesive element
c   Fc  : Cohesive internal force vector
c   T   : Cohesive traction vector
c   T_d : Derivative of the cohesive traction (Tangent matrix)
c   U_l : Nodal displacement in the local coordinate system
c   R   : Coordinate transformation matrix
c   Bc  : Global displacement-separation relation matrix
c   Bct : Transpose of Bc
c   ShapeN: Shape functional matrix
c   del : Normal and tangential separations
c   GP  : Gauss points
c   GP_W: Weight at the Gauss points
c   n_GP: Number of the Gauss points
       DOUBLE PRECISION Gn, Gt, Tn_m, Tt_m, alph, beta, ln, lt, th,
     & dn, dt, m, n, Gam_n, Gam_t, dGnt, dGtn,
     & N1, N2, del1, del2, del3, del4, deln_max, delt_max, el_length,
     & deln_pr,delt_pr
c   Gn, Gt: Fracture energies
c   Tn_m, Tt_m: Cohesive strengths
c   alph, beta: Shape parameters
c   ln, lt: Initial slope indicators
c   th    : Thickness of a cohesive element
c   dn, dt: Final crack opening widths
c   m, n  : Exponents in the PPR potential
c   Gam_n, Gam_t: Energy constants in the PPR potential
c   dGnt  : <Gn - Gt>
c   dGtn  : <Gt - Gn>
c   N1, N2: Linear shape functions
c   del1, del2, del3, del4: Nodal separations
c   deln_max, delt_max: Maximum separations in a loading history
c   el_length: Length of a cohesive element
c   deln_pr,delt_pr: 
c -----------------------------------------------------------------------
c Read input data & Initialize
       Gn   = PROPS(1)
       Gt   = PROPS(2)
       Tn_m = PROPS(3)
       Tt_m = PROPS(4)
       alph = PROPS(5)
       beta = PROPS(6)
       ln   = PROPS(7)
       lt   = PROPS(8)
       th   = PROPS(9)
       n_GP = 2
C      1/sqrt(3) 
       data GP   / 0.577350269189626 , -0.577350269189626 /
       data GP_W / 1.0 , 1.0 /
       call k_Matrix_Zero (RHS,ndofel,nrhs)
       call k_Matrix_Zero (AMATRX,ndofel,ndofel)
c Determine the PPR parameters
c   From Equation (20):
       m = (alph-1)*alph*ln**2/(1-alph*ln**2)
       n = (beta-1)*beta*lt**2/(1-beta*lt**2)
c   From Equation (21):
       dn = alph*Gn/(m*Tn_m)*(1-ln)**(alph-1)
     &     * (alph/m*ln+1)**(m-1)*(alph+m)*ln
C       write(*,*) 'dn',dn
c   From Equation (22):
       dt = beta*Gt/(n*Tt_m)*(1-lt)**(beta-1)
     &     * (beta/n*lt+1)**(n-1)*(beta+n)*lt
C       write(*,*) 'dt',dt
c   From Equation (11):
       if (Gt .GT. Gn) then
          dGnt = 0
          dGtn = Gt - Gn
       elseif (Gt .LT. Gn) then
          dGnt = Gn - Gt
          dGtn = 0
       else
          dGnt = 0
          dGtn = 0
       endif
c   From Equations (18) and (19):
       if (Gn .EQ. Gt) then
          Gam_n = -Gn*(alph/m)**m
          Gam_t = -Gt*(beta/n)**n
       else
          Gam_n = (-Gn)**(dGnt/(Gn-Gt))*(alph/m)**m
          Gam_t = (-Gt)**(dGtn/(Gt-Gn))*(beta/n)**n
       endif
c Change from the global coordinates to the local coordinates
       call k_Coords_Transform (R, el_length, COORDS, U, ndofel,
     & nnode, mcrd)
       do i = 0, nnode-1
          U_l(1+i*mcrd) = R(1,1)*U(1+i*mcrd) + R(1,2)*U(2+i*mcrd)
          U_l(2+i*mcrd) = R(2,1)*U(1+i*mcrd) + R(2,2)*U(2+i*mcrd)
       end do
       del1 = U_l(7) - U_l(1)
       del2 = U_l(8) - U_l(2)
       del3 = U_l(5) - U_l(3)
       del4 = U_l(6) - U_l(4)
c Numerical integration to compute RHS and AMATRX
       do i = 1, n_GP
          N1 = 0.5*(1 - GP(i))
          N2 = 0.5*(1 + GP(i))
          del(1) = N1*del1 + N2*del3
          del(2) = N1*del2 + N2*del4
          delt_max = SVARS(n_GP*(i-1)*2+1)
          deln_max = SVARS(n_GP*(i-1)*2+2)
          deln_pr = SVARS(n_GP*(i-1)*2+3)
          delt_pr = SVARS(n_GP*(i-1)*2+4)
C          write(*,*) 'delt_max',delt_max,del(1)
C          write(*,*) 'deln_max',deln_max,del(2)
          call k_Cohesive_PPR (T, T_d, Gam_n, Gam_t, alph, beta, m, n,
     &  dn, dt, dGtn, dGnt, del, deln_max, delt_max,deln_pr,delt_pr)
          ShapeN(1) = -N1
          ShapeN(2) = -N2
          ShapeN(3) = N2
          ShapeN(4) = N1
          do j = 1, nnode
             do k = 1, mcrd
               do l = 1, mcrd
                  Bc(k,l+(j-1)*mcrd) = ShapeN(j)*R(k,l)
               end do
             end do
          end do
          call k_Matrix_Transpose (Bc,Bct,mcrd,ndofel)
          call k_Matrix_Multiply (Bct,T_d,tmp,ndofel,mcrd,mcrd)
          call k_Matrix_Multiply (tmp,Bc,Sc,ndofel,mcrd,ndofel)
          call k_Matrix_Multiply (Bct,T,Fc,ndofel,mcrd,nrhs)
          thick = 0.5 * el_length * GP_w(i) * th
          call k_Matrix_PlusScalar (AMATRX,Sc,thick,ndofel,ndofel)
          call k_Matrix_PlusScalar (RHS,-Fc,thick,ndofel,nrhs)
c   Update the state variables: SVARS
          if((delt_max.LT.abs(del(1))).AND.(abs(del(1)).GT.lt*dt)) then
             SVARS(n_GP*(i-1)*2+1) = abs(del(1))
          end if
          if ((deln_max .LT. del(2)) .AND. (del(2) .GT. ln*dn)) then
             SVARS(n_GP*(i-1)*2+2) = del(2)
          end if
          deln_pr =abs(del(1))
          delt_pr =del(2)
       end do
       RETURN
       END
c
c ====================================================================
c = Cohesive traction-separation relation of the PPR model ===========
       SUBROUTINE k_Cohesive_PPR (T,T_d,Gam_n,Gam_t,alph,beta,m,n,
     &  dn, dt, dGtn, dGnt, del, deln_max, delt_max,deln_pr,delt_pr)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION T(2,1), T_d(2,2), del(2)
       DOUBLE PRECISION Gam_n, Gam_t, alph, beta, m, n, dn, dt,
     &     dGtn, dGnt, deln_max, delt_max, Tn, Tt, deln, delt, sign_dt
     & ,deln_pr,delt_pr, alphaDeg
       alphaDeg  = 0.9 
       delt = abs(del(1))
       deln = del(2)
       if (del(1) .GE. 0) then
          sign_dt = 1
       else
          sign_dt = -1
       end if
       Tn = 0
c Equation (12):
c Pre-calculation of the normal cohesive traction, Tn
       if (deln .LT. 0) then
          deln = 0
       elseif ((deln .GE. dn) .OR. (delt .GE. dt))  then
          Tn = 0
       elseif (deln .GE. deln_max) then ! softening 
          Tn = (Gam_t*(1-delt/dt)**beta*(delt/dt+n/beta)**n+dGtn) *
     & Gam_n/dn*(m*(1-deln/dn)**alph*(m/alph+deln/dn)**(m-1)
     &       -alph*(1-deln/dn)**(alph-1)*(m/alph+deln/dn)**m)
       elseif (deln_pr .GE. deln) then ! unloading 
          Tn = (Gam_t*(1-delt/dt)**beta*(delt/dt+n/beta)**n+dGtn) *
     & Gam_n/dn*(m*(1-deln_max/dn)**alph*(m/alph+deln_max/dn)**(m-1)
     &   -alph*(1-deln_max/dn)**(alph-1)*(m/alph+deln_max/dn)**m) *
     & deln/deln_max
       else !reloading 
          Tn = ((Gam_t*(1-delt/dt)**beta*(delt/dt+n/beta)**n+dGtn) *
     & Gam_n/dn*(m*(1-deln_max/dn)**alph*(m/alph+deln_max/dn)**(m-1)
     &   -alph*(1-deln_max/dn)**(alph-1)*(m/alph+deln_max/dn)**m) *
     & deln/deln_max)*alphaDeg
       end if
c Pre-calculation of the tangential cohesive traction, Tt
       if ((deln .GE. dn) .OR. (delt .GE. dt))  then
          Tt = 0
       elseif (delt .GE. delt_max) then ! softening 
          Tt = (Gam_n*(1-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt*(n*(1-delt/dt)**beta*(delt/dt+n/beta)**(n-1)
     &       -beta*(1-delt/dt)**(beta-1)*(delt/dt+n/beta)**n)
       elseif (delt_pr .GE. delt) then ! unloading ! 
          Tt = (Gam_n*(1-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt*(n*(1-delt_max/dt)**beta*(delt_max/dt+n/beta)**(n-1)
     &   -beta*(1-delt_max/dt)**(beta-1)*(delt_max/dt+n/beta)**n) *
     & delt/delt_max
       else ! reloading 
          Tt = ((Gam_n*(1-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt*(n*(1-delt_max/dt)**beta*(delt_max/dt+n/beta)**(n-1)
     &   -beta*(1-delt_max/dt)**(beta-1)*(delt_max/dt+n/beta)**n) *
     & delt/delt_max)*alphaDeg
       end if
c Algorithm 1 (description)
c Normal cohesive interaction
c   (1) Contact condition
       if (del(2). LT. 0) then
          T_d(2,2) = -Gam_n/dn**2*(m/alph)**(m-1)*(alph+m)*
     & (Gam_t*(n/beta)**n + dGtn)
          T_d(2,1) = 0
          T(2,1) = T_d(2,2)*del(2)
       else if ((deln.LT.dn).AND.(delt.LT.dt).AND.(Tn.GE.-1.0E-5)) then
          T(2,1) = Tn
c   (2) Softening condition
          if (deln .GE. deln_max) then
             T_d(2,2) =
     & (Gam_t*(1-delt/dt)**beta*(delt/dt+n/beta)**n+dGtn) *
     & Gam_n/dn**2 *
     & ((1-deln/dn)**(alph-2)*(alph-1)*alph*(deln/dn+m/alph)**m -
     & 2*(1-deln/dn)**(alph-1)*alph*(deln/dn+m/alph)**(m-1)*m +
     & (1-deln/dn)**alph*(deln/dn+m/alph)**(m-2)*(m-1)*m)
             T_d(2,1) =
     & Gam_t/dt*(-(1-delt/dt)**(beta-1)*beta*(delt/dt+n/beta)**n +
     & (1-delt/dt)**beta*(delt/dt+n/beta)**(n-1)*n) * sign_dt *
     & Gam_n/dn*(-(1-deln/dn)**(alph-1)*alph*(deln/dn+m/alph)**m +
     & (1-deln/dn)**alph*(deln/dn+m/alph)**(m-1)*m)
c   (3) Unloading/reloading condition
          elseif (deln_pr .GE. deln) then ! unloading 
             T_d(2,2) =
     & (Gam_t*(1-delt/dt)**beta*(delt/dt+n/beta)**n+dGtn) *
     & Gam_n/dn*((1-deln_max/dn)**alph*(deln_max/dn+m/alph)**(m-1)*m
     & -(1-deln_max/dn)**(alph-1)*alph*(deln_max/dn+m/alph)**m)
     & / deln_max
             T_d(2,1) =
     & Gam_t/dt*(-(1-delt/dt)**(beta-1)*beta*(delt/dt+n/beta)**n +
     & (1-delt/dt)**beta*(delt/dt+n/beta)**(n-1)*n) * sign_dt *
     & Gam_n/dn*(m*(1-deln_max/dn)**alph*(m/alph+deln_max/dn)**(m-1)
     &   -alph*(1-deln_max/dn)**(alph-1)*(m/alph+deln_max/dn)**m) *
     & deln/deln_max
          else ! reloading  
             T_d(2,2) =
     & ((Gam_t*(1-delt/dt)**beta*(delt/dt+n/beta)**n+dGtn) *
     & Gam_n/dn*((1-deln_max/dn)**alph*(deln_max/dn+m/alph)**(m-1)*m
     & -(1-deln_max/dn)**(alph-1)*alph*(deln_max/dn+m/alph)**m)
     & / deln_max)*alphaDeg
             T_d(2,1) =
     & (Gam_t/dt*(-(1-delt/dt)**(beta-1)*beta*(delt/dt+n/beta)**n +
     & (1-delt/dt)**beta*(delt/dt+n/beta)**(n-1)*n) * sign_dt *
     & Gam_n/dn*(m*(1-deln_max/dn)**alph*(m/alph+deln_max/dn)**(m-1)
     &   -alph*(1-deln_max/dn)**(alph-1)*(m/alph+deln_max/dn)**m) *
     & deln/deln_max)*alphaDeg

          end if
c   (4) Complete failure condition
       else
          T(2,1) = 0
          T_d(2,2) = 0
          T_d(2,1) = 0
       endif
c Tangential cohesive interaction
       if ((delt.LT.dt) .AND. (deln.LT.dn) .AND. (Tt.GE.-1.0E-5)) then
          T(1,1) = Tt*sign_dt
c   (1) Softening condition
          if (delt .GE. delt_max) then
             T_d(1,1) =
     & (Gam_n*(1-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt**2 *
     & ((1-delt/dt)**(beta-2)*(beta-1)*beta*(delt/dt+n/beta)**n -
     & 2*(1-delt/dt)**(beta-1)*beta*(delt/dt+n/beta)**(n-1)*n +
     & (1-delt/dt)**beta*(delt/dt+n/beta)**(n-2)*(n-1)*n)
             T_d(1,2) =
     & Gam_t/dt*(-(1-delt/dt)**(beta-1)*beta*(delt/dt+n/beta)**n +
     & (1-delt/dt)**beta*(delt/dt+n/beta)**(n-1)*n) * sign_dt *
     & Gam_n/dn*(-(1-deln/dn)**(alph-1)*alph*(deln/dn+m/alph)**m +
     & (1-deln/dn)**alph*(deln/dn+m/alph)**(m-1)*m)
c   (2) Unloading/reloading condition
          elseif (delt_pr .GE. delt) then ! unloading ! 
             T_d(1,1) =
     & (Gam_n*(1-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt*(n*(1-delt_max/dt)**beta*(delt_max/dt+n/beta)**(n-1)
     &      -beta*(1-delt_max/dt)**(beta-1)*(delt_max/dt+n/beta)**n)
     & / delt_max
             T_d(1,2) =
     & Gam_n/dn*(-(1-deln/dn)**(alph-1)*alph*(deln/dn+m/alph)**m +
     & (1-deln/dn)**alph*(deln/dn+m/alph)**(m-1)*m) * sign_dt *
     & Gam_t/dt*(n*(1-delt_max/dt)**beta*(delt_max/dt+n/beta)**(n-1)
     &   -beta*(1-delt_max/dt)**(beta-1)*(delt_max/dt+n/beta)**n) *
     & delt/delt_max
          else !relaoding 
             T_d(1,1) =
     & ((Gam_n*(1-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt*(n*(1-delt_max/dt)**beta*(delt_max/dt+n/beta)**(n-1)
     &      -beta*(1-delt_max/dt)**(beta-1)*(delt_max/dt+n/beta)**n)
     & / delt_max)*alphaDeg
             T_d(1,2) =
     & (Gam_n/dn*(-(1-deln/dn)**(alph-1)*alph*(deln/dn+m/alph)**m +
     & (1-deln/dn)**alph*(deln/dn+m/alph)**(m-1)*m) * sign_dt *
     & Gam_t/dt*(n*(1-delt_max/dt)**beta*(delt_max/dt+n/beta)**(n-1)
     &   -beta*(1-delt_max/dt)**(beta-1)*(delt_max/dt+n/beta)**n) *
     & delt/delt_max)*alphaDeg
          end if
c   (3) Complete failure condition
       else
          T(1,1) = 0
          T_d(1,1) = 0
          T_d(1,2) = 0
       endif
c      if (T_d(1,2) .NE. T_d(2,1)) then
c         T_d(1,2) = 0.5*(T_d(1,2) + T_d(2,1))
c         T_d(2,1) = T_d(1,2)
c      endif
       RETURN
       END
c
c =====================================================================
c = Coordinate Transformation =========================================
c   : Coordinate transformation matrix (R) is obtained on the basis of
c   the deformed configuration
       SUBROUTINE k_Coords_Transform (R, el_length, COORDS, U, ndofel,
     & nnode, mcrd)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION R(mcrd,mcrd), COORDS(mcrd,nnode), U(ndofel)
       DIMENSION Co_de(mcrd,nnode), Co_de_m(2,2)
c Variables used in the k_Coords_Transform subroutine
c   Co_de  : Coord. of a cohesive element in the deformed configuration
c   Co_de_m: Mid-points of a cohesive element to compute the orientation
c   el_length: length of a cohesive element
c
       do i = 1, mcrd
          do j = 1, nnode
             Co_de(i,j) = COORDS(i,j) + U(2*(j-1)+i)
          end do
       end do
       do i = 1, 2
          Co_de_m(i,1) = (Co_de(i,1)+Co_de(i,4))*0.5
          Co_de_m(i,2) = (Co_de(i,2)+Co_de(i,3))*0.5
       end do
c Calculate the directional cosine & the transformation matrix
       d_x = Co_de_m(1,2) - Co_de_m(1,1)
       d_y = Co_de_m(2,2) - Co_de_m(2,1)
       el_length = (d_x**2 + d_y**2)**0.5
       cos_a = d_x / el_length
       sin_a = d_y / el_length
       R(1,1) = cos_a
       R(1,2) = sin_a
       R(2,1) = -sin_a
       R(2,2) = cos_a
       RETURN
       END
c
c =====================================================================
c = Matrix operations =================================================
       SUBROUTINE k_Matrix_Zero (A,n,m)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION A(n,m)
       do i = 1, n
          do j = 1, m
             A(i,j) = 0.0
          end do
       end do
       RETURN
       END
c
       SUBROUTINE k_Matrix_Transpose (A,B,n,m)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION A(n,m), B(m,n)
       call k_Matrix_zero (B,m,n)
       do i = 1, n
          do j = 1, m
             B(j,i) = A(i,j)
          end do
       end do
       RETURN
       END
c
       SUBROUTINE k_Matrix_PlusScalar (A,B,c,n,m)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION A(n,m), B(n,m)
       do i = 1, n
          do j = 1, m
             A(i,j) = A(i,j) + c*B(i,j)
          end do
       end do
       RETURN
       END
c
       SUBROUTINE k_Matrix_Multiply (A,B,C,l,n,m)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION A(l,n), B(n,m), C(l,m)
       call k_Matrix_zero (C,l,m)
       do i = 1, l
          do j = 1, m
             do k = 1, n
               C(i,j) = C(i,j) + A(i,k) * B (k,j)
             end do
          end do
       end do
       RETURN
       END
c =====================================================================
