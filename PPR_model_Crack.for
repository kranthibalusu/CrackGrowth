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
c Kranthi Balusu 
c Arizona State 
c     V.1.1 : sign max approximated , Cs and Ct have to be zero 
c     V.1.2 : failSt variable, to set intial crack 
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
c Basic hysteresis added from [1] J. Pre-proofs, Mixed-Mode Fatigue 
c     Crack Growth using Cohesive Zone Modeling, Eng. Fract. Mech. (2020)
c     107234. https://doi.org/10.1016/j.engfracmech.2020.107234.
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
c Variables updated in the UEL subroutinet
c   SVARS(1) : deln_max
c   SVARS(2) : delt_max
c   SVARS(3) : deln_pr : Separation previous 
c   SVARS(4) : delt_pr : speartion previous 
c   SVARS(5) : Tn_deldot : Traction rate 
c   SVARS(5) : Tt_deldot : Traction rate 
c   The same variables for the 2nd integration point 
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
c      PROPS(10): chi_n_Del = ! speration rate resistance paramaters
c      PROPS(11): chi_t_Del
c      PROPS(12): chi_n_T ! Stress rate resistance paramaters 
c      PROPS(13): chi_t_T 
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
c     ,delnb, deltb
      double precision Tn_pr, Tt_pr, Tn_deldot, Tt_deldot, deln_deldot, 
     & delt_deldot  
c     fatigue model Parameters 
      double precision chi_n_Del,chi_t_Del, chi_n_T, chi_t_T, C_S, C_T,
     & damn_S, damt_S, damn_T, damt_T, deltime, failSt1, failSt
      integer idx
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
c   deln_max, delt_max: Maximum separations in a loading history; or 
c       assuming the same as the  delnb, deltb :the boundary separation
c       (softening - unloading/reloading region) 
c   el_length: Length of a cohesive element
c   deln_del,delt_del: The change of in the previous timestep 
c   chi_n_Del, chi_t_Del :  seperation rate resistance paramaters
c   chi_n_T, chi_t_T :  ! Stress rate resistance paramaters
c   damn_S, damt_S : separation based damage paramater 
c   damn_T, damt_T : Traction based damage paramater 
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
       chi_n_Del = PROPS(10) ! speration rate resistance paramaters
       chi_t_Del = PROPS(11)
       chi_n_T = PROPS(12) ! Stress rate resistance paramaters 
       chi_t_T = PROPS(13)
       failSt1 = PROPS(14) ! failSt, if EQ 1; failed ; EQ 0 not failed
       C_S = 0.0 !PROPS(14) !contact speration ratio;set to zero for now
       C_T =  0.0!PROPS(15)!Traction speration ratio;set to zero for now
       
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
c       write(*,*) 'dn',dn
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
       do 11 i = 0, nnode-1
          U_l(1+i*mcrd) = R(1,1)*U(1+i*mcrd) + R(1,2)*U(2+i*mcrd)
          U_l(2+i*mcrd) = R(2,1)*U(1+i*mcrd) + R(2,2)*U(2+i*mcrd)
   11  continue 
       del1 = U_l(7) - U_l(1)
       del2 = U_l(8) - U_l(2)
       del3 = U_l(5) - U_l(3)
       del4 = U_l(6) - U_l(4)
c Numerical integration to compute RHS and AMATRX
       do 12 i = 1, n_GP

          N1 = 0.5*(1 - GP(i))
          N2 = 0.5*(1 + GP(i))
          del(1) = N1*del1 + N2*del3
          del(2) = N1*del2 + N2*del4

c         get varaibles from previous timesteps           
          idx = (i-1)*15

          delt_max = SVARS( idx+1)
          deln_max = SVARS(idx+2)
          deln_pr = SVARS(idx+3)
          delt_pr = SVARS(idx+4)
          Tn_pr = SVARS(idx+5)
          Tt_pr = SVARS(idx+6)
          Tn_deldot = SVARS(idx+7)
          Tt_deldot = SVARS(idx+8)
          damn_S =  SVARS(idx+9)
          damt_S =  SVARS(idx+10)
          damn_T =  SVARS(idx+11)
          damt_T =  SVARS(idx+12)
          damn_Tdot =  SVARS(idx+13)
          damt_Tdot =  SVARS(idx+14)
          failSt =  SVARS(idx+15)



          deln_deldot = (del(2) - deln_pr) !/DTIME
          delt_deldot = (del(1) - delt_pr)!/DTIME


          if ( TIME(1) .LT. 0.000000001 ) then
             delt_max = 0.0D0
             deln_max = 0.0D0
             deln_pr = 0.0D0
             delt_pr = 0.0D0
             Tn_pr = 0.0D0
             Tt_pr = 0.0D0
             Tn_deldot = 0.0D0
             Tt_deldot = 0.0D0
             damn_S =  0.0D0
             damt_S =  0.0D0
             damn_T =  0.0D0
             damt_T =  0.0D0
             damn_Tdot =  0.0D0
             damt_Tdot =  0.0D0

             deln_deldot = 0.0D0
             delt_deldot = 0.0D0
             failSt = 0
          endif  


c          write(*,*), 'damn_S', damn_S , TIME(1),i
c          write(*,*), 'damt_S', damt_S, TIME(1),i
c          write(*,*), 'damn_T', damn_T , TIME(1),i
c          write(*,*), 'damt_T', damt_T , TIME(1),i
c          write(*,*), 'damn_Tdot', damn_Tdot , TIME(1),i
c          write(*,*), 'damt_Tdot', damt_Tdot, TIME(1),i
c          write(*,*), 'deln_del', (del(2) - deln_pr) , TIME(1),i
c          write(*,*), 'delt_del', (del(1) - delt_pr), TIME(1),i 
c          write(*,*), 'Tn_pr', Tn_pr , TIME(1),i
c          write(*,*), 'Tt_pr', Tt_pr, TIME(1),i          
c          write(*,*), 'Tn_deldot', Tn_deldot , TIME(1),i
c          write(*,*), 'Tt_deldot', Tt_deldot, TIME(1),i    
cccccccccccccccc
C          write(*,*) 'delt_max',delt_max,del(1)
C          write(*,*) 'deln_max',deln_max,del(2)
          deltime = 1.0D0!DTIME ! just to be cautious 
          call k_Cohesive_PPR (T, T_d, Gam_n, Gam_t, alph, beta, m, n,
     &  dn, dt, dGtn, dGnt, del, deln_max, delt_max, deln_deldot,
     & delt_deldot, damn_S, damt_S, damn_T, damt_T, damn_Tdot, 
     & damt_Tdot, deltime, chi_n_Del, chi_t_Del, chi_n_T, chi_t_T,
     & Tn_deldot, Tt_deldot, failSt)

c          if ( TIME(2) .GT. 0.99) then
c             write(*,*), 'JELEM', JELEM,  TIME(2)
c             write(*,*),'damn_S',damn_S,'damt_S',damt_S
c             write(*,*),'damn_T', damn_T, 'damt_T', damt_T  
c          endif 

ccc    Element related; no need to touch this section           
          ShapeN(1) = -N1
          ShapeN(2) = -N2
          ShapeN(3) = N2
          ShapeN(4) = N1
          do 13 j = 1, nnode
             do 14 k = 1, mcrd
               do 15 l = 1, mcrd
                  Bc(k,l+(j-1)*mcrd) = ShapeN(j)*R(k,l)
   15          continue
   14        continue
   13     continue
          call k_Matrix_Transpose (Bc,Bct,mcrd,ndofel)
          call k_Matrix_Multiply (Bct,T_d,tmp,ndofel,mcrd,mcrd)
          call k_Matrix_Multiply (tmp,Bc,Sc,ndofel,mcrd,ndofel)
          call k_Matrix_Multiply (Bct,T,Fc,ndofel,mcrd,nrhs)
          thick = 0.5 * el_length * GP_w(i) * th
          call k_Matrix_PlusScalar (AMATRX,Sc,thick,ndofel,ndofel)
          call k_Matrix_PlusScalar (RHS,-Fc,thick,ndofel,nrhs)

c   Update the state variables: SVARS
          if((delt_max.LT.abs(del(1))).AND.(abs(del(1)).GT.lt*dt)) then
             SVARS(idx+1) = abs(del(1))
          end if
          if ((deln_max .LT. del(2)) .AND. (del(2) .GT. ln*dn)) then
             SVARS(idx+2) = del(2)
          end if


          SVARS(idx+4) =abs(del(1)) ! delt_pr
          SVARS(idx+3) = del(2) ! deln_pr 

          SVARS(idx+5) = T(2,1) !Tn_pr
          SVARS(idx+6) = T(1,1) ! Tt_pr
c          write(*,*), 'Tn_pr', T(2,1) , TIME(1),i
c          write(*,*), 'Tt_pr', T(1,1) , TIME(1),i       
 
          SVARS(idx+7) = (T(2,1) - Tn_pr) !/DTIME ! Tn _deldot 
          SVARS(idx+8) = (T(1,1) - Tt_pr) !/DTIME ! Tt _deldot
          
          SVARS(idx+9) = damn_S
          SVARS(idx+10) = damt_S  
          SVARS(idx+11) = damn_T
          SVARS(idx+12) = damt_T  
          SVARS(idx+13) = damn_Tdot
          SVARS(idx+14) = damt_Tdot
          SVARS(idx+15) = failSt



   12  continue

       RETURN
       END
c
c ====================================================================
c = Cohesive traction-separation relation of the PPR model ===========
c       SUBROUTINE k_Cohesive_PPR (T,T_d,Gam_n,Gam_t,alph,beta,m,n,
c     &  dn, dt, dGtn, dGnt, del, deln_max, delt_max, deln_deldot,
c     & delt_deldot)
      SUBROUTINE k_Cohesive_PPR (T, T_d, Gam_n, Gam_t, alph, beta, m,n,
     &  dn, dt, dGtn, dGnt, del, deln_max, delt_max, deln_deldot,
     & delt_deldot, damn_S, damt_S, damn_T, damt_T, damn_Tdot, 
     & damt_Tdot, deltime, chi_n_Del, chi_t_Del, chi_n_T, chi_t_T, 
     & Tn_deldot, Tt_deldot, failSt)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION T(2,1), T_d(2,2), del(2)
       DOUBLE PRECISION Gam_n, Gam_t, alph, beta, m, n, dn, dt,
     & dGtn, dGnt, deln_max, delt_max, Tn, Tt, deln, delt, sign_dt
     & ,deln_deldot,delt_deldot
       double precision  damn_S, damt_S, damn_T, damt_T, damn_Tdot, 
     & damt_Tdot, deltime
c     local dc  
       double precision damn_Sdot, damt_Sdot, chi_n_Del, chi_t_Del, 
     &  sign_max, sigt_max, Tn_deldot, Tt_deldot, C_T, C_S, failSt
c     chi_n_T, chi_t_T,
c     Run a monotonics case to determine these first
      sign_max = 529e3
      sigt_max = 529e3
ccccccccc
      C_T = 0.0D0
      C_S = 0.0D0

       delt = abs(del(1))
       deln = del(2)
       if (del(1) .GE. 0) then
          sign_dt = 1
       else
          sign_dt = -1
       end if
       Tn = 0

cc     damage stuff 

c     separation damage evolution 







CCCCCCCCCCCCC     Traction calculation        
c Equation (12):
c Pre-calculation of the normal cohesive traction, Tn
       if (deln .LT. 0) then
          deln = 0
c          write(*,*) 'I am in contact' 
       elseif ((deln .GE. 0.98*dn) .OR. (delt .GE. 0.98*dt)) then 
c     & .OR. (failSt .EQ. 1)) then 
          if (failSt .LT. 0.1) then
             failSt = 1
             write(*,*) 'Element ', JELEM, 'failed at time', TIME(2) 
          endif 
          Tn = 0
       elseif (deln .GE. deln_max) then ! softening 
c          write(*,*) 'I am in softening' 
cccccc
          damn_Sdot = deln_deldot/dn
          damn_S = damn_S + damn_Sdot*deltime

          damn_Tdot = 0.0
          damn_T = damn_T+damn_Tdot*deltime

          Tn = ((Gam_t*(1-delt/dt)**beta*(delt/dt+n/beta)**n+dGtn) *
     & Gam_n/dn*(m*(1-deln/dn)**alph*(m/alph+deln/dn)**(m-1)
     & -alph*(1-deln/dn)**(alph-1)*(m/alph+deln/dn)**m))*(1-damn_T)


ccccc
       elseif (deln_deldot .LT. 0.0) then ! unloading 

          damn_Sdot = 0.0
          damn_S = damn_S+damn_Sdot*deltime
c          write(*,*) 'I am in unloading' 

          damn_Tdot = 0.0
          damn_T = damn_T+damn_Tdot*deltime

          Tn = ((Gam_t*(1-delt/dt)**beta*(delt/dt+n/beta)**n+dGtn) *
     & Gam_n/dn*(m*(1-deln_max/dn)**alph*(m/alph+deln_max/dn)**(m-1)
     &   -alph*(1-deln_max/dn)**(alph-1)*(m/alph+deln_max/dn)**m) *
     & ( C_T + (deln/deln_max - C_S)*(1-C_T)/(1-C_S)) )*(1-damn_T)

ccccc
       else !reloading 
c          write(*,*) 'I am in reloading ' 

          damn_Sdot = deln_deldot/dn/chi_n_Del
          damn_S = damn_S+damn_Sdot*deltime
c          write(*,*) 'deln_deldot ' ,deln_deldot
c          write(*,*) 'damn_Sdot ' ,damn_Sdot
c          write(*,*) 'deln_max ' ,deln_max
c          deln_max = damn_S*dn
          
c          write(*,*) 'deln_max ' ,deln_max
c          write(*,*) 'Tn_deldot ' ,Tn_deldot
          damn_Tdot = Tn_deldot/sign_max/chi_n_T
c          write(*,*) 'damn_T ' ,damn_T
          damn_T = damn_T+damn_Tdot*deltime
c          write(*,*) 'damn_T ' ,damn_T
          


          Tn = ((Gam_t*(1-delt/dt)**beta*(delt/dt+n/beta)**n+dGtn) *
     & Gam_n/dn*(m*(1-deln_max/dn)**alph*(m/alph+deln_max/dn)**(m-1)
     &   -alph*(1-deln_max/dn)**(alph-1)*(m/alph+deln_max/dn)**m) *
     & ( C_T + (deln/deln_max - C_S)*(1-C_T)/(1-C_S)) )*(1-damn_T)
c          write(*,*) 'Tn ' ,Tn


       end if

c Pre-calculation of the tangential cohesive traction, Tt
       if ((deln .GE. dn) .OR. (delt .GE. dt) 
     & .OR. (failSt .EQ. 1)) then
          Tt = 0
       elseif (delt .GE. delt_max) then ! softening 
cccc   
          damt_Sdot = delt_deldot/dt
          damt_S = damt_S+damt_Sdot*deltime

          damt_Tdot = 0.0
          damt_T = damt_T + damt_Tdot*deltime

          Tt = ((Gam_n*(1-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt*(n*(1-delt/dt)**beta*(delt/dt+n/beta)**(n-1)
     & -beta*(1-delt/dt)**(beta-1)*(delt/dt+n/beta)**n))*(1-damt_T)

       elseif (delt_deldot .LT. 0.0) then ! unloading ! 
cccc    

          damt_Sdot = 0.0
          damt_S = damt_S+damt_Sdot*deltime

          damt_Tdot = 0.0
          damt_T = damt_T+damt_Tdot*deltime

          Tt = ((Gam_n*(1-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt*(n*(1-delt_max/dt)**beta*(delt_max/dt+n/beta)**(n-1)
     &   -beta*(1-delt_max/dt)**(beta-1)*(delt_max/dt+n/beta)**n) *
     & delt/delt_max)*(1-damt_T)

cccc
       else ! reloading
cccc        
          damt_Sdot = delt_deldot/dt/chi_t_Del
          damt_S = damt_S+damt_Sdot*deltime
c          delt_max = damt_S*dt          

          damt_Tdot = Tt_deldot/sigt_max/chi_t_T
          damt_T = damt_T+damt_Tdot*deltime

          Tt = ((Gam_n*(1-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt*(n*(1-delt_max/dt)**beta*(delt_max/dt+n/beta)**(n-1)
     &   -beta*(1-delt_max/dt)**(beta-1)*(delt_max/dt+n/beta)**n) *
     & delt/delt_max)*(1-damt_T)
       end if

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC derivative of T        
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
     & ((Gam_t*(1-delt/dt)**beta*(delt/dt+n/beta)**n+dGtn) *
     & Gam_n/dn**2 *
     & ((1-deln/dn)**(alph-2)*(alph-1)*alph*(deln/dn+m/alph)**m -
     & 2*(1-deln/dn)**(alph-1)*alph*(deln/dn+m/alph)**(m-1)*m +
     & (1-deln/dn)**alph*(deln/dn+m/alph)**(m-2)*(m-1)*m))*(1-damn_T)

             T_d(2,1) =
     & (Gam_t/dt*(-(1-delt/dt)**(beta-1)*beta*(delt/dt+n/beta)**n +
     & (1-delt/dt)**beta*(delt/dt+n/beta)**(n-1)*n) * sign_dt *
     & Gam_n/dn*(-(1-deln/dn)**(alph-1)*alph*(deln/dn+m/alph)**m +
     & (1-deln/dn)**alph*(deln/dn+m/alph)**(m-1)*m))*(1-damn_T)

c   (3) Unloading/reloading condition
          elseif (deln_deldot .LT. 0.0) then ! unloading 
             T_d(2,2) =
     & ((Gam_t*(1-delt/dt)**beta*(delt/dt+n/beta)**n+dGtn) *
     & Gam_n/dn*((1-deln_max/dn)**alph*(deln_max/dn+m/alph)**(m-1)*m
     & -(1-deln_max/dn)**(alph-1)*alph*(deln_max/dn+m/alph)**m)
     & / deln_max)*(1-damn_T)
             T_d(2,1) =
     & (Gam_t/dt*(-(1-delt/dt)**(beta-1)*beta*(delt/dt+n/beta)**n +
     & (1-delt/dt)**beta*(delt/dt+n/beta)**(n-1)*n) * sign_dt *
     & Gam_n/dn*(m*(1-deln_max/dn)**alph*(m/alph+deln_max/dn)**(m-1)
     &   -alph*(1-deln_max/dn)**(alph-1)*(m/alph+deln_max/dn)**m) *
     & deln/deln_max)*(1-damn_T)
          else ! reloading  
             T_d(2,2) =
     & (((Gam_t*(1-delt/dt)**beta*(delt/dt+n/beta)**n+dGtn) *
     & Gam_n/dn*((1-deln_max/dn)**alph*(deln_max/dn+m/alph)**(m-1)*m
     & -(1-deln_max/dn)**(alph-1)*alph*(deln_max/dn+m/alph)**m)
     & / deln_max))*(1-damn_T)
             T_d(2,1) =
     & ((Gam_t/dt*(-(1-delt/dt)**(beta-1)*beta*(delt/dt+n/beta)**n +
     & (1-delt/dt)**beta*(delt/dt+n/beta)**(n-1)*n) * sign_dt *
     & Gam_n/dn*(m*(1-deln_max/dn)**alph*(m/alph+deln_max/dn)**(m-1)
     &   -alph*(1-deln_max/dn)**(alph-1)*(m/alph+deln_max/dn)**m) *
     & deln/deln_max))*(1-damn_T)

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
     & ((Gam_n*(1-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt**2 *
     & ((1-delt/dt)**(beta-2)*(beta-1)*beta*(delt/dt+n/beta)**n -
     & 2*(1-delt/dt)**(beta-1)*beta*(delt/dt+n/beta)**(n-1)*n +
     & (1-delt/dt)**beta*(delt/dt+n/beta)**(n-2)*(n-1)*n))*(1-damt_T)
             T_d(1,2) =
     & (Gam_t/dt*(-(1-delt/dt)**(beta-1)*beta*(delt/dt+n/beta)**n +
     & (1-delt/dt)**beta*(delt/dt+n/beta)**(n-1)*n) * sign_dt *
     & Gam_n/dn*(-(1-deln/dn)**(alph-1)*alph*(deln/dn+m/alph)**m +
     & (1-deln/dn)**alph*(deln/dn+m/alph)**(m-1)*m))*(1-damt_T)
c   (2) Unloading/reloading condition
          elseif (deln_deldot .LT. 0.0) then ! unloading ! 
             T_d(1,1) =
     & ((Gam_n*(1-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt*(n*(1-delt_max/dt)**beta*(delt_max/dt+n/beta)**(n-1)
     &      -beta*(1-delt_max/dt)**(beta-1)*(delt_max/dt+n/beta)**n)
     & / delt_max)*(1-damt_T)
             T_d(1,2) =
     & (Gam_n/dn*(-(1-deln/dn)**(alph-1)*alph*(deln/dn+m/alph)**m +
     & (1-deln/dn)**alph*(deln/dn+m/alph)**(m-1)*m) * sign_dt *
     & Gam_t/dt*(n*(1-delt_max/dt)**beta*(delt_max/dt+n/beta)**(n-1)
     &   -beta*(1-delt_max/dt)**(beta-1)*(delt_max/dt+n/beta)**n) *
     & delt/delt_max)*(1-damt_T)
          else !relaoding 
             T_d(1,1) =
     & (((Gam_n*(1-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt*(n*(1-delt_max/dt)**beta*(delt_max/dt+n/beta)**(n-1)
     &      -beta*(1-delt_max/dt)**(beta-1)*(delt_max/dt+n/beta)**n)
     & / delt_max))*(1-damt_T)
             T_d(1,2) =
     & ((Gam_n/dn*(-(1-deln/dn)**(alph-1)*alph*(deln/dn+m/alph)**m +
     & (1-deln/dn)**alph*(deln/dn+m/alph)**(m-1)*m) * sign_dt *
     & Gam_t/dt*(n*(1-delt_max/dt)**beta*(delt_max/dt+n/beta)**(n-1)
     &   -beta*(1-delt_max/dt)**(beta-1)*(delt_max/dt+n/beta)**n) *
     & delt/delt_max))*(1-damt_T)
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
       do 16 i = 1, mcrd
          do 17 j = 1, nnode
             Co_de(i,j) = COORDS(i,j) + U(2*(j-1)+i)
   17     continue
   16  continue
       do 18 i = 1, 2
          Co_de_m(i,1) = (Co_de(i,1)+Co_de(i,4))*0.5
          Co_de_m(i,2) = (Co_de(i,2)+Co_de(i,3))*0.5
   18    continue
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
