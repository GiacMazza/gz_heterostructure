MODULE EQS_OF_MOTION
  USE GLOBAL

CONTAINS


  function slab_lead_eom(time,y,Nsys) result(f)
    implicit none
    !inputs
    integer                                     :: Nsys ! nr of equations
    real(8)                                     :: time ! time variable
    complex(8),dimension(Nsys)                  :: y    ! argument array
    complex(8),dimension(Nsys)                  :: f    ! result


    !indeces 
    integer                                     :: ik,ik_sys,ik_sys0
    integer                                     :: dimk,dimLead,dL,dim_hyb
    integer                                     :: ilead,jlead,ihyb

    integer                                     :: side_i,side_j
    integer                                     :: ihyb_i,ihyb_j

    integer                                     :: isum_hyb_R
    integer                                     :: isum_hyb_L

    integer                                     :: iSlab,jSlab,kSlab
    integer                                     :: igz

    integer                                     :: j_time
    integer                                     :: lx,ly

    !auxiliary GF
    complex(8),dimension(:,:),allocatable       :: Glead
    complex(8),dimension(:,:),allocatable       :: Gslab


    !k-sums
    complex(8),dimension(:),allocatable         :: sumGleadR,sumGleadL
    complex(8),dimension(:),allocatable         :: sumGhybL,sumGhybR
    complex(8),dimension(:),allocatable         :: sumGhybL_k,sumGhybR_k

    complex(8),dimension(:),allocatable         :: hop_plus,hop_minus,hop
    complex(8)                                  :: hyb_left,hyb_right

    !gutzwiller
    type(gz_projector),dimension(:),allocatable :: gz_phi
    complex(8),dimension(:),allocatable         :: gz_R
    real(8),dimension(:),allocatable            :: x
    complex(8)                                  :: h00,h01,h22
    complex(8)                                  :: prj_s0,prj_s1,prj_s2
    real(8)                                     :: P0,P1,P2

    !dispersions and momenta
    real(8)                                     :: ek
    real(8)                                     :: kl_i,kl_j,ekl_i,ekl_j
    real(8)                                     :: vk_i,vk_j,chem_lead
    !chemical potential
    real(8)                                     :: chem_i,chem_j
    !coulom potential
    real(8),dimension(L)                        :: coulomb
    real(8) :: tnk

    dimk    = 2*Nk_orth*L + L*L
    dimLead = 2*Nk_orth                       
    dL      = dimLead**2
    dim_hyb = dimLead*L

    j_time = 2*(time+1.d-5)/dt + 1

    allocate(hop_plus(L),hop_minus(L),hop(L))
    !+-------------------------+!
    !+- gutzwiller quantities -+!
    !+-------------------------+!
    allocate(gz_phi(L),gz_R(L),x(L))
    do iSlab = 1,L
       ik_sys = Nk_tot*dimk + 3*(iSlab-1) + 2*Nk_orth*2*Nk_orth
       gz_phi(iSlab)%p0 = y(ik_sys+1) 
       gz_phi(iSlab)%p1 = y(ik_sys+2) 
       gz_phi(iSlab)%p2 = y(ik_sys+3) 
    end do
    gz_R = GZ_hop(gz_phi)
    x = GZ_doping(gz_phi)


    !initialize k sums
    hop_plus  = Z0
    hop_minus = Z0
    hop       = Z0
    hyb_left  = Z0
    hyb_right = Z0


    !+- compute k_sums over hybrid green functions -+!
    allocate(sumGhybL_k(dimLead),sumGhybR_k(dimLead))
    do ilead=1,dimLead
       sumGhybL_k(ilead) = Z0
       sumGhybR_k(ilead) = Z0
       do ik = 1,Nk_tot
          isum_hyb_L = dimLead*dimLead + (ik-1)*dimk + ilead
          isum_hyb_R = dimLead*dimLead + (ik-1)*dimk + (L-1)*dimLead + ilead
          sumGhybL_k(ilead) = sumGhybL_k(ilead) + y(isum_hyb_L)*wt(ik)
          sumGhybR_k(ilead) = sumGhybR_k(ilead) + y(isum_hyb_R)*wt(ik)
       end do
    end do



    allocate(sumGleadR(dimLead),sumGleadL(dimLead),sumGhybL(L),sumGhybR(L))
    allocate(Gslab(L,L))

    !+- initializae global index -+!
    ik_sys = 0

    !+- LEADS -+!

    !+- left-(left/right) -+!
    do ilead=1,Nk_orth
       side_i = 1
       kl_i = k_orth(ilead)
       chem_lead = muL(j_time)
       vk_i = vk_L(ilead,j_time)
       ekl_i = chain_disp(kl_i)*t_lead - chem_lead


       sumGleadL(ilead) = 0.d0
       sumGleadR(ilead) = 0.d0

       !+- LEFT-LEFT-+!
       do jlead = 1,Nk_orth
          !+- update global index -+!
          ik_sys = ik_sys + 1
          !+-----------------------+!
          side_j = 1
          kl_j = k_orth(jlead)
          chem_lead = muL(j_time)
          vk_j = vk_L(jlead,j_time)
          ekl_j = chain_disp(kl_j)*t_lead - chem_lead

          !+- update Glead sums -+!
          sumGleadL(ilead) = sumGleadL(ilead) + vk_j*y(ik_sys)
          !+---------------------+!

          f(ik_sys) = -Zi*(ekl_i-ekl_j)*y(ik_sys) + &
               Zi*vk_j*conjg(gz_R(side_j))*sumGhybL_k(ilead) + &
               Zi*vk_i*gz_R(side_i)*conjg(sumGhybL_k(jlead))
       end do

       !+- LEFT/RIGHT -+!
       do jlead = Nk_orth+1,2*Nk_orth
          !+- update global index -+!
          ik_sys = ik_sys + 1
          !+-----------------------+!
          kl_j = k_orth(jlead-Nk_orth)
          chem_lead = muR(j_time)
          vk_j = vk_R(jlead-Nk_orth,j_time)

          side_j = L

          ekl_j = chain_disp(kl_j)*t_lead - chem_lead

          sumGleadR(ilead) = sumGleadR(ilead) + vk_j*y(ik_sys)   

          f(ik_sys) = -Zi*(ekl_i-ekl_j)*y(ik_sys) + &
               Zi*vk_j*conjg(gz_R(side_j))*sumGhybR_k(ilead) + &
               Zi*vk_i*gz_R(side_i)*conjg(sumGhybL_k(jlead))
       end do
    end do


    !+- right/(left-right) -+!
    do ilead=Nk_orth+1,2*Nk_orth
       kl_i = k_orth(ilead-Nk_orth)
       chem_lead = muR(j_time)
       vk_i = vk_R(ilead-Nk_orth,j_time)
       side_i = L
       ekl_i = chain_disp(kl_i)*t_lead - chem_lead

       sumGleadL(ilead) = 0.d0
       sumGleadR(ilead) = 0.d0

       !+- RIGHT-LEFT -+!
       do jlead = 1,Nk_orth
          !+- update global index -+!
          ik_sys = ik_sys + 1
          !+-----------------------+!
          kl_j = k_orth(jlead)
          chem_lead = muL(j_time)
          vk_j = vk_L(jlead,j_time)
          sumGleadL(ilead) = sumGleadL(ilead) + vk_j*y(ik_sys)
          side_j = 1

          ekl_j = chain_disp(kl_j)*t_lead - chem_lead

          f(ik_sys) = -Zi*(ekl_i-ekl_j)*y(ik_sys) + &
               Zi*vk_j*conjg(gz_R(side_j))*sumGhybL_k(ilead) + &
               Zi*vk_i*gz_R(side_i)*conjg(sumGhybR_k(jlead))
       end do

       !+- RIGHT-RIGHT -+!
       do jlead = Nk_orth+1,2*Nk_orth
          !+- update global index -+!
          ik_sys = ik_sys + 1
          !+-----------------------+!
          kl_j = k_orth(jlead-Nk_orth)
          chem_lead = muR(j_time)
          vk_j = vk_R(jlead-Nk_orth,j_time)
          sumGleadR(ilead) = sumGleadR(ilead) + vk_j*y(ik_sys)
          side_j = L

          ekl_j = chain_disp(kl_j)*t_lead - chem_lead

          f(ik_sys) = -Zi*(ekl_i-ekl_j)*y(ik_sys) + &
               Zi*vk_j*conjg(gz_R(side_j))*sumGhybR_k(ilead) + &
               Zi*vk_i*gz_R(side_i)*conjg(sumGhybR_k(jlead))
       end do
    end do

    !+- LOOP ON (planar) k POINTS -+!
    do ik=1,Nk_tot
       ek = square_lattice_disp(vec_k(ik))
       ik_sys0 = ik_sys

       !+- hybrid green -+!
       do iSlab = 1,L
          sumGhybL(iSlab) = Z0
          sumGhybR(iSlab) = Z0

          !+- LEFT HYBRID -+!
          do ihyb=1,Nk_orth
             !+- update global index -+!
             ik_sys = ik_sys+1
             !+-----------------------+!
             kl_i = k_orth(ihyb)
             chem_lead = muL(j_time)
             vk_i = vk_L(ihyb,j_time)
             side_i = 1
             sumGhybL(iSlab) = sumGhybL(iSlab) + vk_i*y(ik_sys)

             ekl_i = chain_disp(kl_i)*t_lead - chem_lead

             !+--------------------------------+!
             !+- HYBRID GF EQUATION OF MOTION -+!
             !+--------------------------------+!
             f(ik_sys) = -Zi*ekl_i*y(ik_sys) + Zi*ek*ABS(gz_R(iSlab))**2*y(ik_sys) 

             if(iSlab.gt.1) then
                f(ik_sys) = f(ik_sys)-Zi*y(ik_sys-dimLead)*gz_R(iSlab)*conjg(gz_R(iSlab-1)) 
             end if

             if(iSlab.lt.L) then
                f(ik_sys) = f(ik_sys)-Zi*y(ik_sys+dimLead)*gz_R(iSlab)*conjg(gz_R(iSlab+1)) 
             end if

             ihyb_i = ik_sys0 + dimk - L*L +(side_i-1)*L + iSlab ! it's ok upon changing dimk=dimlead*L + L*L
             f(ik_sys) = f(ik_sys)-Zi*vk_i*y(ihyb_i)*gz_R(side_i)

             if(iSlab.eq.1) then
                f(ik_sys) = f(ik_sys) + Zi*sumGleadL(ihyb)*gz_R(iSlab)
             end if

             if(iSlab.eq.L) then
                f(ik_sys) = f(ik_sys) + Zi*sumGleadR(ihyb)*gz_R(iSlab)
             end if

          end do


          !+- RIGHT HYBRID -+!
          do ihyb=Nk_orth+1,2*Nk_orth
             !+- update global index -+!
             ik_sys = ik_sys+1
             !+-----------------------+!
             kl_i = k_orth(ihyb-Nk_orth)
             chem_lead = muR(j_time)
             vk_i = vk_R(ihyb-Nk_orth,j_time)
             side_i = L
             sumGhybR(iSlab) = sumGhybR(iSlab) + vk_i*y(ik_sys)
             ekl_i = chain_disp(kl_i)*t_lead - chem_lead
             !+--------------------------------+!
             !+- HYBRID GF EQUATION OF MOTION -+!
             !+--------------------------------+!
             f(ik_sys) = -Zi*ekl_i*y(ik_sys) + Zi*ek*ABS(gz_R(iSlab))**2*y(ik_sys) 

             if(iSlab.gt.1) then
                f(ik_sys) = f(ik_sys)-Zi*y(ik_sys-dimLead)*gz_R(iSlab)*conjg(gz_R(iSlab-1)) 
             end if

             if(iSlab.lt.L) then
                f(ik_sys) = f(ik_sys)-Zi*y(ik_sys+dimLead)*gz_R(iSlab)*conjg(gz_R(iSlab+1)) 
             end if

             ihyb_i = ik_sys0 + dimk - L*L +(side_i-1)*L + iSlab !it's ok upon changing dimk

             f(ik_sys) = f(ik_sys)-Zi*vk_i*y(ihyb_i)*gz_R(side_i)

             if(iSlab.eq.1) then
                f(ik_sys) = f(ik_sys) + Zi*sumGleadL(ihyb)*gz_R(iSlab)
             end if

             if(iSlab.eq.L) then
                f(ik_sys) = f(ik_sys) + Zi*sumGleadR(ihyb)*gz_R(iSlab)
             end if
          end do


       end do

       !+- dot GSlab -+!       

       do iSlab=1,L
          do jSlab=1,L
             ik_sys = ik_sys+1
             Gslab(iSlab,jSlab) = y(ik_sys)
          end do
       end do
       !+- put back global index -+!
       ik_sys = ik_sys - L*L
       !+-------------------------+!
       do iSlab=1,L
          do jSlab=1,L
             !+- up date global index -+!
             ik_sys = ik_sys+1
             !+-------------------------------+!
             !+- SLAB GF EQUATIONS OF MOTION -+!
             !+-------------------------------+!
             f(ik_sys) = -Zi*ek*( ABS(gz_R(iSlab))**2 - ABS(gz_R(jSlab))**2 )*GSlab(iSlab,jSlab)    
             if(iSlab.gt.1) then
                f(ik_sys) = f(ik_sys) + Zi*conjg(gz_R(iSlab))*gz_R(iSlab-1)*Gslab(iSlab-1,jSlab)
             end if

             if(iSlab.lt.L) then
                f(ik_sys) = f(ik_sys) + Zi*conjg(gz_R(iSlab))*gz_R(iSlab+1)*Gslab(iSlab+1,jSlab)
             end if

             if(jSlab.gt.1) then
                f(ik_sys) = f(ik_sys) - Zi*conjg(gz_R(jSlab-1))*gz_R(jSlab)*Gslab(iSlab,jSlab-1)
             end if

             if(jSlab.lt.L) then
                f(ik_sys) = f(ik_sys) - Zi*conjg(gz_R(jSlab+1))*gz_R(jSlab)*Gslab(iSlab,jSlab+1)
             end if

             !+- hybridization -+!
             if(iSlab.eq.1) then
                f(ik_sys) = f(ik_sys)-Zi*(sumGhybL(jSlab))*conjg(gz_R(iSlab)) 
             end if

             if(iSlab.eq.L) then
                f(ik_sys) = f(ik_sys)-Zi*(sumGhybR(jSlab))*conjg(gz_R(iSlab)) 
             end if

             if(jSlab.eq.1) then
                f(ik_sys) = f(ik_sys)-Zi*conjg(sumGhybL(iSlab))*gz_R(jSlab) 
             end if

             if(jSlab.eq.L) then
                f(ik_sys) = f(ik_sys)-Zi*conjg(sumGhybR(iSlab))*gz_R(jSlab) 
             end if

          end do
       end do


       !+- accumulate k-point sums -+!
       do iSlab = 1,L
          if(iSlab.lt.L) then
             hop_plus(iSlab)  = hop_plus(iSlab)-2.d0*Zi*(Gslab(iSlab,iSlab+1))*wt(ik)
             hop_minus(iSlab) = hop_minus(iSlab)-2.d0*Zi*(Gslab(iSlab+1,iSlab))*wt(ik)
          end if
          hop(iSlab) = hop(iSlab) + 2.d0*(1.d0-Zi*Gslab(iSlab,iSlab))*wt(ik)*ek 
       end do

       hyb_left = hyb_left - 2.d0*Zi*sumGhybL(1)*wt(ik)
       hyb_right = hyb_right - 2.d0*Zi*sumGhybR(L)*wt(ik)

    end do !Nk_tot


    coulomb = 0.d0
    if(dabs(alpha_electrostatic).gt.1.d-5) then
       do lx = 1,L
          coulomb(lx) = 0.d0
          do ly = 1,L
             coulomb(lx) = coulomb(lx) + alpha_electrostatic*abs(dble(lx-ly))*x(ly)
          end do
       end do
    end if

    do iSlab = 1,L

       h00 = 2*ABS(gz_R(iSlab))**2*hop(iSlab)
       h01 = conjg(gz_R(iSlab))*hop(iSlab)

       if(iSlab.gt.1) then
          h00 = h00 - 2.d0*REAL(conjg(gz_R(iSlab-1))*gz_R(iSlab)*hop_minus(iSlab-1))
          h01 = h01 - conjg(gz_R(iSlab-1))*hop_minus(iSlab-1)
       else
          h00 = h00 + 2.d0*REAL(conjg(hyb_left)*gz_R(iSlab))
          h01 = h01 + conjg(hyb_left)
       end if

       if(iSlab.lt.L) then
          h00 = h00 - 2.d0*REAL(conjg(gz_R(iSlab+1))*gz_R(iSlab)*hop_plus(iSlab))
          h01 = h01 - conjg(gz_R(iSlab+1))*hop_plus(iSlab)
       else
          h00 = h00 + 2.d0*REAL(conjg(hyb_right)*gz_R(iSlab))
          h01 = h01 + conjg(hyb_right)
       end if

       h22 = x(iSlab)/(x(iSlab)**2-1.d0)*h00
       h00 = x(iSlab)/(1.d0-x(iSlab)**2)*h00
       h01 = sqrt(2.d0/(1.d0-x(iSlab)**2))*h01 


       h22 = h22 + U*0.5d0  + e_loc(iSlab,j_time) + coulomb(iSlab)
       h00 = h00 + U*0.5d0  - e_loc(iSlab,j_time) - coulomb(iSlab)

       !+- prj_entropy contribution -+!
       P0 = 0.25d0*(1.d0+x(iSlab))**2.d0
       P1 = 0.5d0*(1.d0-x(iSlab)**2.d0)
       P2 = 0.25d0*(1.d0-x(iSlab))**2.d0

       prj_s0 = gz_phi(iSlab)%p0
       prj_s1 = gz_phi(iSlab)%p1 
       prj_s2 = gz_phi(iSlab)%p2

       prj_s0 = prj_s0 - 0.5d0*( abs(gz_phi(iSlab)%p0)**2.d0/P0*(1.d0+x(iSlab)) - abs(gz_phi(iSlab)%p2)**2.d0/P2*(1.d0-x(iSlab)) )*gz_phi(iSlab)%p0
       prj_s0 = prj_s0 + x(iSlab)*abs(gz_phi(iSlab)%p1)**2.d0/P1*gz_phi(iSlab)%p0

       prj_s2 = prj_s2 + 0.5d0*( abs(gz_phi(iSlab)%p0)**2.d0/P0*(1.d0+x(iSlab)) - abs(gz_phi(iSlab)%p2)**2.d0/P2*(1.d0-x(iSlab)) )*gz_phi(iSlab)%p2
       prj_s2 = prj_s2 - x(iSlab)*abs(gz_phi(iSlab)%p1)**2.d0/P1*gz_phi(iSlab)%p2

       if( abs(gz_phi(iSlab)%p0) .gt. 1.d-10) then
          prj_s0 = prj_s0 + log( abs(gz_phi(iSlab)%p0)**2.d0 / P0 )*gz_phi(iSlab)%p0
       end if
       if( abs(gz_phi(iSlab)%p1) .gt. 1.d-10) then
          prj_s1 = prj_s1 + log( abs(gz_phi(iSlab)%p1)**2.d0 / P1 )*gz_phi(iSlab)%p1
       end if
       if( abs(gz_phi(iSlab)%p2) .gt. 1.d-10) then
          prj_s2 = prj_s2 + log( abs(gz_phi(iSlab)%p2)**2.d0 / P2 )*gz_phi(iSlab)%p2
       end if
       !+----------------------------+!

       do igz = 1,3
          !+- update global index -+!
          ik_sys = ik_sys + 1
          !+-----------------------+!
          select case(igz)
          case(1)
             f(ik_sys) = -Zi*(h00*gz_phi(iSlab)%p0 + h01*gz_phi(iSlab)%p1)        
             f(ik_sys) = f(ik_sys) - Zi*temp*prj_s0
          case(2)
             f(ik_sys) = -Zi*(conjg(h01)*gz_phi(iSlab)%p0 + h01*gz_phi(iSlab)%p2) 
             f(ik_sys) = f(ik_sys) - Zi*temp*prj_s1
          case(3)
             f(ik_sys) = -Zi*(h22*gz_phi(iSlab)%p2 + conjg(h01)*gz_phi(iSlab)%p1) 
             f(ik_sys) = f(ik_sys) - Zi*temp*prj_s2
          end select
       end do


    end do

    deallocate(gz_phi)
    deallocate(sumGleadL,sumGleadR)
    deallocate(sumGhybR,sumGhybL)
    deallocate(Gslab)

    return

  end function slab_lead_eom





  function slab_lead_eom_kleads(time,y,Nsys) result(f)
    implicit none
    !inputs
    integer                                     :: Nsys ! nr of equations
    real(8)                                     :: time ! time variable
    complex(8),dimension(Nsys)                  :: y    ! argument array
    complex(8),dimension(Nsys)                  :: f    ! result
    !indeces 
    integer                                     :: ik,ik_sys,ik_sys0
    integer                                     :: dimk,dimLead,dL
    integer                                     :: ilead,jlead,ihyb
    integer                                     :: side_i,side_j
    integer                                     :: ihyb_i,ihyb_j
    integer                                     :: iSlab,jSlab,kSlab
    integer                                     :: igz
    integer                                     :: j_time
    integer                                     :: lx,ly
    !auxiliary GF
    complex(8),dimension(:,:),allocatable       :: Glead
    complex(8),dimension(:,:),allocatable       :: Gslab
    !k-sums
    complex(8),dimension(:),allocatable         :: sumGleadR,sumGleadL
    complex(8),dimension(:),allocatable         :: sumGhybL,sumGhybR
    complex(8),dimension(:),allocatable         :: hop_plus,hop_minus,hop
    complex(8)                                  :: hyb_left,hyb_right
    !gutzwiller
    type(gz_projector),dimension(:),allocatable :: gz_phi
    complex(8),dimension(:),allocatable         :: gz_R
    real(8),dimension(:),allocatable            :: x
    complex(8)                                  :: h00,h01,h22
    complex(8)                                  :: prj_s0,prj_s1,prj_s2
    real(8)                                     :: P0,P1,P2
    !dispersions and momenta
    real(8)                                     :: ek
    real(8)                                     :: kl_i,kl_j,ekl_i,ekl_j
    real(8)                                     :: vk_i,vk_j,chem_lead
    !chemical potential
    real(8)                                     :: chem_i,chem_j
    !coulom potential
    real(8),dimension(L)                        :: coulomb

    dimk    = 2*Nk_orth*(2*Nk_orth+L) + L*L   !nr of GF for each k(2d) point
    dimLead = 2*Nk_orth                       
    dL      = dimLead**2
    j_time = 2*(time+1.d-5)/dt + 1
    allocate(hop_plus(L),hop_minus(L),hop(L))
    !+-------------------------+!
    !+- gutzwiller quantities -+!
    !+-------------------------+!
    allocate(gz_phi(L),gz_R(L),x(L))
    do iSlab = 1,L
       ik_sys = Nk_tot*dimk + 3*(iSlab-1)
       gz_phi(iSlab)%p0 = y(ik_sys+1) 
       gz_phi(iSlab)%p1 = y(ik_sys+2) 
       gz_phi(iSlab)%p2 = y(ik_sys+3) 
    end do
    gz_R = GZ_hop(gz_phi)
    x = GZ_doping(gz_phi)
    !initialize k sums
    hop_plus  = Z0
    hop_minus = Z0
    hop       = Z0
    hyb_left  = Z0
    hyb_right = Z0
    ik_sys = 0
    !+- LOOP ON (planar) k POINTS -+!
    do ik=1,Nk_tot
       ek = square_lattice_disp(vec_k(ik))
       ik_sys = (ik-1)*dimk
       ik_sys0 = ik_sys
       allocate(sumGleadR(dimLead),sumGleadL(dimLead),sumGhybL(L),sumGhybR(L))
       !+- LEADS -+!
       do ilead=1,Nk_orth
          kl_i = k_orth(ilead)
          chem_lead = muL(j_time)
          vk_i = vk_L(ilead,j_time)
          side_i = 1
          ekl_i = chain_disp(kl_i)*t_lead - chem_lead
          sumGleadL(ilead) = 0.d0
          sumGleadR(ilead) = 0.d0
          do jlead = 1,Nk_orth
             ik_sys = ik_sys + 1
             kl_j = k_orth(jlead)
             chem_lead = muL(j_time)
             vk_j = vk_L(jlead,j_time)
             sumGleadL(ilead) = sumGleadL(ilead) + vk_j*y(ik_sys)
             side_j = 1
             ekl_j = chain_disp(kl_j)*t_lead - chem_lead
             ihyb_i = ik_sys0 + dL + (side_i-1)*dimLead + jlead
             ihyb_j = ik_sys0 + dL + (side_j-1)*dimLead + ilead
             !+- LEADS GF EQUATIONS OF MOTION  (kp L, kp' L) -+!
             f(ik_sys) = -Zi*(ekl_i-ekl_j)*y(ik_sys) + &
                  Zi*( vk_j*y(ihyb_j)*conjg(gz_R(side_j))+vk_i*conjg(y(ihyb_i))*gz_R(side_i) )  
          end do
          do jlead = Nk_orth+1,2*Nk_orth
             ik_sys = ik_sys + 1
             kl_j = k_orth(jlead-Nk_orth)
             chem_lead = muR(j_time)
             vk_j = vk_R(jlead-Nk_orth,j_time)
             sumGleadR(ilead) = sumGleadR(ilead) + vk_j*y(ik_sys)
             side_j = L
             ekl_j = chain_disp(kl_j)*t_lead - chem_lead
             ihyb_i = ik_sys0 + dL + (side_i-1)*dimLead + jlead
             ihyb_j = ik_sys0 + dL + (side_j-1)*dimLead + ilead
             !+- LEADS GF EQUATIONS OF MOTION (kp L, kp' R)-+!
             f(ik_sys) = -Zi*(ekl_i-ekl_j)*y(ik_sys) + &
                  Zi*( vk_j*y(ihyb_j)*conjg(gz_R(side_j))+vk_i*conjg(y(ihyb_i))*gz_R(side_i) )  
          end do
       end do

       do ilead=Nk_orth+1,2*Nk_orth
          kl_i = k_orth(ilead-Nk_orth)
          chem_lead = muR(j_time)
          vk_i = vk_R(ilead-Nk_orth,j_time)
          side_i = L
          ekl_i = chain_disp(kl_i)*t_lead - chem_lead
          sumGleadL(ilead) = 0.d0
          sumGleadR(ilead) = 0.d0
          do jlead = 1,Nk_orth
             ik_sys = ik_sys + 1
             kl_j = k_orth(jlead)
             chem_lead = muL(j_time)
             vk_j = vk_L(jlead,j_time)
             sumGleadL(ilead) = sumGleadL(ilead) + vk_j*y(ik_sys)
             side_j = 1
             ekl_j = chain_disp(kl_j)*t_lead - chem_lead
             ihyb_i = ik_sys0 + dL + (side_i-1)*dimLead + jlead
             ihyb_j = ik_sys0 + dL + (side_j-1)*dimLead + ilead
             !+- LEADS GF EQUATIONS OF MOTION (kp R, kp' L) -+!
             f(ik_sys) = -Zi*(ekl_i-ekl_j)*y(ik_sys) + &
                  Zi*( vk_j*y(ihyb_j)*conjg(gz_R(side_j))+vk_i*conjg(y(ihyb_i))*gz_R(side_i) )  
          end do
          do jlead = Nk_orth+1,2*Nk_orth
             ik_sys = ik_sys + 1
             kl_j = k_orth(jlead-Nk_orth)
             chem_lead = muR(j_time)
             vk_j = vk_R(jlead-Nk_orth,j_time)
             sumGleadR(ilead) = sumGleadR(ilead) + vk_j*y(ik_sys)
             side_j = L
             ekl_j = chain_disp(kl_j)*t_lead - chem_lead
             ihyb_i = ik_sys0 + dL + (side_i-1)*dimLead + jlead
             ihyb_j = ik_sys0 + dL + (side_j-1)*dimLead + ilead
             !+- LEADS GF EQUATIONS OF MOTION (kp R, kp' R) -+!
             f(ik_sys) = -Zi*(ekl_i-ekl_j)*y(ik_sys) + &
                  Zi*( vk_j*y(ihyb_j)*conjg(gz_R(side_j))+vk_i*conjg(y(ihyb_i))*gz_R(side_i) )  
          end do
       end do
       !+- hybrid green -+!
       do iSlab = 1,L
          sumGhybL(iSlab) = Z0
          sumGhybR(iSlab) = Z0
          do ihyb=1,Nk_orth
             ik_sys = ik_sys+1
             kl_i = k_orth(ihyb)
             chem_lead = muL(j_time)
             vk_i = vk_L(ihyb,j_time)
             side_i = 1
             sumGhybL(iSlab) = sumGhybL(iSlab) + vk_i*y(ik_sys)
             ekl_i = chain_disp(kl_i)*t_lead - chem_lead
             !+--------------------------------+!
             !+- HYBRID GF EQUATION OF MOTION -+!
             !+--------------------------------+!
             select case(lead_type)
             case('3d_tb')
                f(ik_sys) = -Zi*ekl_i*y(ik_sys) - Zi*ek*(t_lead - ABS(gz_R(iSlab))**2)*y(ik_sys) 
             case('generic_bath')
                f(ik_sys) = -Zi*ekl_i*y(ik_sys) + Zi*ek*ABS(gz_R(iSlab))**2*y(ik_sys) 
             end select
             if(iSlab.gt.1) then
                f(ik_sys) = f(ik_sys)-Zi*y(ik_sys-dimLead)*gz_R(iSlab)*conjg(gz_R(iSlab-1)) 
             end if
             if(iSlab.lt.L) then
                f(ik_sys) = f(ik_sys)-Zi*y(ik_sys+dimLead)*gz_R(iSlab)*conjg(gz_R(iSlab+1)) 
             end if
             ihyb_i = ik_sys0 + dimk - L*L +(side_i-1)*L + iSlab
             f(ik_sys) = f(ik_sys)-Zi*vk_i*y(ihyb_i)*gz_R(side_i)
             if(iSlab.eq.1) then
                f(ik_sys) = f(ik_sys) + Zi*sumGleadL(ihyb)*gz_R(iSlab)
             end if
             if(iSlab.eq.L) then
                f(ik_sys) = f(ik_sys) + Zi*sumGleadR(ihyb)*gz_R(iSlab)
             end if
          end do
          do ihyb=Nk_orth+1,2*Nk_orth
             ik_sys = ik_sys+1
             kl_i = k_orth(ihyb-Nk_orth)
             chem_lead = muR(j_time)
             vk_i = vk_R(ihyb-Nk_orth,j_time)
             side_i = L
             sumGhybR(iSlab) = sumGhybR(iSlab) + vk_i*y(ik_sys)
             ekl_i = chain_disp(kl_i)*t_lead - chem_lead
             !+--------------------------------+!
             !+- HYBRID GF EQUATION OF MOTION -+!
             !+--------------------------------+!
             select case(lead_type)
             case('3d_tb')
                f(ik_sys) = -Zi*ekl_i*y(ik_sys) - Zi*ek*(t_lead - ABS(gz_R(iSlab))**2)*y(ik_sys) 
             case('generic_bath')
                f(ik_sys) = -Zi*ekl_i*y(ik_sys) + Zi*ek*ABS(gz_R(iSlab))**2*y(ik_sys) 
             end select
             if(iSlab.gt.1) then
                f(ik_sys) = f(ik_sys)-Zi*y(ik_sys-dimLead)*gz_R(iSlab)*conjg(gz_R(iSlab-1)) 
             end if
             if(iSlab.lt.L) then
                f(ik_sys) = f(ik_sys)-Zi*y(ik_sys+dimLead)*gz_R(iSlab)*conjg(gz_R(iSlab+1)) 
             end if
             ihyb_i = ik_sys0 + dimk - L*L +(side_i-1)*L + iSlab
             f(ik_sys) = f(ik_sys)-Zi*vk_i*y(ihyb_i)*gz_R(side_i)
             if(iSlab.eq.1) then
                f(ik_sys) = f(ik_sys) + Zi*sumGleadL(ihyb)*gz_R(iSlab)
             end if
             if(iSlab.eq.L) then
                f(ik_sys) = f(ik_sys) + Zi*sumGleadR(ihyb)*gz_R(iSlab)
             end if
          end do
       end do
       deallocate(sumGleadL,sumGleadR)
       !+- dot GSlab -+!       
       allocate(Gslab(L,L))
       do iSlab=1,L
          do jSlab=1,L
             ik_sys = ik_sys+1
             Gslab(iSlab,jSlab) = y(ik_sys)
          end do
       end do
       ik_sys = ik_sys - L*L
       do iSlab=1,L
          do jSlab=1,L
             ik_sys = ik_sys+1
             !+-------------------------------+!
             !+- SLAB GF EQUATIONS OF MOTION -+!
             !+-------------------------------+!
             f(ik_sys) = -Zi*ek*( ABS(gz_R(iSlab))**2 - ABS(gz_R(jSlab))**2 )*GSlab(iSlab,jSlab)    
             if(iSlab.gt.1) then
                f(ik_sys) = f(ik_sys) + Zi*conjg(gz_R(iSlab))*gz_R(iSlab-1)*Gslab(iSlab-1,jSlab)
             end if
             if(iSlab.lt.L) then
                f(ik_sys) = f(ik_sys) + Zi*conjg(gz_R(iSlab))*gz_R(iSlab+1)*Gslab(iSlab+1,jSlab)
             end if
             if(jSlab.gt.1) then
                f(ik_sys) = f(ik_sys) - Zi*conjg(gz_R(jSlab-1))*gz_R(jSlab)*Gslab(iSlab,jSlab-1)
             end if
             if(jSlab.lt.L) then
                f(ik_sys) = f(ik_sys) - Zi*conjg(gz_R(jSlab+1))*gz_R(jSlab)*Gslab(iSlab,jSlab+1)
             end if
             !+- hybridization -+!
             if(iSlab.eq.1) then
                f(ik_sys) = f(ik_sys)-Zi*(sumGhybL(jSlab))*conjg(gz_R(iSlab)) 
             end if
             if(iSlab.eq.L) then
                f(ik_sys) = f(ik_sys)-Zi*(sumGhybR(jSlab))*conjg(gz_R(iSlab)) 
             end if
             if(jSlab.eq.1) then
                f(ik_sys) = f(ik_sys)-Zi*conjg(sumGhybL(iSlab))*gz_R(jSlab) 
             end if
             if(jSlab.eq.L) then
                f(ik_sys) = f(ik_sys)-Zi*conjg(sumGhybR(iSlab))*gz_R(jSlab) 
             end if
          end do
       end do
       !+- accumulate k-point sums -+!
       do iSlab = 1,L
          if(iSlab.lt.L) then
             hop_plus(iSlab)  = hop_plus(iSlab)-2.d0*Zi*(Gslab(iSlab,iSlab+1))*wt(ik)
             hop_minus(iSlab) = hop_minus(iSlab)-2.d0*Zi*(Gslab(iSlab+1,iSlab))*wt(ik)
          end if
          hop(iSlab) = hop(iSlab) + 2.d0*(1.d0-Zi*Gslab(iSlab,iSlab))*wt(ik)*ek 
       end do
       hyb_left = hyb_left - 2.d0*Zi*sumGhybL(1)*wt(ik)
       hyb_right = hyb_right - 2.d0*Zi*sumGhybR(L)*wt(ik)
       deallocate(Gslab)
       deallocate(sumGhybR,sumGhybL)
    end do !Nk_tot

    coulomb = 0.d0
    if(dabs(alpha_electrostatic).gt.1.d-5) then
       do lx = 1,L
          coulomb(lx) = 0.d0
          do ly = 1,L
             coulomb(lx) = coulomb(lx) + alpha_electrostatic*abs(dble(lx-ly))*x(ly)
          end do
       end do
    end if
    do iSlab = 1,L
       h00 = 2*ABS(gz_R(iSlab))**2*hop(iSlab)
       h01 = conjg(gz_R(iSlab))*hop(iSlab)
       if(iSlab.gt.1) then
          h00 = h00 - 2.d0*REAL(conjg(gz_R(iSlab-1))*gz_R(iSlab)*hop_minus(iSlab-1))
          h01 = h01 - conjg(gz_R(iSlab-1))*hop_minus(iSlab-1)
       else
          h00 = h00 + 2.d0*REAL(conjg(hyb_left)*gz_R(iSlab))
          h01 = h01 + conjg(hyb_left)
       end if
       if(iSlab.lt.L) then
          h00 = h00 - 2.d0*REAL(conjg(gz_R(iSlab+1))*gz_R(iSlab)*hop_plus(iSlab))
          h01 = h01 - conjg(gz_R(iSlab+1))*hop_plus(iSlab)
       else
          h00 = h00 + 2.d0*REAL(conjg(hyb_right)*gz_R(iSlab))
          h01 = h01 + conjg(hyb_right)
       end if
       h22 = x(iSlab)/(x(iSlab)**2-1.d0)*h00
       h00 = x(iSlab)/(1.d0-x(iSlab)**2)*h00
       h01 = sqrt(2.d0/(1.d0-x(iSlab)**2))*h01 
       h22 = h22 + U*0.5d0  + e_loc(iSlab,j_time) + coulomb(iSlab)
       h00 = h00 + U*0.5d0  - e_loc(iSlab,j_time) - coulomb(iSlab)
       !+- prj_entropy contribution -+!
       P0 = 0.25d0*(1.d0+x(iSlab))**2.d0
       P1 = 0.5d0*(1.d0-x(iSlab)**2.d0)
       P2 = 0.25d0*(1.d0-x(iSlab))**2.d0

       prj_s0 = gz_phi(iSlab)%p0
       prj_s1 = gz_phi(iSlab)%p1 
       prj_s2 = gz_phi(iSlab)%p2

       prj_s0 = prj_s0 - 0.5d0*( abs(gz_phi(iSlab)%p0)**2.d0/P0*(1.d0+x(iSlab)) - abs(gz_phi(iSlab)%p2)**2.d0/P2*(1.d0-x(iSlab)) )*gz_phi(iSlab)%p0
       prj_s0 = prj_s0 + x(iSlab)*abs(gz_phi(iSlab)%p1)**2.d0/P1*gz_phi(iSlab)%p0

       prj_s2 = prj_s2 + 0.5d0*( abs(gz_phi(iSlab)%p0)**2.d0/P0*(1.d0+x(iSlab)) - abs(gz_phi(iSlab)%p2)**2.d0/P2*(1.d0-x(iSlab)) )*gz_phi(iSlab)%p2
       prj_s2 = prj_s2 - x(iSlab)*abs(gz_phi(iSlab)%p1)**2.d0/P1*gz_phi(iSlab)%p2

       if( abs(gz_phi(iSlab)%p0) .gt. 1.d-10) then
          prj_s0 = prj_s0 + log( abs(gz_phi(iSlab)%p0)**2.d0 / P0 )*gz_phi(iSlab)%p0
       end if
       if( abs(gz_phi(iSlab)%p1) .gt. 1.d-10) then
          prj_s1 = prj_s1 + log( abs(gz_phi(iSlab)%p1)**2.d0 / P1 )*gz_phi(iSlab)%p1
       end if
       if( abs(gz_phi(iSlab)%p2) .gt. 1.d-10) then
          prj_s2 = prj_s2 + log( abs(gz_phi(iSlab)%p2)**2.d0 / P2 )*gz_phi(iSlab)%p2
       end if
       !+----------------------------+!
       do igz = 1,3
          ik_sys = ik_sys + 1
          select case(igz)
          case(1)
             f(ik_sys) = -Zi*(h00*gz_phi(iSlab)%p0 + h01*gz_phi(iSlab)%p1)        
             f(ik_sys) = f(ik_sys) - Zi*temp*prj_s0
          case(2)
             f(ik_sys) = -Zi*(conjg(h01)*gz_phi(iSlab)%p0 + h01*gz_phi(iSlab)%p2) 
             f(ik_sys) = f(ik_sys) - Zi*temp*prj_s1
          case(3)
             f(ik_sys) = -Zi*(h22*gz_phi(iSlab)%p2 + conjg(h01)*gz_phi(iSlab)%p1) 
             f(ik_sys) = f(ik_sys) - Zi*temp*prj_s2
          end select
       end do
    end do
    deallocate(gz_phi)
    return
  end function slab_lead_eom_kleads











  function slab_lead_eom_real_space(time,y,Nsys) result(f)
    implicit none
    !inputs
    integer                                     :: Nsys ! nr of equations
    real(8)                                     :: time ! time variable
    complex(8),dimension(Nsys)                  :: y    ! argument array
    complex(8),dimension(Nsys)                  :: f    ! result


    !indeces 
    integer                                     :: ik,ik_sys,ik_sys0
    integer                                     :: dimk,dimLead,dL,dim_hyb
    integer                                     :: ilead,jlead,ihyb

    integer                                     :: side_i,side_j
    integer                                     :: ihyb_i,ihyb_j,ihyb_left,ihyb_right

    integer                                     :: isum_hyb_R
    integer                                     :: isum_hyb_L

    integer                                     :: iSlab,jSlab,kSlab
    integer                                     :: igz

    integer                                     :: j_time
    integer                                     :: lx,ly

    !auxiliary GF
    complex(8),dimension(:,:),allocatable       :: Glead
    complex(8),dimension(:,:),allocatable       :: Gslab


    !k-sums
    complex(8),dimension(:),allocatable         :: sumGleadR,sumGleadL
    complex(8),dimension(:),allocatable         :: sumGhybL,sumGhybR
    complex(8),dimension(:),allocatable         :: sumGhybL_k,sumGhybR_k

    complex(8),dimension(:),allocatable         :: hop_plus,hop_minus,hop
    complex(8)                                  :: hyb_left,hyb_right

    !gutzwiller
    type(gz_projector),dimension(:),allocatable :: gz_phi
    complex(8),dimension(:),allocatable         :: gz_R
    real(8),dimension(:),allocatable            :: x
    complex(8)                                  :: h00,h01,h22
    complex(8)                                  :: prj_s0,prj_s1,prj_s2
    real(8)                                     :: P0,P1,P2

    !dispersions and momenta
    real(8)                                     :: ek
    real(8)                                     :: kl_i,kl_j,ekl_i,ekl_j
    real(8)                                     :: vk_i,vk_j,chem_lead
    !chemical potential
    real(8)                                     :: chem_i,chem_j
    !coulom potential
    real(8),dimension(L)                        :: coulomb
    real(8) :: tnk
    real(8) :: left_t,right_t
    !
    dimk    = 2*Nk_orth*L + L*L
    dimLead = 2*Nk_orth                       
    dL      = dimLead**2
    dim_hyb = dimLead*L
    !
    j_time  = 2*(time+1.d-5)/dt + 1
    left_t  = get_real_space_coupling(time,'L')
    right_t = get_real_space_coupling(time,'R')
    ihyb_left  = Nk_orth
    ihyb_right = Nk_orth + Slab
    !
    allocate(hop_plus(L),hop_minus(L),hop(L))
    !+-------------------------+!
    !+- gutzwiller quantities -+!
    !+-------------------------+!
    allocate(gz_phi(L),gz_R(L),x(L))
    do iSlab = 1,L
       ik_sys = Nk_tot*L*L + 3*(iSlab-1) 
       gz_phi(iSlab)%p0 = y(ik_sys+1) 
       gz_phi(iSlab)%p1 = y(ik_sys+2) 
       gz_phi(iSlab)%p2 = y(ik_sys+3) 
    end do
    gz_R = GZ_hop(gz_phi)
    x = GZ_doping(gz_phi)
    !initialize k sums
    hop_plus  = Z0
    hop_minus = Z0
    hop       = Z0
    hyb_left  = Z0
    hyb_right = Z0
    allocate(sumGleadR(dimLead),sumGleadL(dimLead),sumGhybL(L),sumGhybR(L))
    allocate(Gslab(L,L))
    !+- initializae global index -+!
    ik_sys = 0
    !+- LOOP ON (planar) k POINTS -+!
    do ik=1,Nk_tot
       ek = square_lattice_disp(vec_k(ik))
       do iSlab=1,L
          do jSlab=1,L
             ik_sys = ik_sys+1
             Gslab(iSlab,jSlab) = y(ik_sys)
          end do
       end do
       !+- put back global index -+!
       ik_sys = ik_sys - L*L
       !+-------------------------+!
       do iSlab=1,L
          do jSlab=1,L
             !+- up date global index -+!
             ik_sys = ik_sys+1
             !+-------------------------------+!
             !+- SLAB GF EQUATIONS OF MOTION -+!
             !+-------------------------------+!
             f(ik_sys) = -Zi*ek*( ABS(gz_R(iSlab))**2 - ABS(gz_R(jSlab))**2 )*GSlab(iSlab,jSlab)    

             !+- left hoppings -+!
             if(iSlab.gt.1) then 
                if(iSlab.eq.ihyb_left+1) then
                   !left lead coupling
                   f(ik_sys) = f(ik_sys) + Zi*conjg(gz_R(iSlab))*gz_R(iSlab-1)*Gslab(iSlab-1,jSlab)*left_t
                else
                   if(iSlab.eq.ihyb_right+1) then
                      !right lead coupling
                      f(ik_sys) = f(ik_sys) + Zi*conjg(gz_R(iSlab))*gz_R(iSlab-1)*Gslab(iSlab-1,jSlab)*right_t
                   else
                      f(ik_sys) = f(ik_sys) + Zi*conjg(gz_R(iSlab))*gz_R(iSlab-1)*Gslab(iSlab-1,jSlab)
                   end if
                end if
             end if
             !
             if(jSlab.gt.1) then
                if(jSlab.eq.ihyb_left+1) then
                   ! left led coupling
                   f(ik_sys) = f(ik_sys) - Zi*conjg(gz_R(jSlab-1))*gz_R(jSlab)*Gslab(iSlab,jSlab-1)*left_t
                else
                   if(jSlab.eq.ihyb_right+1) then
                      ! rigth lead coupling
                      f(ik_sys) = f(ik_sys) - Zi*conjg(gz_R(jSlab-1))*gz_R(jSlab)*Gslab(iSlab,jSlab-1)*right_t
                   else
                      f(ik_sys) = f(ik_sys) - Zi*conjg(gz_R(jSlab-1))*gz_R(jSlab)*Gslab(iSlab,jSlab-1)
                   end if
                end if
             end if

             !+- right hoppings -+!
             if(iSlab.lt.L) then
                if(iSlab.eq.ihyb_left) then
                   ! left lead coupling
                   f(ik_sys) = f(ik_sys) + Zi*conjg(gz_R(iSlab))*gz_R(iSlab+1)*Gslab(iSlab+1,jSlab)*left_t
                else
                   if(iSlab.eq.ihyb_right) then
                      ! right lead coupling
                      f(ik_sys) = f(ik_sys) + Zi*conjg(gz_R(iSlab))*gz_R(iSlab+1)*Gslab(iSlab+1,jSlab)*right_t
                   else
                      f(ik_sys) = f(ik_sys) + Zi*conjg(gz_R(iSlab))*gz_R(iSlab+1)*Gslab(iSlab+1,jSlab)
                   end if
                end if
             end if
             if(jSlab.lt.L) then
                if(jSlab.eq.ihyb_left)  then
                   ! left lead coupling
                   f(ik_sys) = f(ik_sys) - Zi*conjg(gz_R(jSlab+1))*gz_R(jSlab)*Gslab(iSlab,jSlab+1)*left_t
                else
                   if(jSlab.eq.ihyb_right) then
                      ! right lead coupling
                      f(ik_sys) = f(ik_sys) - Zi*conjg(gz_R(jSlab+1))*gz_R(jSlab)*Gslab(iSlab,jSlab+1)*right_t
                   else
                      f(ik_sys) = f(ik_sys) - Zi*conjg(gz_R(jSlab+1))*gz_R(jSlab)*Gslab(iSlab,jSlab+1)
                   end if
                end if
             end if

          end do
       end do
       !+- accumulate k-point sums -+!
       do iSlab = 1,L
          if(iSlab.lt.L) then
             hop_plus(iSlab)  = hop_plus(iSlab)-2.d0*Zi*(Gslab(iSlab,iSlab+1))*wt(ik)
             hop_minus(iSlab) = hop_minus(iSlab)-2.d0*Zi*(Gslab(iSlab+1,iSlab))*wt(ik)
          end if
          hop(iSlab) = hop(iSlab) + 2.d0*(1.d0-Zi*Gslab(iSlab,iSlab))*wt(ik)*ek 
       end do
    end do !Nk_tot


    coulomb = 0.d0
    if(dabs(alpha_electrostatic).gt.1.d-5) then
       do lx = 1,L
          coulomb(lx) = 0.d0
          do ly = 1,L
             coulomb(lx) = coulomb(lx) + alpha_electrostatic*abs(dble(lx-ly))*x(ly)
          end do
       end do
    end if

    do iSlab = 1,L
       h00 = 2*ABS(gz_R(iSlab))**2*hop(iSlab)
       h01 = conjg(gz_R(iSlab))*hop(iSlab)
       !+- left hopping -+!
       if(iSlab.gt.1) then
          if(iSlab.eq.ihyb_left+1) then
             ! left lead coupling
             h00 = h00 - 2.d0*REAL(conjg(gz_R(iSlab-1))*gz_R(iSlab)*hop_minus(iSlab-1))*left_t             
             h01 = h01 - conjg(gz_R(iSlab-1))*hop_minus(iSlab-1)*left_t
          else
             if(iSlab.eq.ihyb_right+1) then
                h00 = h00 - 2.d0*REAL(conjg(gz_R(iSlab-1))*gz_R(iSlab)*hop_minus(iSlab-1))*right_t
                h01 = h01 - conjg(gz_R(iSlab-1))*hop_minus(iSlab-1)*right_t
             else
                h00 = h00 - 2.d0*REAL(conjg(gz_R(iSlab-1))*gz_R(iSlab)*hop_minus(iSlab-1))
                h01 = h01 - conjg(gz_R(iSlab-1))*hop_minus(iSlab-1)
             end if
          end if
       end if
       !+- right hopping -+!
       if(iSlab.lt.L) then
          if(iSlab.eq.ihyb_left) then
             h00 = h00 - 2.d0*REAL(conjg(gz_R(iSlab+1))*gz_R(iSlab)*hop_plus(iSlab))*left_t
             h01 = h01 - conjg(gz_R(iSlab+1))*hop_plus(iSlab)*left_t
          else
             if(iSlab.eq.ihyb_right) then
                h00 = h00 - 2.d0*REAL(conjg(gz_R(iSlab+1))*gz_R(iSlab)*hop_plus(iSlab))*right_t
                h01 = h01 - conjg(gz_R(iSlab+1))*hop_plus(iSlab)*right_t
             else
                h00 = h00 - 2.d0*REAL(conjg(gz_R(iSlab+1))*gz_R(iSlab)*hop_plus(iSlab))
                h01 = h01 - conjg(gz_R(iSlab+1))*hop_plus(iSlab)
             end if
          end if
       end if
       h22 = x(iSlab)/(x(iSlab)**2-1.d0)*h00
       h00 = x(iSlab)/(1.d0-x(iSlab)**2)*h00
       h01 = sqrt(2.d0/(1.d0-x(iSlab)**2))*h01 
       !
       h22 = h22 + U*0.5d0  + e_loc(iSlab,j_time) + coulomb(iSlab)
       h00 = h00 + U*0.5d0  - e_loc(iSlab,j_time) - coulomb(iSlab)
       !+- prj_entropy contribution -+!
       P0 = 0.25d0*(1.d0+x(iSlab))**2.d0
       P1 = 0.5d0*(1.d0-x(iSlab)**2.d0)
       P2 = 0.25d0*(1.d0-x(iSlab))**2.d0
       !
       prj_s0 = gz_phi(iSlab)%p0
       prj_s1 = gz_phi(iSlab)%p1 
       prj_s2 = gz_phi(iSlab)%p2
       !
       prj_s0 = prj_s0 - 0.5d0*( abs(gz_phi(iSlab)%p0)**2.d0/P0*(1.d0+x(iSlab)) - abs(gz_phi(iSlab)%p2)**2.d0/P2*(1.d0-x(iSlab)) )*gz_phi(iSlab)%p0
       prj_s0 = prj_s0 + x(iSlab)*abs(gz_phi(iSlab)%p1)**2.d0/P1*gz_phi(iSlab)%p0
       !
       prj_s2 = prj_s2 + 0.5d0*( abs(gz_phi(iSlab)%p0)**2.d0/P0*(1.d0+x(iSlab)) - abs(gz_phi(iSlab)%p2)**2.d0/P2*(1.d0-x(iSlab)) )*gz_phi(iSlab)%p2
       prj_s2 = prj_s2 - x(iSlab)*abs(gz_phi(iSlab)%p1)**2.d0/P1*gz_phi(iSlab)%p2
       !
       if( abs(gz_phi(iSlab)%p0) .gt. 1.d-10) then
          prj_s0 = prj_s0 + log( abs(gz_phi(iSlab)%p0)**2.d0 / P0 )*gz_phi(iSlab)%p0
       end if
       if( abs(gz_phi(iSlab)%p1) .gt. 1.d-10) then
          prj_s1 = prj_s1 + log( abs(gz_phi(iSlab)%p1)**2.d0 / P1 )*gz_phi(iSlab)%p1
       end if
       if( abs(gz_phi(iSlab)%p2) .gt. 1.d-10) then
          prj_s2 = prj_s2 + log( abs(gz_phi(iSlab)%p2)**2.d0 / P2 )*gz_phi(iSlab)%p2
       end if
       !+----------------------------+!
       do igz = 1,3
          !+- update global index -+!
          ik_sys = ik_sys + 1
          !+-----------------------+!
          select case(igz)
          case(1)
             f(ik_sys) = -Zi*(h00*gz_phi(iSlab)%p0 + h01*gz_phi(iSlab)%p1)        
             f(ik_sys) = f(ik_sys) - Zi*temp*prj_s0
          case(2)
             f(ik_sys) = -Zi*(conjg(h01)*gz_phi(iSlab)%p0 + h01*gz_phi(iSlab)%p2) 
             f(ik_sys) = f(ik_sys) - Zi*temp*prj_s1
          case(3)
             f(ik_sys) = -Zi*(h22*gz_phi(iSlab)%p2 + conjg(h01)*gz_phi(iSlab)%p1) 
             f(ik_sys) = f(ik_sys) - Zi*temp*prj_s2
          end select
       end do
    end do
    deallocate(gz_phi)
    deallocate(sumGleadL,sumGleadR)
    deallocate(sumGhybR,sumGhybL)
    deallocate(Gslab)
    return
  end function slab_lead_eom_real_space




  function get_real_space_coupling(t,bath) result(vkp)
    real(8) :: kp,t,r
    real(8) :: vkp
    character(len=1) :: bath
    real(8) :: switch,ekp,nu

    if(t.lt.ramp_time) then
       r = (1.d0 - 1.5d0*cos(pi*t/ramp_time) + 0.5d0*(cos(pi*t/ramp_time))**3)*0.5d0
    else
       r = 1.d0
    end if

    select case(bath)
    case('L')
       vkp = vL
    case('R')
       vkp = vR
    end select
    !vkp = vkp*r
    return
  end function get_real_space_coupling








  !+----------------------------------------+!
  !+- k,time dependent leads-slab coupling -+!
  !+----------------------------------------+!
  function get_bath_coupling(kp,bath,t) result(vkp)
    real(8) :: kp,t,r
    real(8) :: vkp
    character(len=1) :: bath
    real(8) :: switch,ekp,nu

    nu=1.0

    ekp = chain_disp(kp)*t_lead

    if(t.lt.cut_time) then
       if(t.lt.ramp_time) then
          r = (1.d0 - 1.5d0*cos(pi*t/ramp_time) + 0.5d0*(cos(pi*t/ramp_time))**3)*0.5d0
       else
          r = 1.d0
       end if
    else
       r = 0.2d0
    end if

    select case(lead_type)
    case('3d_tb')
       vkp = -1.d0*sqrt(2.d0/dble(Nk_orth+1))*sin(kp)*vL
    case('cubic')
       vkp=1.d0/sqrt(dble(Nk_orth))
    case('generic_bath')
       !+- insert here the possibility of different baths -+!
       select case(i_bath)
       case('flat')
          vkp = 1.d0
       case('flat_smooth')
          vkp = 1.d0/(1.d0+exp(nu*(ekp - W_bath)))/(1.d0+exp(-nu*(ekp + W_bath)))
       case('circular')
          vkp = 4.d0/pi*sqrt(1.d0-ekp**2/W_bath**2)
       end select
       vkp = sqrt(0.5d0*pi*W_bath*Gamma)/sqrt(dble(NK_orth))*sqrt(abs(sin(kp)))*vkp       
    end select
    vkp = vkp*r
    return
  end function get_bath_coupling

  !+-------------------------------------+!
  !+- time dependent chemical potential -+!
  !+-------------------------------------+!
  function chem_bath(bath,t) result(mu)
    real(8)          :: t
    character(len=1) :: bath
    real(8)          :: switch
    real(8)          :: mu

    select case(bath)
    case('L')
       mu = mu_L
    case('R')
       mu = mu_R
    end select

    if(t.lt.time_bias) then
       switch = 0.d0
    else

       if(t.lt.time_bias + ramp_bias) then
          switch = (1.d0 - 1.5d0*cos(pi*(t-time_bias)/ramp_bias) + 0.5d0*(cos(pi*(t-time_bias)/ramp_bias))**3)*0.5d0
       else
          switch = 1.d0
       end if
    end if

    mu = mu*switch


  end function chem_bath




END MODULE EQS_OF_MOTION
