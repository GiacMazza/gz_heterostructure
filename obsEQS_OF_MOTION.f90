MODULE OBS_EQS_OF_MOTION
  USE GLOBAL

CONTAINS

  function slab_lead_obs_eom(time,y,Nsys) result(f)
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


    real(8),dimension(:),allocatable            :: doublons,holons,ph_doublons,ph_holons
    complex(8),dimension(:),allocatable         :: dr_dh,dr_dd,dr_dpd,dr_dph 
    real(8)                                     :: PH_D,PH_H


    dimk    = 2*Nk_orth*L + L*L
    dimLead = 2*Nk_orth                       
    dL      = dimLead**2
    dim_hyb = dimLead*L

    j_time = 2*(time+1.d-5)/dt + 1

    allocate(hop_plus(L),hop_minus(L),hop(L))
    allocate(doublons(L),holons(L),ph_doublons(L),ph_holons(L))
    allocate(dr_dd(L),dr_dh(L),dr_dpd(L),dr_dph(L))

    !+-------------------------+!
    !+- gutzwiller quantities -+!
    !+-------------------------+!
    allocate(gz_phi(L),gz_R(L),x(L))

    PH_D = 0.d0
    PH_H = 0.d0
    do iSlab = 1,L
       ik_sys = Nk_tot*dimk + dimLead*dimLead + 4*(iSlab-1)
       doublons(iSlab)    = dreal(y(ik_sys+1))
       ph_doublons(iSlab) = dreal(y(ik_sys+2))
       holons(iSlab)      = dreal(y(ik_sys+3))
       ph_holons(iSlab)   = dreal(y(ik_sys+4))
       PH_D = PH_D + ph_doublons(iSlab)
       PH_H = PH_H + ph_holons(iSlab)
    end do

    gz_R = GZ_hop_hd(doublons,holons,ph_doublons,ph_holons)
    x    = GZ_doping_hd(doublons,holons,ph_doublons,ph_holons)

    dr_dd  = deriveR_doublon(doublons,holons,ph_doublons,ph_holons)
    dr_dh  = deriveR_holon(doublons,holons,ph_doublons,ph_holons)
    dr_dpd = deriveR_ph_doublon(doublons,holons,ph_doublons,ph_holons)
    dr_dph = deriveR_ph_holon(doublons,holons,ph_doublons,ph_holons)

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
          !+- update Glead sums -+!
          sumGleadR(ilead) = sumGleadR(ilead) + vk_j*y(ik_sys)   
          !+---------------------+!

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
       !+- update Glead sums -+!
       sumGleadL(ilead) = 0.d0
       sumGleadR(ilead) = 0.d0
       !+---------------------+!

       !+- RIGHT-LEFT -+!
       do jlead = 1,Nk_orth
          !+- update global index -+!
          ik_sys = ik_sys + 1
          !+-----------------------+!
          kl_j = k_orth(jlead)
          chem_lead = muL(j_time)
          vk_j = vk_L(jlead,j_time)
          side_j = 1
          ekl_j = chain_disp(kl_j)*t_lead - chem_lead
          !+- update Glead sums -+!
          sumGleadL(ilead) = sumGleadL(ilead) + vk_j*y(ik_sys)
          !+---------------------+!

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
          side_j = L
          ekl_j = chain_disp(kl_j)*t_lead - chem_lead
          !+- update Gleads sums -+!
          sumGleadR(ilead) = sumGleadR(ilead) + vk_j*y(ik_sys)
          !+----------------------+!

          f(ik_sys) = -Zi*(ekl_i-ekl_j)*y(ik_sys) + &
               Zi*vk_j*conjg(gz_R(side_j))*sumGhybR_k(ilead) + &
               Zi*vk_i*gz_R(side_i)*conjg(sumGhybR_k(jlead))
       end do
    end do

    !+- LOOP ON (planar) k POINTS -+!
    do ik=1,Nk_tot
       ek = epsik(ik)
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

       !+- doublons -+!
       ik_sys = ik_sys + 1
       f(ik_sys) = 2.d0*dreal(conjg(gz_R(iSlab))*hop(iSlab)*dr_dpd(iSlab))
       if(iSlab.gt.1) then
          f(ik_sys) = f(ik_sys) - 2.d0*real(conjg(gz_R(iSlab-1))*conjg(hop_plus(iSlab-1))*dr_dpd(iSlab))
       else
          !+- hybridize with left lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*real(conjg(hyb_left)*dr_dpd(iSlab))
          !+----------------------------+!
       end if
       if(iSlab.lt.L) then
          f(ik_sys) = f(ik_sys) - 2.d0*real(conjg(gz_R(iSlab+1))*hop_plus(iSlab)*dr_dpd(iSlab))
       else
          !+- hybridize with right lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*real(conjg(hyb_right)*dr_dpd(iSlab))
          !+-----------------------------+!
       end if
       !f(ik_sys) = -1.d0*f(ik_sys)
       f(ik_sys) = 0.d0

       !+- ph_doublons -+!
       ik_sys = ik_sys + 1
       f(ik_sys) = 2.d0*dreal(conjg(gz_R(iSlab))*hop(iSlab)*dr_dd(iSlab)) + Uz_time(iSlab,j_time)*0.5d0 - e_loc(iSlab,j_time) 
       f(ik_sys) = f(ik_sys) + coulomb(iSlab)
       if(iSlab.gt.1) then
          f(ik_sys) = f(ik_sys) - 2.d0*real(conjg(gz_R(iSlab-1))*conjg(hop_plus(iSlab-1))*dr_dd(iSlab))
       else
          !+- hybridize with left lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*real(conjg(hyb_left)*dr_dd(iSlab))
          !+----------------------------+!
       end if
       if(iSlab.lt.L) then
          f(ik_sys) = f(ik_sys) - 2.d0*real(conjg(gz_R(iSlab+1))*hop_plus(iSlab)*dr_dd(iSlab))
       else 
          !+- hybridize with right lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*real(conjg(hyb_right)*dr_dd(iSlab))
          !+-----------------------------+!
       end if
       !+- INSERT DISSIPATIVE TERM HERE -+!
       !f(ik_sys) = f(ik_sys) - eta_bath*y(ik_sys)
       !+--------------------------------+!
       f(ik_sys) = 0.d0

       !+- holons -+!
       ik_sys = ik_sys + 1
       f(ik_sys) = 2.d0*dreal(conjg(gz_R(iSlab))*hop(iSlab)*dr_dph(iSlab))
       if(iSlab.gt.1) then
          f(ik_sys) = f(ik_sys) - 2.d0*real(conjg(gz_R(iSlab-1))*conjg(hop_plus(iSlab-1))*dr_dph(iSlab))
       else
          !+- hybridize with left lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*real(conjg(hyb_left)*dr_dph(iSlab))
          !+----------------------------+!
       end if
       if(iSlab.lt.L) then
          f(ik_sys) = f(ik_sys) - 2.d0*real(conjg(gz_R(iSlab+1))*hop_plus(iSlab)*dr_dph(iSlab))
       else
          !+- hybridize with right lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*real(conjg(hyb_right)*dr_dph(iSlab))
          !+-----------------------------+!
       end if
       f(ik_sys) = -1.d0*f(ik_sys)
       f(ik_sys) = 0.d0

       !+- ph_holon -+!
       ik_sys = ik_sys + 1
       f(ik_sys) = 2.d0*dreal(conjg(gz_R(iSlab))*hop(iSlab)*dr_dh(iSlab)) + Uz_time(iSlab,j_time)*0.5d0 + e_loc(iSlab,j_time)  
       f(ik_sys) = f(ik_sys) - coulomb(iSlab)
       if(iSlab.gt.1) then
          f(ik_sys) = f(ik_sys) - 2.d0*real(conjg(gz_R(iSlab-1))*conjg(hop_plus(iSlab-1))*dr_dh(iSlab))
       else
          !+- hybridize with left lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*real(conjg(hyb_left)*dr_dh(iSlab))
          !+----------------------------+!
       end if
       if(iSlab.lt.L) then
          f(ik_sys) = f(ik_sys) - 2.d0*real(conjg(gz_R(iSlab+1))*hop_plus(iSlab)*dr_dh(iSlab))
       else
          !+- hybridize with right lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*real(conjg(hyb_right)*dr_dh(iSlab))
          !+-----------------------------+!
       end if
       !+- INSERT DISSIPATIVE TERM HERE -+!
       !f(ik_sys) = f(ik_sys) - eta_bath*y(ik_sys)
       !+--------------------------------+!
       f(ik_sys) = 0.d0

    end do


    deallocate(gz_phi)
    deallocate(sumGleadL,sumGleadR)
    deallocate(sumGhybR,sumGhybL)
    deallocate(Gslab)

    return

  end function slab_lead_obs_eom









  function slab_lead_obs_eom_k(time,y,Nsys) result(f)
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
    real(8),dimension(L)                        :: coulomb,phase


    real(8),dimension(:),allocatable            :: doublons,holons,ph_doublons,ph_holons
    complex(8),dimension(:),allocatable         :: dr_dh,dr_dd,dr_dpd,dr_dph 
    real(8)                                     :: PH_D,PH_H



    dimk    = 2*Nk_orth*(2*Nk_orth+L) + L*L   !nr of GF for each k(2d) point
    dimLead = 2*Nk_orth                       
    dL      = dimLead**2

    j_time = 2*(time+1.d-5)/dt + 1

    !+-------------------------+!
    !+- gutzwiller quantities -+!
    !+-------------------------+!
    allocate(gz_phi(L),gz_R(L),x(L))

    allocate(doublons(L),holons(L),ph_doublons(L),ph_holons(L))
    allocate(dr_dd(L),dr_dh(L),dr_dpd(L),dr_dph(L))

    phase=0.d0
    PH_D = 0.d0
    PH_H = 0.d0
    do iSlab = 1,L
       ik_sys = Nk_tot*dimk +  4*(iSlab-1)
       doublons(iSlab)    = abs(dreal(y(ik_sys+1)))
       ph_doublons(iSlab) = dreal(y(ik_sys+2))
       holons(iSlab)      = abs(dreal(y(ik_sys+3)))
       ph_holons(iSlab)   = dreal(y(ik_sys+4))
       PH_D = PH_D + ph_doublons(iSlab)
       PH_H = PH_H + ph_holons(iSlab)
       if(iSlab.lt.L) then
          phase(iSlab) = (e_loc(iSlab,j_time) - e_loc(iSlab+1,j_time))*time*0.d0
       end if
    end do

    gz_R = GZ_hop_hd(doublons,holons,ph_doublons,ph_holons)
    x    = GZ_doping_hd(doublons,holons,ph_doublons,ph_holons)

    dr_dd  = deriveR_doublon(doublons,holons,ph_doublons,ph_holons)
    dr_dh  = deriveR_holon(doublons,holons,ph_doublons,ph_holons)
    dr_dpd = deriveR_ph_doublon(doublons,holons,ph_doublons,ph_holons)
    dr_dph = deriveR_ph_holon(doublons,holons,ph_doublons,ph_holons)

    allocate(hop_plus(L),hop_minus(L),hop(L))
    !initialize k sums
    hop_plus  = 0.d0
    hop_minus = 0.d0
    hop       = 0.d0
    hyb_left  = 0.d0
    hyb_right = 0.d0

    ik_sys = 0
    !+- LOOP ON (planar) k POINTS -+!
    do ik=1,Nk_tot
       ek =epsik(ik)
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

          sumGhybL(iSlab) = 0.d0
          sumGhybR(iSlab) = 0.d0

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
                f(ik_sys) = f(ik_sys) + Zi*conjg(gz_R(iSlab))*gz_R(iSlab-1)*Gslab(iSlab-1,jSlab)*Exp(Zi*phase(iSlab-1)) !
             end if

             if(iSlab.lt.L) then
                f(ik_sys) = f(ik_sys) + Zi*conjg(gz_R(iSlab))*gz_R(iSlab+1)*Gslab(iSlab+1,jSlab)*Exp(-Zi*phase(iSlab))
             end if

             if(jSlab.gt.1) then
                f(ik_sys) = f(ik_sys) - Zi*conjg(gz_R(jSlab-1))*gz_R(jSlab)*Gslab(iSlab,jSlab-1)*Exp(-Zi*phase(jSlab-1))  !
             end if

             if(jSlab.lt.L) then
                f(ik_sys) = f(ik_sys) - Zi*conjg(gz_R(jSlab+1))*gz_R(jSlab)*Gslab(iSlab,jSlab+1)*Exp(Zi*phase(jSlab))
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

       hyb_left  = hyb_left - 2.d0*Zi*sumGhybL(1)*wt(ik)
       hyb_right = hyb_right - 2.d0*Zi*sumGhybR(L)*wt(ik)

       deallocate(Gslab)
       deallocate(sumGhybR,sumGhybL)

    end do !Nk_tot

    hyb_right=conjg(hyb_right)
    hyb_left=conjg(hyb_left)

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
       !+- doublons -+!
       ik_sys = ik_sys + 1
       f(ik_sys) = 2.d0*dreal(conjg(gz_R(iSlab))*hop(iSlab)*dr_dpd(iSlab))
       if(iSlab.gt.1) then
          f(ik_sys) = f(ik_sys) - 2.d0*dreal(conjg(gz_R(iSlab-1))*conjg(hop_plus(iSlab-1))*dr_dpd(iSlab)*Exp(-Zi*phase(iSlab-1))) !
       else
          !+- hybridize with left lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*dreal(hyb_left*dr_dpd(iSlab))
          !+----------------------------+!
       end if
       !
       if(iSlab.lt.L) then
          f(ik_sys) = f(ik_sys) - 2.d0*dreal(conjg(gz_R(iSlab+1))*hop_plus(iSlab)*dr_dpd(iSlab)*Exp(Zi*phase(iSlab)))
       else
          !+- hybridize with right lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*dreal(hyb_right*dr_dpd(iSlab))
          !+-----------------------------+!
       end if
       f(ik_sys) = -1.d0*f(ik_sys)
       f(ik_sys) = cmplx(dreal(f(ik_sys)),0.d0)

       !+- ph_doublons -+!
       ik_sys = ik_sys + 1
       f(ik_sys) = 2.d0*dreal(conjg(gz_R(iSlab))*hop(iSlab)*dr_dd(iSlab)) + Uz_time(iSlab,j_time)*0.5d0 - e_loc(iSlab,j_time) 
       f(ik_sys) = f(ik_sys) + coulomb(iSlab)
       if(iSlab.gt.1) then
          f(ik_sys) = f(ik_sys) - 2.d0*dreal(conjg(gz_R(iSlab-1))*conjg(hop_plus(iSlab-1))*dr_dd(iSlab)*Exp(-Zi*phase(iSlab-1)))
       else
          !+- hybridize with left lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*dreal(hyb_left*dr_dd(iSlab))
          !+----------------------------+!
       end if
       if(iSlab.lt.L) then
          f(ik_sys) = f(ik_sys) - 2.d0*dreal(conjg(gz_R(iSlab+1))*hop_plus(iSlab)*dr_dd(iSlab)*Exp(Zi*phase(iSlab)))
       else 
          !+- hybridize with right lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*dreal(hyb_right*dr_dd(iSlab))
          !+-----------------------------+!
       end if
       f(ik_sys) = cmplx(dreal(f(ik_sys)),0.d0)

       !+- holons -+!
       ik_sys = ik_sys + 1
       f(ik_sys) = 2.d0*dreal(conjg(gz_R(iSlab))*hop(iSlab)*dr_dph(iSlab))
       if(iSlab.gt.1) then
          f(ik_sys) = f(ik_sys) - 2.d0*dreal(conjg(gz_R(iSlab-1))*conjg(hop_plus(iSlab-1))*dr_dph(iSlab)*Exp(-Zi*phase(iSlab-1)))
       else
          !+- hybridize with left lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*dreal(hyb_left*dr_dph(iSlab))
          !+----------------------------+!
       end if
       if(iSlab.lt.L) then
          f(ik_sys) = f(ik_sys) - 2.d0*dreal(conjg(gz_R(iSlab+1))*hop_plus(iSlab)*dr_dph(iSlab)*Exp(Zi*phase(iSlab)))
       else
          !+- hybridize with right lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*dreal(hyb_right*dr_dph(iSlab))
          !+-----------------------------+!
       end if
       f(ik_sys) = -1.d0*f(ik_sys)
       f(ik_sys) = cmplx(dreal(f(ik_sys)),0.d0)

       !+- ph_holon -+!
       ik_sys = ik_sys + 1
       f(ik_sys) = 2.d0*dreal(conjg(gz_R(iSlab))*hop(iSlab)*dr_dh(iSlab)) + Uz_time(iSlab,j_time)*0.5d0 + e_loc(iSlab,j_time)  
       f(ik_sys) = f(ik_sys) - coulomb(iSlab)
       if(iSlab.gt.1) then
          f(ik_sys) = f(ik_sys) - 2.d0*dreal(conjg(gz_R(iSlab-1))*conjg(hop_plus(iSlab-1))*dr_dh(iSlab)*Exp(-Zi*phase(iSlab-1)))
       else
          !+- hybridize with left lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*dreal(hyb_left*dr_dh(iSlab))
          !+----------------------------+!
       end if
       if(iSlab.lt.L) then
          f(ik_sys) = f(ik_sys) - 2.d0*dreal(conjg(gz_R(iSlab+1))*hop_plus(iSlab)*dr_dh(iSlab)*Exp(Zi*phase(iSlab)))
       else
          !+- hybridize with right lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*dreal(hyb_right*dr_dh(iSlab))
          !+-----------------------------+!
       end if
       f(ik_sys) = cmplx(dreal(f(ik_sys)),0.d0)

    end do
    deallocate(gz_phi)
    return

  end function slab_lead_obs_eom_k






  function slab_lead_obs_eom_Hij(time,y,Nsys) result(f)
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
    real(8),dimension(L)                        :: coulomb,phase


    real(8),dimension(:),allocatable            :: doublons,holons,ph_doublons,ph_holons
    complex(8),dimension(:),allocatable         :: dr_dh,dr_dd,dr_dpd,dr_dph 
    real(8)                                     :: PH_D,PH_H



    dimk    = 2*Nk_orth*(2*Nk_orth+L) + L*L   !nr of GF for each k(2d) point
    dimLead = 2*Nk_orth                       
    dL      = dimLead**2

    j_time = 2*(time+1.d-5)/dt + 1

    !+-------------------------+!
    !+- gutzwiller quantities -+!
    !+-------------------------+!
    allocate(gz_phi(L),gz_R(L),x(L))

    allocate(doublons(L),holons(L),ph_doublons(L),ph_holons(L))
    allocate(dr_dd(L),dr_dh(L),dr_dpd(L),dr_dph(L))

    phase=0.d0
    PH_D = 0.d0
    PH_H = 0.d0
    do iSlab = 1,L
       ik_sys = Nk_tot*dimk +  4*(iSlab-1)
       doublons(iSlab)    = abs(dreal(y(ik_sys+1)))
       ph_doublons(iSlab) = dreal(y(ik_sys+2))
       holons(iSlab)      = abs(dreal(y(ik_sys+3)))
       ph_holons(iSlab)   = dreal(y(ik_sys+4))
       PH_D = PH_D + ph_doublons(iSlab)
       PH_H = PH_H + ph_holons(iSlab)
       if(iSlab.lt.L) then
          phase(iSlab) = (e_loc(iSlab,j_time) - e_loc(iSlab+1,j_time))*time*0.d0
       end if
    end do

    gz_R = GZ_hop_hd(doublons,holons,ph_doublons,ph_holons)
    x    = GZ_doping_hd(doublons,holons,ph_doublons,ph_holons)

    dr_dd  = deriveR_doublon(doublons,holons,ph_doublons,ph_holons)
    dr_dh  = deriveR_holon(doublons,holons,ph_doublons,ph_holons)
    dr_dpd = deriveR_ph_doublon(doublons,holons,ph_doublons,ph_holons)
    dr_dph = deriveR_ph_holon(doublons,holons,ph_doublons,ph_holons)

    allocate(hop_plus(L),hop_minus(L),hop(L))
    !initialize k sums
    hop_plus  = 0.d0
    hop_minus = 0.d0
    hop       = 0.d0
    hyb_left  = 0.d0
    hyb_right = 0.d0

    ik_sys = 0
    !+- LOOP ON (planar) k POINTS -+!
    do ik=1,Nk_tot
       ek =epsik(ik)
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

          sumGhybL(iSlab) = 0.d0
          sumGhybR(iSlab) = 0.d0

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
                f(ik_sys) = f(ik_sys)-Zi*y(ik_sys-dimLead)*gz_R(iSlab)*conjg(gz_R(iSlab-1))*(-1.*Hslab(iSlab-1,iSlab)) 
             end if

             if(iSlab.lt.L) then
                f(ik_sys) = f(ik_sys)-Zi*y(ik_sys+dimLead)*gz_R(iSlab)*conjg(gz_R(iSlab+1))*(-1.*Hslab(iSlab+1,iSlab)) 
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
                f(ik_sys) = f(ik_sys)-Zi*y(ik_sys-dimLead)*gz_R(iSlab)*conjg(gz_R(iSlab-1))*(-1.*Hslab(iSlab-1,iSlab)) 
             end if

             if(iSlab.lt.L) then
                f(ik_sys) = f(ik_sys)-Zi*y(ik_sys+dimLead)*gz_R(iSlab)*conjg(gz_R(iSlab+1))*(-1.*Hslab(iSlab+1,iSlab)) 
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
                f(ik_sys) = f(ik_sys) + Zi*conjg(gz_R(iSlab))*gz_R(iSlab-1)*Gslab(iSlab-1,jSlab)*Exp(Zi*phase(iSlab-1))*(-1.*Hslab(iSlab,iSlab-1)) !
             end if

             if(iSlab.lt.L) then
                f(ik_sys) = f(ik_sys) + Zi*conjg(gz_R(iSlab))*gz_R(iSlab+1)*Gslab(iSlab+1,jSlab)*Exp(-Zi*phase(iSlab))*(-1.*Hslab(iSlab,iSlab+1))
             end if

             if(jSlab.gt.1) then
                f(ik_sys) = f(ik_sys) - Zi*conjg(gz_R(jSlab-1))*gz_R(jSlab)*Gslab(iSlab,jSlab-1)*Exp(-Zi*phase(jSlab-1))*(-1.*Hslab(jSlab-1,jSlab))  !
             end if

             if(jSlab.lt.L) then
                f(ik_sys) = f(ik_sys) - Zi*conjg(gz_R(jSlab+1))*gz_R(jSlab)*Gslab(iSlab,jSlab+1)*Exp(Zi*phase(jSlab))*(-1.*Hslab(jSlab+1,jSlab))
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

       hyb_left  = hyb_left - 2.d0*Zi*sumGhybL(1)*wt(ik)
       hyb_right = hyb_right - 2.d0*Zi*sumGhybR(L)*wt(ik)

       deallocate(Gslab)
       deallocate(sumGhybR,sumGhybL)

    end do !Nk_tot

    hyb_right=conjg(hyb_right)
    hyb_left=conjg(hyb_left)

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
       !+- doublons -+!
       ik_sys = ik_sys + 1
       f(ik_sys) = 2.d0*dreal(conjg(gz_R(iSlab))*hop(iSlab)*dr_dpd(iSlab))
       if(iSlab.gt.1) then
          f(ik_sys) = f(ik_sys) - 2.d0*dreal(conjg(gz_R(iSlab-1))*conjg(hop_plus(iSlab-1))*dr_dpd(iSlab)*Exp(-Zi*phase(iSlab-1))*(-1.*Hslab(iSlab-1,iSlab))) !
       else
          !+- hybridize with left lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*dreal(hyb_left*dr_dpd(iSlab))
          !+----------------------------+!
       end if
       !
       if(iSlab.lt.L) then
          f(ik_sys) = f(ik_sys) - 2.d0*dreal(conjg(gz_R(iSlab+1))*hop_plus(iSlab)*dr_dpd(iSlab)*Exp(Zi*phase(iSlab))*(-1.*Hslab(iSlab+1,iSlab)))
       else
          !+- hybridize with right lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*dreal(hyb_right*dr_dpd(iSlab))
          !+-----------------------------+!
       end if
       f(ik_sys) = -1.d0*f(ik_sys)
       f(ik_sys) = cmplx(dreal(f(ik_sys)),0.d0)

       !+- ph_doublons -+!
       ik_sys = ik_sys + 1
       f(ik_sys) = 2.d0*dreal(conjg(gz_R(iSlab))*hop(iSlab)*dr_dd(iSlab)) + Uz_time(iSlab,j_time)*0.5d0 - e_loc(iSlab,j_time) 
       f(ik_sys) = f(ik_sys) + coulomb(iSlab)
       if(iSlab.gt.1) then
          f(ik_sys) = f(ik_sys) - 2.d0*dreal(conjg(gz_R(iSlab-1))*conjg(hop_plus(iSlab-1))*dr_dd(iSlab)*Exp(-Zi*phase(iSlab-1))*(-1.*Hslab(iSlab-1,iSlab)))
       else
          !+- hybridize with left lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*dreal(hyb_left*dr_dd(iSlab))
          !+----------------------------+!
       end if
       if(iSlab.lt.L) then
          f(ik_sys) = f(ik_sys) - 2.d0*dreal(conjg(gz_R(iSlab+1))*hop_plus(iSlab)*dr_dd(iSlab)*Exp(Zi*phase(iSlab))*(-1.*Hslab(iSlab+1,iSlab)))
       else 
          !+- hybridize with right lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*dreal(hyb_right*dr_dd(iSlab))
          !+-----------------------------+!
       end if
       f(ik_sys) = cmplx(dreal(f(ik_sys)),0.d0)

       !+- holons -+!
       ik_sys = ik_sys + 1
       f(ik_sys) = 2.d0*dreal(conjg(gz_R(iSlab))*hop(iSlab)*dr_dph(iSlab))
       if(iSlab.gt.1) then
          f(ik_sys) = f(ik_sys) - 2.d0*dreal(conjg(gz_R(iSlab-1))*conjg(hop_plus(iSlab-1))*dr_dph(iSlab)*Exp(-Zi*phase(iSlab-1))*(-1.*Hslab(iSlab-1,iSlab)))
       else
          !+- hybridize with left lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*dreal(hyb_left*dr_dph(iSlab))
          !+----------------------------+!
       end if
       if(iSlab.lt.L) then
          f(ik_sys) = f(ik_sys) - 2.d0*dreal(conjg(gz_R(iSlab+1))*hop_plus(iSlab)*dr_dph(iSlab)*Exp(Zi*phase(iSlab))*(-1.*Hslab(iSlab+1,iSlab)))
       else
          !+- hybridize with right lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*dreal(hyb_right*dr_dph(iSlab))
          !+-----------------------------+!
       end if
       f(ik_sys) = -1.d0*f(ik_sys)
       f(ik_sys) = cmplx(dreal(f(ik_sys)),0.d0)

       !+- ph_holon -+!
       ik_sys = ik_sys + 1
       f(ik_sys) = 2.d0*dreal(conjg(gz_R(iSlab))*hop(iSlab)*dr_dh(iSlab)) + Uz_time(iSlab,j_time)*0.5d0 + e_loc(iSlab,j_time)  
       f(ik_sys) = f(ik_sys) - coulomb(iSlab)
       if(iSlab.gt.1) then
          f(ik_sys) = f(ik_sys) - 2.d0*dreal(conjg(gz_R(iSlab-1))*conjg(hop_plus(iSlab-1))*dr_dh(iSlab)*Exp(-Zi*phase(iSlab-1))*(-1.*Hslab(iSlab-1,iSlab)))
       else
          !+- hybridize with left lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*dreal(hyb_left*dr_dh(iSlab))
          !+----------------------------+!
       end if
       if(iSlab.lt.L) then
          f(ik_sys) = f(ik_sys) - 2.d0*dreal(conjg(gz_R(iSlab+1))*hop_plus(iSlab)*dr_dh(iSlab)*Exp(Zi*phase(iSlab))*(-1.*Hslab(iSlab+1,iSlab)))
       else
          !+- hybridize with right lead -+!
          f(ik_sys) = f(ik_sys) + 2.d0*dreal(hyb_right*dr_dh(iSlab))
          !+-----------------------------+!
       end if
       f(ik_sys) = cmplx(dreal(f(ik_sys)),0.d0)
       
    end do
    deallocate(gz_phi)
    return

  end function slab_lead_obs_eom_Hij


  

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
       r = 1.d0
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
    real(8)          :: pulse_ramp
    
    pulse_ramp=0.5d0
    select case(bath)
    case('L')
       if(t.lt.pulse_time.or.t.gt.pulse_time+cut_time) then
          mu = mu_L
       else
          if(t.lt.pulse_time+pulse_ramp) then
             mu = mu_L + (mu_pulseL-mu_L)*ramp(t,pulse_time,pulse_ramp)
          else
             if(t.lt.pulse_time+cut_time-pulse_ramp) then
                mu = mu_pulseL
             else
                mu = mu_pulseL + (mu_L - mu_pulseL)*ramp(t,pulse_time+cut_time-pulse_ramp,pulse_ramp)
             end if
          end if
       end if
       ! mu = mu_L
       ! if(t.gt.pulse_time.and.t.lt.pulse_time+cut_time) mu=mu_pulseL
    case('R')
       if(t.lt.pulse_time.or.t.gt.pulse_time+cut_time) then
          mu = mu_R
       else
          if(t.lt.pulse_time+pulse_ramp) then
             mu = mu_R + (mu_pulseR-mu_R)*ramp(t,pulse_time,pulse_ramp)
          else
             if(t.lt.pulse_time+cut_time-pulse_ramp) then
                mu = mu_pulseR
             else
                mu = mu_pulseR + (mu_R - mu_pulseR)*ramp(t,pulse_time+cut_time-pulse_ramp,pulse_ramp)
             end if
          end if
       end if
       ! mu = mu_R
       ! if(t.gt.pulse_time.and.t.lt.pulse_time+cut_time) mu=mu_pulseR
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
  
  function ramp(t,t_star,t_ramp) result(r)
    real(8) :: t,t_star,t_ramp,r
    r = 0.5d0*(1.d0 - 1.5d0*cos(pi*(t-t_star)/t_ramp) + 0.5d0*(cos(pi*(t-t_star)/t_ramp))**3)
    return
  end function ramp


END MODULE OBS_EQS_OF_MOTION
