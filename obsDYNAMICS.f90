!+-------------------------------------------------+!
!+- Time evolution of the Gutzwiller wavefunction -+!
! Gz projectors      ---> effective Schr. eq.       !
! Slater determinant ---> Kadanoff-Baym eqs.        !
! Author G Mazza                                    !
!+-------------------------------------------------+!
MODULE OBS_DYNAMICS
  USE GLOBAL
  use RK_IDE
  USE OBS_EQS_OF_MOTION
  implicit none
  private


  integer                                       ::  ik,jk,igrid
  type(vec2D)                                   ::  k
  real(8)                                       ::  ek

  complex(8),dimension(:),allocatable           :: psi_t
  integer                                       :: Nsys
  real(8)                                       :: t


  !+- time dependent observables -+!
  real(8),dimension(:,:,:),allocatable          :: nk
  real(8),dimension(:,:,:,:),allocatable        :: nk_bz
  complex(8),dimension(:,:),allocatable         :: hop_plus,hop_minus
  real(8),dimension(:,:),allocatable            :: j_layer
  real(8),dimension(:,:),allocatable            :: n_dot
  complex(8),dimension(:,:),allocatable         :: hyb_left,hyb_right
  real(8),dimension(:,:),allocatable            :: nSlab
  real(8),dimension(:,:),allocatable            :: eSlab
  real(8),dimension(:,:),allocatable            :: x_gz
  complex(8),dimension(:,:),allocatable         :: r_gz
  real(8),dimension(:,:),allocatable            :: Docc
  real(8),dimension(:,:),allocatable            :: norm_gz
  type(gz_projector),dimension(:,:),allocatable :: gz_phi
  real(8),dimension(:,:),allocatable            :: doublons,holons,ph_doublons,ph_holons

  public :: solve_observables_dynamics

CONTAINS


  SUBROUTINE solve_observables_dynamics

    call allocate_dynamics
    call initial_condition
    call evolve

  CONTAINS

    !+- allocate dynamics -+!
    subroutine allocate_dynamics
      character(len=10) :: fileout
      integer           :: iL
      !+- system -+!
      allocate(t_grid(Nt),t_finer(2*Nt+1))
      !+- time grid -+!
      if(Nt.gt.1) then
         t_grid = linspace(0.d0,dt*real(Nt-1,8),Nt)
      else
         t_grid = linspace(0.d0,dt*real(Nt,8),Nt)
      end if
      t_finer = linspace(0.d0,0.5d0*dt*real(2*Nt,8),2*Nt+1)
      

      !+- THIS IS THE BIG CHANGE -+!
      Nsys = Nk_tot*(2*Nk_orth*(2*Nk_orth+L)+L*L) + 4*L
      !+--------------------------+!

      !write(*,*) Nk_tot,Nk_orth
      
      allocate(psi_t(Nsys))
      write(*,*) 'allocated solution',Nsys,'(',2*Nk_orth*(2*Nk_orth+L)+L*L,'x',Nk_tot,')'      
      !+- time dependent observables -+!
      allocate(nk(Nk_tot,Nt,L),nSlab(Nt,L),hop_plus(Nt,L),hop_minus(Nt,L),hyb_left(Nt,L),hyb_right(Nt,L),eSlab(Nt,L))
      allocate(j_layer(Nt,L),n_dot(Nt,L))
      !+- gutzwiller observables -+!
      allocate(gz_phi(Nt,L),x_gz(Nt,L),r_gz(Nt,L),Docc(Nt,L),norm_gz(Nt,L))
      allocate(doublons(Nt,L),holons(Nt,L),ph_doublons(Nt,L),ph_holons(Nt,L))
      
      
      allocate(vk_L(Nk_orth,2*Nt+1),vk_R(Nk_orth,2*Nt+1))
      allocate(muL(2*Nt+1),muR(2*Nt+1))
      allocate(e_loc(L,2*Nt+1),Uz_time(L,2*Nt+1))
      

      !+- OPEN OUTPUT FILES -+!
      open(50,file="columns_Layer_info.out")
      write(50,"(A1,A17,16A18)")"#","1t","2n","3x","4n_dot","5E_i","6|R_i|**2","7-8R_i","9Docc_i","10-11d^+_{i+1}d_i","12-13R^+(i+1)*R(i)","14-17L-RHyb"
      close(50)
      
      do iL=1,L
         write(fileout,'(I3.3)') iL
         open(unit=200+iL,file='Layer_'//trim(fileout)//'.out')
         open(unit=500+iL,file='Gz_layer_'//trim(fileout)//'.out')
      end do
      
      open(unit=10,file='Slab.out')
      
    end subroutine allocate_dynamics



    subroutine evolve
      integer :: it
      integer :: dim,ik,iks,islab,ilayer
      logical :: iprint
      integer :: ik_sys,dimk

      integer :: t0,t_run

      iprint = .false.

      !+- START TIME EVOLUTION -+!
      call system_clock(count=t0)
      open(100,file='time_loop.out')
      
      call start_timer
      
      t = -dt
      do it = 1,Nt
         t = t + dt
         call test_obs(it) !+- compute observables -+!
         if(mod(it-1,1).eq.0)  call print_dynamics(it)

         !+- extrema ratio -+!
         dimk    = 2*Nk_orth*(2*Nk_orth+L) + L*L
         do iSlab = 1,L
            ik_sys = Nk_tot*dimk +  4*(iSlab-1)
            psi_t(ik_sys+1)    = cmplx(abs(dreal(psi_t(ik_sys+1))),0.d0)
            psi_t(ik_sys+3)    = cmplx(abs(dreal(psi_t(ik_sys+3))),0.d0) 
            ! psi_t(ik_sys+1)    = abs(dreal(psi_t(ik_sys+1)))
            ! psi_t(ik_sys+3)    = abs(dreal(psi_t(ik_sys+3)))
         end do

         psi_t = RK_step(Nsys,mrk,dt,t,psi_t,slab_lead_obs_eom_Hij)
         !psi_t = RK_step(Nsys,mrk,dt,t,psi_t,slab_lead_obs_eom_extreme)
         call system_clock(count=t_run)
         write(100,'(I4,F18.10)') it,log(dble(t_run-t0)/10000.d0)
         call eta(it,Nt)
      end do
      !+- STOP TIME EVOLUTION -+!
      call stop_timer
      call print_dynamics(it)
    end subroutine evolve




    subroutine test_obs(it)
      integer            :: it
      integer            :: ilayer
      integer            :: dim,iSlab,iL
      integer            :: ihyb,ik_sys0,ik_sys
      integer            :: ik_sys_hop
      integer            :: igz
      integer            :: ileads
      integer            :: ikx,iky,j_time
      complex(8)         :: hyb_k
      complex(8),dimension(L)         :: hyb_kL,hyb_kR
      complex(8)         :: hop_k
      real(8)            :: vk,kp,ek
      real(8)            :: time
      real(8),dimension(L) :: phase 

      time=t_grid(it)
      j_time = 2*(time+1.d-5)/dt + 1

      igz = Nsys - 4*L

      !+- change stuff here -+!
      do ilayer = 1,L
         doublons(it,ilayer)    = dreal(psi_t(igz+1))
         ph_doublons(it,ilayer) = dreal(psi_t(igz+2))
         holons(it,ilayer)      = dreal(psi_t(igz+3))
         ph_holons(it,ilayer)   = dreal(psi_t(igz+4))
         igz = igz + 4
      end do

      !tmp
      ! if(iSlab.lt.L) then
      !    phase(iSlab) = (e_loc(iSlab,j_time) - e_loc(iSlab+1,j_time))*time
      ! end if

      r_gz(it,:)    = GZ_hop_hd(doublons(it,:),holons(it,:),ph_doublons(it,:),ph_holons(it,:))
      x_gz(it,:)    = GZ_doping_hd(doublons(it,:),holons(it,:),ph_doublons(it,:),ph_holons(it,:))

      !+- THIS IS THE BIG CHANGE
      dim = 2*Nk_orth*(2*Nk_orth+L)+L*L

      !+- compute one-body opreators on the uncorrelated wave-function -+!
      nSlab(it,:)     = 0.d0
      eSlab(it,:)     = 0.d0
      hyb_right(it,:) = 0.d0
      hyb_left(it,:)  = 0.d0
      hop_plus(it,:)  = 0.d0
      j_layer(it,:)   = 0.d0
      n_dot(it,:)     = 0.d0
      time = (it-1)*dt
      do ik = 1,Nk_tot
         ek = epsik(ik)
         !+- HYBRIDIZATIONS -+!
         ik_sys0 = (ik-1)*dim + (2*Nk_orth)**2   
         do iL=1,L
            !+- left hybridization -+!
            ik_sys = ik_sys0 + (iL-1)*2*Nk_orth
            hyb_k = 0.d0
            hyb_kL(iL) = 0.d0            
            do ihyb=1,Nk_orth
               ik_sys = ik_sys + 1
               kp = k_orth(ihyb)
               vk = get_bath_coupling(kp,'L',time)
               hyb_kL(iL) = hyb_kL(iL) - vk*psi_t(ik_sys)
            end do
            hyb_left(it,iL) = hyb_left(it,iL) + 2.d0*Zi*hyb_kL(iL)*wt(ik)
            !+- right hybridization -+!
            ik_sys = ik_sys0 + (iL-1)*2*Nk_orth + Nk_orth
            hyb_k = 0.d0
            hyb_kR(iL) = 0.d0            
            do ihyb=1,Nk_orth
               ik_sys = ik_sys + 1
               kp = k_orth(ihyb)
               vk = get_bath_coupling(kp,'R',time)
               hyb_kR(iL) = hyb_kR(iL) - vk*psi_t(ik_sys)
            end do
            hyb_right(it,iL) = hyb_right(it,iL) + 2.d0*Zi*hyb_kR(iL)*wt(ik)
         end do

         !+- SLAB GREEN FUNCTION -+!
         ik_sys0 = ik_sys0 + 2*Nk_orth*L
         do ilayer = 1,L
            ik_sys = ik_sys0 + (ilayer-1)*L + ilayer
            ik_sys_hop = ik_sys + 1

            nk(ik,it,ilayer) = 1-Zi*psi_t(ik_sys)
            nSlab(it,ilayer) = nSlab(it,ilayer) + nk(ik,it,ilayer)*wt(ik)
            eSlab(it,ilayer) = eSlab(it,ilayer) + nk(ik,it,ilayer)*wt(ik)*ek
            if(ilayer.lt.L) then
               hop_plus(it,ilayer) = hop_plus(it,ilayer) + 2.d0*Zi*psi_t(ik_sys_hop)*wt(ik)
               j_layer(it,ilayer) = j_layer(it,ilayer) + &
                    2.d0*AIMAG(conjg(r_gz(it,ilayer+1))*r_gz(it,ilayer)*Zi*psi_t(ik_sys_hop))*wt(ik)
               n_dot(it,ilayer) = n_dot(it,ilayer) - &
                    2.d0*AIMAG(conjg(r_gz(it,ilayer+1))*r_gz(it,ilayer)*Zi*psi_t(ik_sys_hop))*wt(ik)
            else
               n_dot(it,ilayer) = n_dot(it,ilayer) - 2.d0*AIMAG(Zi*r_gz(it,ilayer)*conjg(hyb_kR(ilayer)))*wt(ik)
            end if
            ik_sys_hop = ik_sys_hop - 2
            if(ilayer.gt.1) then
               n_dot(it,ilayer) = n_dot(it,ilayer) - &
                    2.d0*AIMAG(conjg(r_gz(it,ilayer-1))*r_gz(it,ilayer)*Zi*psi_t(ik_sys_hop))*wt(ik)
            else
               n_dot(it,ilayer) = n_dot(it,ilayer) - 2.d0*AIMAG( r_gz(it,ilayer)*conjg(Zi*hyb_kL(ilayer)))*wt(ik)
            end if
         end do
      end do

    end subroutine test_obs


    subroutine print_dynamics(it)
      integer,intent(in) :: it
      integer           :: iL,ik,itau
      real(8)           :: slab_ene,hyb_ene,kin_ene,int_ene
      real(8)           :: slab_current
      character(len=10) :: fileout
      real(8)           :: t
      real(8),dimension(L)    :: xSlab
      real(8),dimension(Nk_tot)    :: kx,ky

      slab_ene = 0.d0
      slab_current = 0.d0
      hyb_ene = 0.d0
      kin_ene = 0.d0
      int_ene = 0.d0


      if(it.le.Nt) then


         do iL=1,L-1

            slab_ene     = slab_ene + 2.d0*eSlab(it,iL)*abs(r_gz(it,iL))**2 + U*doublons(it,iL)
            kin_ene     = kin_ene + 2.d0*eSlab(it,iL)*abs(r_gz(it,iL))**2 
            !
            slab_ene     = slab_ene + 2.d0*dREAL(hop_plus(it,iL)*conjg(r_gz(it,iL+1))*r_gz(it,iL))
            kin_ene     = kin_ene + 2.d0*dREAL(hop_plus(it,iL)*conjg(r_gz(it,iL+1))*r_gz(it,iL))
            !
            int_ene = int_ene + U*doublons(it,iL)
            !
            slab_current = slab_current + 2.d0*dIMAG(hop_plus(it,iL)*conjg(r_gz(it,iL+1))*r_gz(it,iL))

            if(iL.eq.1) then
               hyb_ene  = hyb_ene  + 2.d0*dREAL(hyb_left(it,iL)*conjg(r_gz(it,iL)))
            end if


            write(200+iL,'(20(F18.10))')           & 
                 t_grid(it)                        & !1     current time
                 ,nSlab(it,iL)                     & !2     n(t) 
                 ,0.5-x_gz(it,iL)*0.5              & !3     x(t) ---> check for GZ constraint
                 ,n_dot(it,iL)                     & !4     n_dot(i)
                 ,eSlab(it,iL)                     & !5     E(i)
                 ,abs(r_gz(it,iL))**2              & !6     |R(i)|^2
                 ,r_gz(it,iL)                      & !7/8   R(i)
                 ,doublons(it,iL)                  & !9     D=<n_up n_dw>(i)
                 ,hop_plus(it,iL)                  & !10/11 d^+_{i+1}d_i 
                 ,conjg(r_gz(it,iL+1))*r_gz(it,iL) & !12/13 R^+(i+1)*R(i)
                 ,eSlab(it,iL)*abs(r_gz(it,iL))**2 & !14
                 ,dreal(hop_plus(it,iL)*conjg(r_gz(it,iL+1))*r_gz(it,iL)) & !15
                 ,hyb_left(it,iL)                  & !16/17 left hybridization                                    
                 ,hyb_right(it,iL)                             !17/18 right hybridization


            write(500+iL,'(15(F18.10))')           &
                 t_grid(it)                        & !1     current time               
                 ,doublons(it,iL)                  &
                 ,holons(it,iL)                    &
                 ,ph_doublons(it,iL)               &
                 ,ph_holons(it,iL)

         end do

         iL = L

         slab_ene = slab_ene + 2.d0*eSlab(it,iL)*abs(r_gz(it,iL))**2 + U*doublons(it,iL)
         !
         kin_ene = kin_ene + 2.d0*eSlab(it,iL)*abs(r_gz(it,iL))**2
         int_ene = int_ene + U*doublons(it,iL)
         !
         
         hyb_ene  = hyb_ene  + 2.d0*dREAL(hyb_right(it,iL)*conjg(r_gz(it,iL)))


         write(200+iL,'(20(F18.10))')              & 
              t_grid(it)                           & !1     :: current time
              ,nSlab(it,iL)                        & !2     :: n(t) 
              ,0.5-x_gz(it,iL)*0.5d0               & !3     ::  x(t) ---> check for GZ constraint
              ,n_dot(it,iL)                        & !4     ::   n_dot(i)
              ,eSlab(it,iL)                        & !5     ::    E(i)
              ,abs(r_gz(it,iL))**2                 & !6     ::    |R(i)|^2
              ,r_gz(it,iL)                         & !7/8   ::   R(i)
              ,doublons(it,iL)                     & !9     ::  
              ,Z0                                  & !10/11 :: d^+_{i+1,i}*R^+(i+1)*R(i)
              ,Z0                                  & !12/13 :: d^+_{i+1,i}
              ,eSlab(it,iL)*abs(r_gz(it,iL))**2+U*doublons(it,iL) & !14
              ,0.d0                             & !15
              ,hyb_left(it,iL)                     & !16/17 :: left hybridization                                   
              ,hyb_right(it,iL)                            !17/18 :: left hybridization                                    

         write(500+iL,'(15(F18.10))')              &
              t_grid(it)                           & !1     current time               
              ,doublons(it,iL)                  &
              ,holons(it,iL)                    &
              ,ph_doublons(it,iL)               &
              ,ph_holons(it,iL)

         slab_ene = slab_ene/dble(L)
         kin_ene = kin_ene/dble(L)
         int_ene = int_ene/dble(L)
         slab_current = slab_current/dble(L)

         write(10,'(10(F18.10))')                            &
              t_grid(it)                                     &
              ,2.d0*dIMAG(hyb_left(it,1)*conjg(r_gz(it,1)))  &
              ,2.d0*dIMAG(hyb_right(it,L)*conjg(r_gz(it,L))) &
              ,slab_current                                  &
              ,slab_ene                                      &
              ,kin_ene                                       &
              ,int_ene                                       &
              ,hyb_ene
      else  !3d plots

         xSlab = linspace(1.d0,dble(L),L)


         call splot3d("density_3d.out",t_grid,xSlab,nSlab)
         call splot3d("Rgz_3d.out",t_grid,xSlab,r_gz)
         call splot3d("Zgz_3d.out",t_grid,xSlab,ABS(r_gz)**2)
         call splot3d("eSlab.out",t_grid,xSlab,eSlab)

         do iL=1,L
            write(fileout,'(I3.3)') iL
            open(unit=400+iL,file='nk_'//trim(fileout)//'.out')
            do itau = 1,Nt

               if(mod(itau,20).eq.0) then

                  if(.not.dos_plane) then
                     do ik=1,Nk_tot
                        write(400+iL,'(5(F18.10))') vec_k(ik)%x,vec_k(ik)%y,nk(ik,itau,iL)
                     end do
                  else
                     do ik=1,Nk_tot
                        write(400+iL,'(5(F18.10))') ene_dos(ik),nk(ik,itau,iL)
                     end do
                  end if

                  write(400+iL,*)
                  write(400+iL,*)

               end if
            end do
         end do

      end if


    end subroutine print_dynamics



    subroutine initial_condition
      implicit none
      integer     :: ik,ik_sys
      integer     :: ialpha,jalpha
      real(8)     :: kl_j,kl_i,kp
      real(8)     :: ekl_j,ekl_i
      real(8)     :: chem_lead
      integer     :: dimk,idimk
      integer     :: dimkp,dimLead
      integer     :: dimHyb,dimSlab
      integer     :: ilead,jlead
      integer     :: ihyb,islab,jslab
      integer     :: ikp,NpSlab,ikpSlab
      integer     :: igz
      real(8)     :: kpSlab,dkpSlab,occSlab
      real(8)     :: fSlab
      integer     :: itest
      integer     :: it_finer
      integer     :: ik_orth

      complex(8),dimension(:,:),allocatable :: Gleads
      real(8),dimension(:),allocatable :: nSlab
      real(8),dimension(:),allocatable :: eSlab
      type(gz_projector),dimension(:),allocatable :: gz_phi
      real(8),dimension(:),allocatable :: x
      real(8) :: eps_v

      dimk = 2*Nk_orth*(2*Nk_orth+L)+L*L
      dimLead = (2*Nk_orth)
      dimHyb = 2*Nk_orth
      dimSlab = L*L
      allocate(nSlab(L),x(L),eSlab(L))
      nSlab = 0.d0
      do ik=1,Nk_tot
         ek = epsik(ik)
         ik_sys = (ik-1)*(dimk)
         !+----------+!
         !+- Gleads -+!
         !+----------+!
         allocate(Gleads(dimLead,dimLead))
         do ilead = 1,dimLead
            if(ilead.le.Nk_orth) then
               kl_i = k_orth(ilead)
            else
               kl_i = k_orth(ilead-Nk_orth)
            end if
            ekl_i = chain_disp(kl_i)*t_lead
            do jlead=1,dimLead
               if(jlead.le.Nk_orth) then
                  kl_j = k_orth(jlead)
               else
                  kl_j = k_orth(jlead-Nk_orth)
               end if
               ekl_j = chain_disp(kl_j)*t_lead
               ik_sys = ik_sys + 1               
               if(jlead.eq.ilead) then
                  select case(lead_type)
                  case('3d_tb')
                     
                     if(jlead.le.Nk_orth) then
                        psi_t(ik_sys) = -Zi*(1.d0 - fermi(ekl_j+ek,beta_left))
                     else
                        psi_t(ik_sys) = -Zi*(1.d0 - fermi(ekl_j+ek,beta_right))
                     end if

                  case('generic_bath')
                     if(jlead.le.Nk_orth) then
                        psi_t(ik_sys) = -Zi*(1.d0 - fermi(ekl_j,beta_left))
                     else
                        psi_t(ik_sys) = -Zi*(1.d0 - fermi(ekl_j,beta_right))
                     end if
                     ! if(jlead.le.Nk_orth) then
                     !    psi_t(ik_sys) = -Zi*(1.d0 - fermi(ekl_j,beta))
                     ! else
                     !    psi_t(ik_sys) = -Zi*(1.d0 - fermi(ekl_j,beta))
                     ! end if
                  end select
               else
                  psi_t(ik_sys) = Z0
               end if
            end do
         end do
         deallocate(Gleads)
         !+-----------+!
         !+- Ghybrid -+!
         !+-----------+!
         do iSlab = 1,L         
            do ihyb = 1,dimHyb
               ik_sys = ik_sys + 1
               psi_t(ik_sys) = Z0
            end do
         end do
         !+---------+!
         !+- GSlab -+!
         !+---------+!
         itest = 0
         do islab = 1,L
            do jslab = 1,L
               itest = itest +1
               ik_sys = ik_sys + 1
               psi_t(ik_sys) = Gslab_equ(ik,iSlab,jSlab)
            end do
         end do
      end do
      !+-------------------------+!
      !+- GUTZWILLER PROJECTORS -+!
      !+-------------------------+!
      allocate(gz_phi(L))
      itest = 0
      do iSlab = 1,L
         if(abs(eqPhi(iSlab)%p0).gt.init) then
            gz_phi(iSlab)%p0 = eqPhi(iSlab)%p0
            gz_phi(iSlab)%p1 = eqPhi(iSlab)%p1
            gz_phi(iSlab)%p2 = eqPhi(iSlab)%p2
         else
            gz_phi(iSlab)%p0 = init
            gz_phi(iSlab)%p1 = sqrt(1.d0 - 2.d0*init*init)
            gz_phi(iSlab)%p2 = init
         end if

         do igz=1,4
            ik_sys = ik_sys + 1
            select case(igz)
            case(1)
               psi_t(ik_sys) = GZ_double(gz_phi(iSlab))
            case(2)
               psi_t(ik_sys) = 0.d0
            case(3)
               psi_t(ik_sys) = GZ_double(gz_phi(iSlab))
            case(4)
               psi_t(ik_sys) = 0.d0
            end select
         end do
      end do

      !+- initialize useful stuff -+!
      open(unit=20,file='inner_field_profile.out')
      do it_finer=1,2*Nt+1
         do ik_orth=1,Nk_orth
            kp = k_orth(ik_orth)
            vk_L(ik_orth,it_finer) = get_bath_coupling(kp,'L',t_finer(it_finer))
            vk_R(ik_orth,it_finer) = get_bath_coupling(kp,'R',t_finer(it_finer))
         end do
         if(t_finer(it_finer).ge.0.d0) then
            muL(it_finer) = chem_bath('L',t_finer(it_finer))
            muR(it_finer) = chem_bath('R',t_finer(it_finer))
         else
            muL(it_finer) = 0.d0
            muR(it_finer) = 0.d0
         end if
         select case(trim(adjustl(trim(inner_field_type))))
         case("none")
            e_loc(:,it_finer) = 0.d0
         case("linear")
            do iSlab=1,L
               e_loc(iSlab,it_finer) = 0.d0
               if(t_finer(it_finer).ge.time_bias) then
                  e_loc(iSlab,it_finer) = (muL(it_finer) - muR(it_finer))/dble(L+1)/2.d0*dble(L+1-2*iSlab)
               else
                  e_loc(iSlab,it_finer) = 0.d0
               end if
               write(20,*) dble(iSlab),t_finer(it_finer),e_loc(iSlab,it_finer)
            end do
            write(20,*)
         case("linear_asymm")
            eps_v = (muL(it_finer) - muR(it_finer))/dble(L+1)
            do iSlab=1,L
               e_loc(iSlab,it_finer) = 0.d0
               if(t_finer(it_finer).ge.time_bias) then
                  e_loc(iSlab,it_finer) = muL(it_finer) - eps_v - (muL(it_finer)-muR(it_finer)-2.d0*eps_v)/dble(L-1)*dble(iSlab-1)
               else
                  e_loc(iSlab,it_finer) = 0.d0
               end if
               write(20,*) dble(iSlab),t_finer(it_finer),e_loc(iSlab,it_finer)
            end do
            write(20,*)
         case("exp_drop")
            e_loc(:,it_finer) = 0.d0
            do iSlab=1,L/2
               if(t_finer(it_finer).ge.time_bias) then
                  e_loc(iSlab,it_finer) = (muL(it_finer) - muR(it_finer))/dble(L+1)*dble(L-1)/2.d0*exp(-dble(iSlab)/dble(2))   !/(exp(dble(iSlab)-L/4)+1.d0)
                  e_loc(L+1-iSlab,it_finer) = -(muL(it_finer)-muR(it_finer))/dble(L+1)*dble(L-1)/2.d0*exp(-dble(iSlab)/dble(2))!/(exp(dble(iSlab)-L/4)+1.d0)
               else
                  e_loc(iSlab,it_finer) = 0.d0
               end if
            end do
            do iSlab=1,L
               write(20,*) dble(iSlab),t_finer(it_finer),e_loc(iSlab,it_finer)
            end do
            write(20,*)
         end select
         !
         if(.not.lead_chem) then
            muL(it_finer) = 0.d0
            muR(it_finer) = 0.d0
         end if
         !
         do iSlab=1,L
            Uz_time(iSlab,it_finer) = U 
         end do
         !
      end do
      close(20)
    end subroutine initial_condition




  END SUBROUTINE solve_observables_dynamics
  


END MODULE OBS_DYNAMICS
