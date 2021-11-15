PROGRAM GZ_SLAB_TEMP
  use dynamics
  use obs_dynamics
  ! use slab_equ_temp
  ! use slab_equ_temp_dop
  use slab_equ
  use global
  implicit none
  real(8),allocatable :: test_dop(:)
  integer :: ik,jk,i,j,iL,iLL,iz,ik_min,ik_min_
  real(8) :: chem_pulse
  real(8),dimension(2) :: temps
  real(8),dimension(:),allocatable :: Uz_meta
  complex(8),dimension(:,:,:),allocatable :: Gk1,Gk2
  type(gz_projector),dimension(:),allocatable :: phi1,phi2

  integer :: Nsys
  complex(8),dimension(:),allocatable :: psi_left,psi_full
  
  logical :: f1,f2
  integer :: L_left,L0,Lbarrier,Ltot
  integer :: Nbarrier
  real(8) :: U0,Ubarrier
  real(8) :: T_left,T_right
  integer :: Nt_diss,Nt_dyn
  real(8) :: k_diss_qp
  real(8) :: wband
  real(8) :: ekmin
  integer,dimension(:),allocatable :: ik_ord_
  real(8),dimension(:),allocatable :: epsik_ord
  logical :: print_equ

  real(8) :: Efield_
  !
  call parse_input_variable(dos_plane,"DOS",'input_slab.in',default=.false.)
  call parse_input_variable(wband,"wband",'input_slab.in',default=2.d0)
  call parse_input_variable(cut_time,"CUT",'input_slab.in',default=1000.d0)
  call parse_input_variable(pulse_time,"PULSE",'input_slab.in',default=0.d0)
  call parse_input_variable(chem_pulse,"CHEM_PULSE",'input_slab.in',default=chem_shift)
  call parse_input_variable(L_left,"L_left",'input_slab.in',default=10)
  call parse_input_variable(L0,"L_zero",'input_slab.in',default=10)
  call parse_input_variable(Lbarrier,"Lbarrier",'input_slab.in',default=1)
  call parse_input_variable(Nbarrier,"Nbarrier",'input_slab.in',default=5)
  call parse_input_variable(U0,"U_zero",'input_slab.in',default=5.d0)
  call parse_input_variable(Ubarrier,"Ubarrier",'input_slab.in',default=12.d0)
  call parse_input_variable(T_left,"T_left",'input_slab.in',default=0.001d0)
  call parse_input_variable(T_right,"T_right",'input_slab.in',default=0.001d0)
  call parse_input_variable(Nt_diss,"NT_DISS",'input_slab.in',default=5000)
  call parse_input_variable(Nt_dyn,"NT_DYN",'input_slab.in',default=10000)
  call parse_input_variable(k_diss_qp,"K_DISS",'input_slab.in',default=0.01d0)
  call parse_input_variable(print_equ,"print_equ",'input_slab.in',default=.false.)
  !
  call read_input("input_slab.in")
  call save_input_file("input_slab.in")
  !
  !+- all this stuff one day will be replaced with scifor.square_lattice -+!
  if(.not.dos_plane) then
     Nk_tot = nx_grid*(nx_grid+1)/2
     allocate(epsik(Nk_tot))
     call reduced_bz(nx_grid,nk_orth,off_set)
     do ik=1,Nk_tot
        epsik(ik) = square_lattice_disp(vec_k(ik))
     end do
     call build_orth_BZ(Nk_orth)

     allocate(ik_ord(Nk_tot))
     allocate(ik_ord_(Nk_tot+1))
     allocate(epsik_ord(Nk_tot+1))
     ik_ord_=0
     epsik_ord=0.d0
     do ik=1,Nk_tot
        ik_ord(ik)=ik
        ik_ord_(ik)=ik
        epsik_ord(ik)=epsik(ik)
     end do

     ik_min_=0
     i=0
     do ik=1,Nk_tot
        
        ik_min=ik
        ekmin=epsik_ord(ik)
        do jk=ik,Nk_tot
           if(ekmin.gt.epsik_ord(jk)) then
              ekmin=epsik_ord(jk)
              ik_min=jk
           end if           
        end do
        
        ! write(38,*) ik_min,ik_min_,ik,i

        ! if(ik_min.le.ik_min_) then
        !    ik_ord(ik)=ik_min-i+1
        ! else
        !    ik_ord(ik)=ik_min
        ! end if
        ! ik_min_=ik_min        
        !
        epsik_ord(Nk_tot+1) = ekmin       
        if(ik_min.gt.ik) i = i+1
        do jk=1,ik_min-ik
           epsik_ord(ik_min+1-jk)=epsik_ord(ik_min-jk)
        end do
        epsik_ord(ik) = epsik_ord(Nk_tot+1)
        !
        !
        ik_ord_(Nk_tot+1) = ik_ord_(ik_min)
        do jk=1,ik_min-ik
           ik_ord_(ik_min+1-jk)=ik_ord_(ik_min-jk)
        end do
        ik_ord_(ik)=ik_ord_(Nk_tot+1)
        !
        !
     end do
     !

     do ik=1,Nk_tot
        ik_ord(ik)=ik_ord_(ik)
        write(39,*) epsik(ik_ord_(ik)),ik_ord_(ik)
        write(40,*) epsik_ord(ik),epsik(ik)
     end do



!     stop
  else
     Nk_tot = nx_grid
     allocate(epsik(Nk_tot))
     call build_layer_dos(nx_grid,wband=wband)
     if(allocated(wt)) deallocate(wt)
     allocate(wt(Nk_tot))
     wt=wt_dos
     epsik=ene_dos
     call build_orth_BZ(Nk_orth)
  end if
  !write(*,*) t_perp;stop

  Efield_=Efield

  mu_L = 0.d0
  mu_R = 0.d0
  mu_pulseL=0.d0
  mu_pulseR=0.d0

  if(left) then
     mu_L= chem_shift*0.5d0
     mu_pulseL= chem_pulse*0.5d0
  end if
  if(right) then
     mu_R=-chem_shift*0.5d0
     mu_pulseR=-chem_pulse*0.5d0
  end if


  Ltot=L0*(Nbarrier-1)+Lbarrier*Nbarrier
  L=L_left
  if(L.gt.0.and.L.le.Ltot) then
     allocate(Uz(L));

     f1=.false.
     f2=.false.
     do iLL=1,10
        if(f1.and.f2) exit        
        do iL=1,Lbarrier
           iz=(iLL-1)*(Lbarrier+L0)+iL
           if(iz.gt.L) then
              f1=.true.
              exit
           end if
           Uz(iz)=Ubarrier
        end do
        do iL=Lbarrier+1,Lbarrier+1+L0
           iz=(iLL-1)*(Lbarrier+L0)+iL
           if(iz.gt.L) then
              f2=.true.
              exit
           end if
           Uz(iz)=U0
        end do
     end do
     open(unit=20,file='Uprofile_left.out')
     do iL=1,L
        write(20,*) iL,Uz(iL)
     end do
     close(20)
     !
     allocate(Hslab(L,L));  Hslab=0.d0
     do iL=1,L
        Hslab(iL,iL) =1.d0
     end do
     do iL=1,L-1
        Hslab(iL,iL+1) = -t_perp
        Hslab(iL+1,iL) = -t_perp
     end do
     ! 
     temps=T_left
     call gz_equilibrium(temps,Gk1,phi1)
  end if
  allocate(Gslab_equ(Nk_tot,L,L),eqPhi(L))
  Gslab_equ=Gk1
  eqPhi(1:L)%p0=phi1(1:L)%p0
  eqPhi(1:L)%p1=phi1(1:L)%p1
  eqPhi(1:L)%p2=phi1(1:L)%p2
  !
  !+- evolve "Gleft" with strongly damped dynamics -+!
  !
  Efield=0.0d0
  Nt=Nt_diss
  k_diss=0.2
  beta_diss = 1.d0/T_left
  if(Nt.gt.5) then
     call allocate_dynamics(iprint=.false.,Nsys_=Nsys)
     call initial_condition
     allocate(psi_left(Nsys))
     call evolve(iprint=.false.,psi_out=psi_left)
     call deallocate_dynamics
     deallocate(Gslab_equ,eqPhi)
  end if
  !stop
  !+
  call get_GK(psi_left,Gk1)
  !+
  !
  Ltot=L0*(Nbarrier-1)+Lbarrier*Nbarrier
  L=Ltot
  deallocate(Uz)
  deallocate(Hslab)
  if(L.gt.0) then
     allocate(Uz(L));
     f1=.false.
     f2=.false.
     do iLL=1,10
        if(f1.and.f2) exit        
        do iL=1,Lbarrier
           iz=(iLL-1)*(Lbarrier+L0)+iL
           if(iz.gt.L) then
              f1=.true.
              exit
           end if
           Uz(iz)=Ubarrier
        end do
        do iL=Lbarrier+1,Lbarrier+1+L0
           iz=(iLL-1)*(Lbarrier+L0)+iL
           if(iz.gt.L) then
              f2=.true.
              exit
           end if
           Uz(iz)=U0
        end do
     end do
     open(unit=20,file='Uprofile_tot.out')
     do iL=1,L
        write(20,*) iL,Uz(iL)
     end do
     close(20)

     allocate(Hslab(L,L));  Hslab=0.d0
     do iL=1,L
        Hslab(iL,iL) =1.d0
     end do
     do iL=1,L-1
        Hslab(iL,iL+1) = -t_perp
        Hslab(iL+1,iL) = -t_perp
     end do
     ! 
     temps=T_right
     call gz_equilibrium(temps,Gk2,phi2)
  end if  
  allocate(Gslab_equ(Nk_tot,L,L),eqPhi(L))
  Gslab_equ=Gk2
  eqPhi(1:L)%p0=phi2(1:L)%p0
  eqPhi(1:L)%p1=phi2(1:L)%p1
  eqPhi(1:L)%p2=phi2(1:L)%p2
  !+- evolve "Gright" with strongly damped dynamics -+!
  beta_diss = 1.d0/T_right
  Nt=Nt_diss
  k_diss=0.2
  Efield=0.d0
  if(Nt.gt.5) then
     call allocate_dynamics(iprint=print_equ,Nsys_=Nsys)
     call initial_condition
     allocate(psi_full(Nsys))
     call evolve(iprint=print_equ,psi_out=psi_full)
     call deallocate_dynamics
  end if
  if(print_equ) stop
  !stop
  
  !+- copy Gk1 inside psi_full -+!
  if(Efield_.eq.0.d0) then
     call get_Gk(psi_full,Gk2)
     Gk2(:,1:L_left,1:L_left) = Gk1
     call overwrite_Gkl(psi_full,Gk2)
  end if
  
  k_diss=k_diss_qp
  Nt=Nt_dyn
  Efield=Efield_
  if(Nt.gt.5) then
     call allocate_dynamics(iprint=.true.)
     call initial_condition(psi_in=psi_full)
     call evolve(iprint=.true.)
  end if
  !
    
  stop

contains

  

  ! L2=Lmeta*Nmeta+Nmeta-1
  ! L=L2
  ! allocate(Uz(L));
  ! iL=0
  ! do j=1,Nmeta
  !    do i=1,Lmeta
  !       iL=iL+1
  !       Uz(iL) = U2
  !    end do
  !    if(j.lt.Nmeta) then
  !       iL=iL+1
  !       Uz(iL) = Ubarrier
  !    end if
  ! end do
  ! allocate(Uz_meta(L)); Uz_meta=Uz
  ! write(iL),Lmeta,Nmeta
  ! do iL=1,L
  !    write(444,*) iL,Uz(iL)
  ! end do
  
  ! allocate(Hslab(L,L));  Hslab=0.d0
  ! do iL=1,L
  !    Hslab(iL,iL) =1.d0
  ! end do
  ! do iL=1,L-1
  !    Hslab(iL,iL+1) = -t_perp
  !    Hslab(iL+1,iL) = -t_perp
  ! end do
  ! !
  ! temps(1)=T2;  temps(2)=T2n;
  ! call gz_equilibrium(temps,Gk2,phi2)
  ! deallocate(Hslab,Uz)


  ! L=L1+L2
  ! allocate(Gslab_equ(Nk_tot,L,L),eqPhi(L))
  ! if(L1.gt.0) Gslab_equ(:,1:L1,1:L1)=Gk1
  ! Gslab_equ(:,L1+1:L1+L2,L1+1:L1+L2)=Gk2

  ! eqPhi(1:L1)%p0=phi1(:)%p0
  ! eqPhi(1:L1)%p1=phi1(:)%p1
  ! eqPhi(1:L1)%p2=phi1(:)%p2
  ! do iL=1,L2
  !    eqPhi(iL+L1)%p0=phi2(iL)%p0
  !    eqPhi(iL+L1)%p1=phi2(iL)%p1
  !    eqPhi(iL+L1)%p2=phi2(iL)%p2
  ! end do
  ! allocate(Hslab(L,L));  Hslab=0.d0
  ! do iL=1,L
  !    Hslab(iL,iL)=1.d0
  ! end do
  ! do iL=1,L-1
  !    Hslab(iL,iL+1)=-t_perp
  !    Hslab(iL+1,iL)=-t_perp
  ! end do
  ! allocate(Uz(L))
  ! Uz(1:L1)=U1
  ! Uz(L1+1:L)=Uz_meta
  ! do il=1,L
  !    write(445,*) iL,Uz(iL)
  ! end do



  
  ! if(Nt.gt.5) then
  !    call allocate_dynamics(iprint=.true.,Nsys_=Nsys)
  !    ! call initial_condition(psi_in=psi_dyn)
  !    ! allocate(psi_full(Nsys))
  !    ! call evolve(iprint=.true.,psi_out=psi_full)

  !    call initial_condition!(psi_in=psi_dyn)
  !    !allocate(psi_full(Nsys))
  !    call evolve(iprint=.true.)!,psi_out=psi_full)

  ! end if
  ! !
  

  !
END PROGRAM GZ_SLAB_TEMP
