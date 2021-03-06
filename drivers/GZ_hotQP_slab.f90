PROGRAM GZ_SLAB_TEMP
  use dynamics
  use obs_dynamics
  ! use slab_equ_temp
  ! use slab_equ_temp_dop
  use slab_equ
  use global
  implicit none
  real(8),allocatable :: test_dop(:)
  integer :: iL,ik,L1,L2
  real(8) :: chem_pulse,U1,U2,T1,T2,T1n,T2n
  real(8),dimension(2) :: temps
  complex(8),dimension(:,:,:),allocatable :: Gk1,Gk2
  type(gz_projector),dimension(:),allocatable :: phi1,phi2
  !
  !call read_input

  call parse_input_variable(dos_plane,"DOS",'input_slab.in',default=.false.)
  call parse_input_variable(cut_time,"CUT",'input_slab.in',default=1000.d0)
  call parse_input_variable(pulse_time,"PULSE",'input_slab.in',default=0.d0)
  call parse_input_variable(chem_pulse,"CHEM_PULSE",'input_slab.in',default=chem_shift)
  call parse_input_variable(L1,"L_left",'input_slab.in',default=10)
  call parse_input_variable(L2,"L_right",'input_slab.in',default=10)
  call parse_input_variable(U1,"U_left",'input_slab.in',default=5.d0)
  call parse_input_variable(U2,"U_right",'input_slab.in',default=5.d0)
  call parse_input_variable(T1,"T_left",'input_slab.in',default=0.0001d0)
  call parse_input_variable(T2,"T_right",'input_slab.in',default=0.0001d0)
  call parse_input_variable(T1n,"Tn_left",'input_slab.in',default=0.0001d0)
  call parse_input_variable(T2n,"Tn_right",'input_slab.in',default=0.0001d0)
  !
  call read_input("input_slab.in")
  call save_input_file("input_slab.in")
  
  !+- all this stuff one day will be replaced with scifor.square_lattice -+!
  if(.not.dos_plane) then
     Nk_tot = nx_grid*(nx_grid+1)/2
     allocate(epsik(Nk_tot))
     call reduced_bz(nx_grid,nk_orth,off_set)
     do ik=1,Nk_tot
        epsik(ik) = square_lattice_disp(vec_k(ik))
     end do
     call build_orth_BZ(Nk_orth)
  else
     Nk_tot = nx_grid
     allocate(epsik(Nk_tot))
     call build_layer_dos(nx_grid)
     if(allocated(wt)) deallocate(wt)
     allocate(wt(Nk_tot))
     wt=wt_dos
     epsik=ene_dos
     call build_orth_BZ(Nk_orth)
  end if


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
  !
  L = L1 + L2
  !
  allocate(Uz(L));Uz=U1
  allocate(Hslab(L,L));  Hslab=0.d0
  do iL=1,L
     Hslab(iL,iL) =1.d0
  end do  
  do iL=1,L-1
     Hslab(iL,iL+1) = -t_perp
     Hslab(iL+1,iL) = -t_perp
  end do  
  !
  temps(1)=T1;  temps(2)=T1n;
  call gz_equilibrium(temps,Gk1,phi1,L1)
  allocate(Gslab_equ(Nk_tot,L,L),eqPhi(L))
  Gslab_equ=Gk1
  eqPhi(1:L)%p0=phi1(:)%p0
  eqPhi(1:L)%p1=phi1(:)%p1
  eqPhi(1:L)%p2=phi1(:)%p2
  !
  write(*,*) Nt
  if(Nt.gt.5) call solve_observables_dynamics
  !
END PROGRAM GZ_SLAB_TEMP
