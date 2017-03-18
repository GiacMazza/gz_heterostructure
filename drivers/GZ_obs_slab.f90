PROGRAM GZ_SLAB_TEMP
  use dynamics
  use obs_dynamics
  ! use slab_equ_temp
  ! use slab_equ_temp_dop
  use slab_equ
  use global
  implicit none
  real(8),allocatable :: test_dop(:)
  integer :: iL,ik
  real(8) :: chem_pulse
  
  call read_input

  call parse_input_variable(dos_plane,"DOS",'input_slab.in',default=.false.)
  call parse_input_variable(cut_time,"CUT",'input_slab.in',default=1000.d0)
  call parse_input_variable(pulse_time,"PULSE",'input_slab.in',default=0.d0)
  call parse_input_variable(chem_pulse,"CHEM_PULSE",'input_slab.in',default=chem_shift)
  
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
  allocate(Gslab_equ(Nk_tot,L,L))

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
  
  call gz_equilibrium
  
  write(*,*) Nt
  if(Nt.gt.5) call solve_observables_dynamics

END PROGRAM GZ_SLAB_TEMP
