PROGRAM GZ_SLAB_TEMP
  use dynamics
  use slab_equ_temp
  use slab_equ_temp_dop
  use slab_equ
  use global
  !use PARSE_INPUT
  implicit none
  real(8),allocatable :: test_dop(:)
  integer :: iL
  real(8) :: Uins
  real(8) :: chem_pulse

  call read_input
  call parse_input_variable(Uins,"Uins",'input_slab.in',default=16.005d0)
  call parse_input_variable(chem_equ,"CHEM_EQU",'input_slab.in',default=0.d0)  
  call parse_input_variable(cut_time,"CUT",'input_slab.in',default=1000.d0)
  call parse_input_variable(pulse_time,"PULSE",'input_slab.in',default=0.d0)

  call parse_input_variable(chem_pulse,"CHEM_PULSE",'input_slab.in',default=chem_shift)

  mu_L = 0.d0
  mu_R = 0.d0
  mu_pulseL=0.d0
  mu_pulseR=0.d0
    
  if(left) then
     mu_L= chem_shift*0.5d0 + chem_equ
     mu_pulseL= chem_pulse*0.5d0 + chem_equ
  end if
  if(right) then
     mu_R=-chem_shift*0.5d0 + chem_equ
     mu_pulseR=-chem_pulse*0.5d0 + chem_equ
  end if
  

  call reduced_bz(nx_grid,nk_orth,off_set)
  call gz_equilibrium
  
  
  write(*,*) Nt
  if(Nt.gt.5) call solve_dynamics
  
END PROGRAM GZ_SLAB_TEMP
