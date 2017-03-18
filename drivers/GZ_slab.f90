PROGRAM GZ_SLAB_TEMP
  use dynamics
  use slab_equ_temp
  use slab_equ_temp_dop
  use slab_equ
  use global
  implicit none
  real(8),allocatable :: test_dop(:)
  integer :: iL

  call read_input
   
  call reduced_bz(nx_grid,nk_orth,off_set)
  
  if(temp.lt.1.d-3) then
     call gz_equilibrium
  else
     call gz_equilibrium_temp
  end if

  write(*,*) Nt
  if(Nt.gt.5) call solve_dynamics
  
END PROGRAM GZ_SLAB_TEMP
