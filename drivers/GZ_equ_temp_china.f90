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
     allocate(local_field(L))
     local_field = 0.d0
     call gz_equilibrium_temp_china
     deallocate(local_field)
  end if
  
  ! allocate(test_dop(L))
  ! test_dop = dop_layer
  ! ! do iL=1,L
  ! !    test_dop(iL) = dop_layer*exp(-abs(dble(iL)-1.d0)*2.d0)
  ! ! end do
  ! test_dop = 0.d0
  ! test_dop(1) = dop_layer
  ! test_dop(2) = dop_layer*0.5d0

  ! call equ_fixed_doping(test_dop)
  ! deallocate(test_dop)
  write(*,*) Nt
  if(Nt.gt.5) call solve_dynamics
  
END PROGRAM GZ_SLAB_TEMP
