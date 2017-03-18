!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SLAB AT EQUILIBRIUM 
!Author G Mazza
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
include 'free_energy_subroutines.f90'
MODULE SLAB_EQU_TEMP
  use global
  USE LANCELOT_simple_double

  !##################################################################################
  implicit none
  private
  real(8)                                     :: energy
  real(8),dimension(:),allocatable            :: eZ,eZL,eZR
  complex(8),dimension(:),allocatable         :: eZnn,eZnnL,eZnnR
  complex(8),dimension(:),allocatable         :: hop_plus,hop_minus
  complex(8)                                  :: Delta_1L,Delta_NR

  complex(8),dimension(:),allocatable         :: R
  real(8),dimension(:),allocatable            :: dop
  real(8)                                     :: hop_ene,hubb_ene
  real(8)                                     :: e_layers,e_orth
  real(8),dimension(:),allocatable            :: hopping_layer
  real(8)                                     :: iter_free_energy
  real(8)                                     :: free_energy_old

  public :: GZ_equilibrium_temp

CONTAINS

  !+--------------------------------------------------------------+!
  !   FULL MINIMIZATION OF THE GUTZWILLER FREE ENERGY ESTIMATION   !
  !+--------------------------------------------------------------+!
  SUBROUTINE GZ_equilibrium_temp
    real (8), allocatable, dimension(:) :: x
    real (8), allocatable, dimension(:) :: x_start_ins,x_start_met
    real(8) :: f_start_ins,f_start_met
    integer :: iL,jL,i_theta
    integer :: iter
    character(len=10) :: fileout
    allocate( x(L),x_start_ins(L),x_start_met(L)) 
    allocate(hop_intra(L),delta_plus(L),delta_minus(L))
    
    !+---------------------------------+!
    !+- TWO STEPS CYCLIC MINIMIZATION -+!
    !+---------------------------------+!
    
    !+- choose the initial values for the Gutzwiller parameters -+!
    x_start_ins = 0.001d0
    call minimize_free_energy(x_start_ins)

    x_start_met = pi*0.25
    call minimize_free_energy(x_start_met)

    f_start_ins = free_energy_check(x_start_ins)
    f_start_met = free_energy_check(x_start_met)

    write(*,*) 'ins',f_start_ins
    write(*,*) 'met',f_start_met

    if(f_start_ins.ge.f_start_met) then
       x = x_start_met
       write(*,*) 'metallic like'
    else
       x = x_start_ins
       write(*,*) 'insulating like'
    end if
    !+-----------------------------------------------------------+!

    !+-------------------------------+!
    !+- iterative minimization step -+!
    !+-------------------------------+!
    open(15,file='free_energy_iter.out')
    free_energy_old = 0.d0
    do iter=1,n_max

       write(*,*) '! iteration nr !',iter
       call minimize_free_energy(x)
       !+-------------------------------------------------------------+!

       !+- writeout iteration results -+!
       write(fileout,'(I3.3)') iter
       open(unit=40+iter,file='iter_'//trim(fileout)//'.out')
       do iL=1,L
          write(40+iter,'(10(F18.10))') dble(iL),sin(x(iL))**2,(1.0-cos(x(iL)))/4.0
       end do
       close(40+iter)
       !+------------------------------+!

       !+- check convergence -+!
       iter_free_energy = free_energy_check(x,iter)
       write(15,'(6(F18.10))') dble(iter),iter_free_energy
       if(abs(free_energy_old-iter_free_energy).lt.conv_treshold) then
          write(*,*) 'CONVERGENCE AFTER',iter,'ITERATIONS'
          exit
       else
          free_energy_old = iter_free_energy
       end if
       !+---------------------+!

    end do
    !+---------------------------------+!
    close(15)

  end SUBROUTINE GZ_equilibrium_temp




  !+------------------------------------------------------------+!
  ! function f_test(x)
  !   real(8),dimension(:) :: x
  !   real(8) :: f_test
  !   f_test = (x(1) - 3.d0)**2 + (x(2) + 1.5d0)**2
  ! end function f_test
  !+------------------------------------------------------------+!


  
  subroutine minimize_free_energy(x)
    real (8),dimension(:),intent(inout) :: x        
    !+- GALAHAD ROUTINE VARIABLES -+!
    integer :: n, neq, nin ! n=#of optimization variables; neq=#of equality constraints; ninl=#of inequality constraints
    integer :: exit_code   !final status of the minimization (0 = success; \=0 error code )
    integer :: iters       !(optional) on output the number of iterations performed by the routine
    integer :: maxit       !(optional) max number of iteration performed by the routine; (default = 1000)
    integer :: print_level !(optional) amount of information at each iteration; (default = 1)

    real(8) :: gradtol     !(optional) maximum permitted norm of the projected gradient of the lagrangian function; (default = 1.d-5)
    real(8) :: feastol     !(optional) maximum permitted violation of the constraints; (default = 1.d-5)
    real(8) :: fx          !on output holds the value of function at the optimized parameter

    character(len= 10),allocatable, dimension(:) :: vnames, cnames !names of variables and constraints

    real (8),     allocatable, dimension(:) :: bl                  !lower bound
    real (8),     allocatable, dimension(:) :: bu                  !upper bound

    real (8),     allocatable, dimension(:) :: cx                  !estimate of  constraints 
    real (8),     allocatable, dimension(:) :: y                   !estimate of lagrange multipliers 
    external :: fun,slab_free_energy,free_energy
    !+---------------------------------+!


    complex(8),allocatable,dimension(:,:) :: H_star
    real(8),allocatable,dimension(:) :: w
    real(8),allocatable,dimension(:) :: s_star
    complex(8),allocatable,dimension(:) :: Rhop

    integer :: iL,jL,i_theta

    integer :: ik
    real(8) :: ek,gm,phi
    real(8) :: n_temp
    type(vec2D) :: k
    type(gz_projector),allocatable,dimension(:) :: phi_gz




    n   = size(x)                   ! number of variables
    neq = 0                   ! number of equality constraints, excluding fixed variables
    nin = 0                   ! number of inequality (<= 0) constraints, excluding bounds

    maxit       = 1000
    gradtol     = 0.0000001d0 ! identical to default
    feastol     = 0.0000001d0 ! identical to default
    print_level = 1           ! identical to default


    allocate( bl(n), bu(n), cx(neq+nin), y(neq+nin) )
    allocate( vnames(n),cnames(neq+nin) )

    allocate(H_star(L,L),W(L),Rhop(L),phi_gz(L),s_star(L))

    !+- define boundaries on optimization parameters -+!
    do iL=1,L
       bl(iL) = 0.d0
       bu(iL) = pi*0.5d0
    end do
    !+------------------------------------------------+!

    
    !+- reconstruction of GZ projectors -+!
    do iL=1,L
          i_theta = iL
          gm=pi*0.25d0
          phi=0.d0
          phi_gz(iL)%p0 = sin(x(i_theta)*0.5d0)*sin(gm)*Exp(Zi*phi)
          phi_gz(iL)%p2 = sin(x(i_theta)*0.5d0)*cos(gm)*Exp(Zi*phi)
          phi_gz(iL)%p1 = cos(x(i_theta)*0.5d0)
       end do
       Rhop = GZ_hop(phi_gz)

       ! free_electrons distribution optimization !
       hop_intra    = 0.d0
       delta_plus   = 0.d0
       delta_minus  = 0.d0

       do ik=1,Nk_tot !+- k points loop -+!

          !+- build and diagonalize real space hamiltonian -+!
          k = vec_k(ik)
          ek = square_lattice_disp(k)
          H_star=0.d0
          w = 0.d0
          do iL=1,L
             H_star(iL,iL) =  ek*abs(Rhop(iL))**2
             if(iL.lt.L) then
                H_star(iL,iL+1) = -1.d0*conjg(Rhop(iL))*Rhop(iL+1)
                H_star(iL+1,iL) = -1.d0*Rhop(iL)*conjg(Rhop(iL+1))
             end if
          end do
          if(pbc) then
             H_star(1,L) = -1.d0*conjg(Rhop(1))*Rhop(L)
             H_star(L,1) = -1.d0*conjg(Rhop(L))*Rhop(1)
          end if
          call  matrix_diagonalize(H_star,w,'V','L')
          !+------------------------------------------------+!

          !+- compute observables tracing over free_electrons distribuzion function -+!
          do iL=1,L
             do jL=1,L
                hop_intra(jL) = hop_intra(jL) + 2.d0*ek*abs(H_star(jL,iL))**2*wt(ik)*fermi(w(iL),beta)
                if(jL.lt.L) then
                   delta_plus(jL)  = delta_plus(jL) + 2.d0*conjg(H_star(jL+1,iL))*H_star(jL,iL)*wt(ik)*fermi(w(iL),beta)
                   delta_minus(jL) = delta_minus(jL) + 2.d0*conjg(H_star(jL,iL))*H_star(jL+1,iL)*wt(ik)*fermi(w(iL),beta)
                end if
             end do
          end do
          !+-------------------------------------------------------------------------+!

       end do !+- end k points loop -+! 

       !+- gutzwiller parameter optimization -+!
       print_level=0
       call lancelot_simple( n,  x, fx, exit_code, my_fun = free_energy,&
            bl = bl, bu = bu,                              &
            neq = neq, nin = nin,                          &
            cx = cx, y = y, iters  = iters, maxit = maxit, &
            gradtol = gradtol, feastol = feastol,          &
            print_level = print_level )
       !+-------------------------------------+!



       deallocate( bl, bu, cx, y )
       deallocate( vnames, cnames )
       deallocate(H_star,W,Rhop,phi_gz,s_star)
       

  end subroutine minimize_free_energy





  function free_energy_check(x,iter) result(free)
    !+- input variables -+!
    real(8),dimension(:),intent(in)       :: x
    integer,intent(in),optional           :: iter
    !+-------------------+!

    complex(8),dimension(size(x),size(x)) :: H_star
    real(8),dimension(size(x))            :: w

    type(gz_projector),dimension(size(x)) :: phi_gz
    complex(8),dimension(size(x))         :: Rhop
    real(8),dimension(size(x))            :: Docc

    real(8)                               :: free,energy
    real(8)                               :: ek,gm,phi,n_temp
    real(8)                               :: free_entropy
    real(8)                               :: prj_entropy
    real(8)                               :: entropy
    type(vec2D)                           :: k
    
    integer                               :: Lslab
    integer                               :: iL,jL
    integer                               :: i_theta
    integer                               :: ik

    !+- reconstruction of GZ parameters -+!
    do iL=1,L
       i_theta = iL
       gm=pi*0.25d0
       phi=0.d0
       phi_gz(iL)%p0 = sin(x(i_theta)*0.5d0)*sin(gm)*Exp(Zi*phi)
       phi_gz(iL)%p2 = sin(x(i_theta)*0.5d0)*cos(gm)*Exp(Zi*phi)
       phi_gz(iL)%p1 = cos(x(i_theta)*0.5d0)
    end do
    Rhop = GZ_hop(phi_gz)
    Docc = GZ_double(phi_gz)
    Lslab = size(x)
    !+-----------------------------------+!

    !+- free_electrons free energy -+!
    hop_intra = 0.d0
    delta_plus = 0.d0
    delta_minus= 0.d0
    free_entropy = 0.d0
    free = 0.d0

    do ik=1,Nk_tot
       k = vec_k(ik)
       ek = square_lattice_disp(k)
       H_star=0.d0
       do iL=1,Lslab
          H_star(iL,iL) =  ek*abs(Rhop(iL))**2
          if(iL.lt.Lslab) then
             H_star(iL,iL+1) = -1.d0*conjg(Rhop(iL))*Rhop(iL+1)
             H_star(iL+1,iL) = -1.d0*Rhop(iL)*conjg(Rhop(iL+1))
          end if
       end do
       if(pbc) then
          H_star(1,Lslab) = -1.d0*conjg(Rhop(1))*Rhop(Lslab)
          H_star(Lslab,1) = -1.d0*conjg(Rhop(Lslab))*Rhop(1)
       end if
       call  matrix_diagonalize(H_star,w,'V','L')

       !+- trace over free electrons thermal distribution -+!
       do iL=1,Lslab
          do jL=1,Lslab

             hop_intra(jL) = hop_intra(jL) + 2.d0*ek*abs(H_star(jL,iL))**2*wt(ik)*fermi(w(iL),beta)
             if(jL.lt.Lslab) then
                delta_plus(jL)  = delta_plus(jL) + 2.d0*conjg(H_star(jL+1,iL))*H_star(jL,iL)*wt(ik)*fermi(w(iL),beta)
                delta_minus(jL) = delta_minus(jL) + 2.d0*conjg(H_star(jL,iL))*H_star(jL+1,iL)*wt(ik)*fermi(w(iL),beta)
             end if

             !+- free electrons entropy -+!
             if(jL.eq.1) then
                n_temp = fermi(w(iL),beta)
                if(n_temp.gt.1.d-10) then
                   free_entropy = free_entropy - ( n_temp*log(n_temp))*wt(ik)*2.d0 
                end if
                if(abs(1.d0-n_temp).gt.1.d-10) then
                   free_entropy = free_entropy - (1.d0-n_temp)*log(1.d0-n_temp )*wt(ik)*2.d0
                end if
             end if
             !+--------------------------+!
          end do
       end do
    end do !k-sum


    prj_entropy = 0.d0
    do iL=1,L
       !+- energy contribution to free energy -+!
       free = free + hop_intra(iL)*ABS(Rhop(iL))**2
       if(iL.gt.1) then
          free = free - conjg(Rhop(iL-1))*Rhop(iL)*delta_minus(iL-1)
       end if
       if(iL.lt.Lslab) then
          free = free - conjg(Rhop(iL+1))*Rhop(iL)*delta_plus(iL)
       end if
       free = free + Uz(iL)*Docc(iL)
       !+--------------------------------------+!

       !+- projector entropy -+!
       if(abs(phi_gz(iL)%p0).gt.1.d-10) then
          prj_entropy = prj_entropy - abs(phi_gz(iL)%p0)**2*log(4.d0*phi_gz(iL)%p0) 
       end if

       if(abs(phi_gz(iL)%p2).gt.1.d-10) then
          prj_entropy = prj_entropy - abs(phi_gz(iL)%p2)**2*log(4.d0*phi_gz(iL)%p2) 
       end if

       if(abs(phi_gz(iL)%p1).gt.1.d-10) then
          prj_entropy = prj_entropy - abs(phi_gz(iL)%p1)**2*log(2.d0*phi_gz(iL)%p1) 
       end if
       !+---------------------+!
    end do

    !+- internal energy -+!
    energy = free

    !+- free energy estimation -+!
    free = free -temp*(prj_entropy + free_entropy)
    
    if(present(iter)) then
       write(15,'(6(F18.10))') dble(iter),free,energy,free_entropy,prj_entropy,free_entropy+prj_entropy
    end if

  end function free_energy_check

END MODULE SLAB_EQU_TEMP




!#####################################################################################################!
!#####################################################################################################!
!#####################################################################################################!
!#####################################################################################################!
!#####################################################################################################!
!#####################################################################################################!
!#####################################################################################################!
!#####################################################################################################!
!#####################################################################################################!
!#####################################################################################################!
!#####################################################################################################!
!#####################################################################################################!
!#####################################################################################################!
!#####################################################################################################!
!#####################################################################################################!
!#####################################################################################################!



