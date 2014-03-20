  !+-----------------------------------------------------------------------------------------+!
  !+- SLAB AT EQUILIBRIUM USING GUTZWILLER ESTIMATION OF THE FINITE TEMPERATURE FREE ENERGY -+!
  !+- Author: G. Mazza @ SISSA                                                              -+!
  !+-----------------------------------------------------------------------------------------+!
  MODULE SLAB_EQU_TEMP_DOP
!    USE MIN_AMOEBA
    USE GLOBAL
    USE LANCELOT_simple_double
    implicit none
    private

    integer :: amoeba_funct_calls
    logical :: print_iter
    real(8),allocatable :: dop_ext(:)
    
    public  :: equ_fixed_doping
    
  CONTAINS

    SUBROUTINE EQU_FIXED_DOPING(dop_profile)
      real(8),intent(in)          :: dop_profile(:)
      real(8) :: f_fixed_doping

      !+- inner variables -+!
      real(8),dimension(:),allocatable :: opt_gz,gz_met,gz_ins,dop_gz
      complex(8),dimension(:),allocatable :: R_gz

      real(8) :: iter_free_energy
      real(8) :: free_energy_old
      real(8) :: free_met,free_init
      real(8) :: free_ins
      real(8) :: theta,gm,phi


      integer :: Lsize
      integer :: n_opt_gz,iL
      integer :: out_flag
      integer :: iter

      type(gz_projector),dimension(:),allocatable :: phi_met,phi_ins
      real(8),dimension(:),allocatable :: mu_met,mu_ins

      !+- GALAHAD ROUTINE VARIABLES -+!
      integer                                      :: n_min, neq, nin  ! n=#of optimization variables; neq=#of equality constraints; ninl=#of inequality constraints
      integer                                      :: exit_code        !final status of the minimization (0 = success; \=0 error code )
      integer                                      :: iters            !(optional) on output the number of iterations performed by the routine
      integer                                      :: maxit            !(optional) max number of iteration performed by the routine; (default = 1000)
      integer                                      :: print_level      !(optional) amount of information at each iteration; (default = 1)

      real(8)                                      :: gradtol          !(optional) maximum permitted norm of the projected gradient of the lagrangian function; (default = 1.d-5)
      real(8)                                      :: feastol          !(optional) maximum permitted violation of the constraints; (default = 1.d-5)
      real(8)                                      :: fx               !on output holds the value of function at the optimized parameter
      character(len= 10),allocatable, dimension(:) :: vnames, cnames   !names of variables and constraints

      real (8),     allocatable, dimension(:)      :: bl               !lower bound
      real (8),     allocatable, dimension(:)      :: bu               !upper bound

      real (8),     allocatable, dimension(:)      :: cx               !estimate of  constraints 
      real (8),     allocatable, dimension(:)      :: y                !estimate of lagrange multipliers 

      
      external :: free_energy_GZ_parameters,fun

      interface
         subroutine minimize_free_energy_fixed_doping(x,x_half,out_flag)
           real(8),dimension(:),intent(inout) :: x        
           real(8),dimension(:),intent(in)    :: x_half
           integer                            :: out_flag
         end subroutine minimize_free_energy_fixed_doping
         subroutine free_energy_check(x,free) 
           use global
           real(8),dimension(:),intent(in)       :: x
           real(8) :: free
         end subroutine free_energy_check
      end interface

      n_opt_gz = 3*L
      allocate(opt_gz(n_opt_gz),phi_met(L),phi_ins(L),dop_gz(L),R_gz(L))

      allocate(mu_star(L),hop_intra(L),delta_plus(L),delta_minus(L),tmp_dop(L),phi_out(L))
      allocate(local_field(L),mu_met(L),mu_ins(L))
      
      local_field=0.d0

      !+- MINIMIZATION FROM METALLIC-LIKE INITIAL VALUES -+!
      do iL=1,L
         opt_gz((iL-1)*3+1) = 0.5d0 + dop_profile(iL)*0.5d0
         opt_gz((iL-1)*3+2) = 0.5d0 - dop_profile(iL)*0.5d0
         opt_gz((iL-1)*3+3) = sqrt( 1.d0 - dop_profile(iL)**2.d0 )/sqrt(2.d0)
      end do
      !+-------------------------------------------------+!

      call free_energy_check(opt_gz,free_init)

      tmp_dop = dop_profile
      
      !+- lancelot_routine parameters -+!
      n_min = 3*L
      neq   = 2*L                   
      nin   = 0                     
      maxit       = 10000
      gradtol     = 0.0000001d0 
      feastol     = 0.0000001d0 
      print_level = 1
      !+-------------------------------+!   

      allocate( bl(n_min), bu(n_min), cx(neq+nin), y(neq+nin) )
      allocate( vnames(n_min),cnames(neq+nin) )
      !+- boundaries on optimization parameters -+!
      do iL=1,L
         ! theta boundaries
         bl((iL-1)*3+1) = tmp_dop(iL)
         bu((iL-1)*3+1) = 0.5d0 + tmp_dop(iL)/2.d0

         bl((iL-1)*3+2) = -tmp_dop(iL)
         bu((iL-1)*3+2) = 0.5d0 - tmp_dop(iL)/2.d0
         
         bl((iL-1)*3+3) = sqrt(0.5d0 - tmp_dop(iL)**2.d0/2.d0) 
         bu((iL-1)*3+3) = sqrt(1.d0 - tmp_dop(iL)**2.d0/2.d0) 
      end do
      !+-----------------------------------------+!

      open(11,file='min_met.out')
      !+- gutzwiller parameter optimization -+!
      call lancelot_simple( n_min,  opt_gz, fx, exit_code, my_fun = free_energy_GZ_parameters, &
           bl = bl, bu = bu,                                                              &
           neq = neq, nin = nin,                                                          &
           cx = cx, y = y, iters  = iters, maxit = maxit,                                 &
           gradtol = gradtol, feastol = feastol,                                          &
           print_level = print_level )
      !+-------------------------------------+!
      close(11)
      
      mu_met = mu_star
      mu_star = 0.d0

      call free_energy_check(opt_gz,free_met)      
      open(12,file='out_met.out')
      do iL=1,L
         phi_met(iL)%p0 = opt_gz((iL-1)*3+1)
         phi_met(iL)%p2 = opt_gz((iL-1)*3+2)
         phi_met(iL)%p1 = opt_gz((iL-1)*3+3)
         R_gz(iL) = GZ_hop(phi_met(iL))
         write(*,'(8f18.10)') dble(iL),opt_gz((iL-1)*3+1),opt_gz((iL-1)*3+2),opt_gz((iL-1)*3+3),R_gz(iL),GZ_doping(phi_met(iL))
         write(12,'(8f18.10)') dble(iL),R_gz(iL),GZ_doping(phi_met(iL)),opt_gz((iL-1)*3+1),opt_gz((iL-1)*3+2),opt_gz((iL-1)*3+3)
      end do
      close(12)

      open(11,file='min_ins.out')
      !+- MINIMIZATION FROM INSULATING-LIKE INITIAL VALUES -+!
      do iL=1,L
         opt_gz((iL-1)*3+1) = 0.25 + dop_profile(iL)*0.5d0
         opt_gz((iL-1)*3+2) = 0.25 -  dop_profile(iL)*0.5d0
         opt_gz((iL-1)*3+3) = sqrt(1.d0 - opt_gz((iL-1)*3+1)**2.d0 - opt_gz((iL-1)*3+2)**2.d0 )!sqrt( 1.99d0 - dop_profile(iL)**2.d0 )/sqrt(2.d0)
      end do
      call lancelot_simple( n_min,  opt_gz, fx, exit_code, my_fun = free_energy_GZ_parameters, &
           bl = bl, bu = bu,                                                              &
           neq = neq, nin = nin,                                                          &
           cx = cx, y = y, iters  = iters, maxit = maxit,                                 &
           gradtol = gradtol, feastol = feastol,                                          &
           print_level = print_level )
      !+----------------------------------------------------+!
      close(11)
      
      mu_ins = mu_star
      open(12,file='out_ins.out')
      do iL=1,L
         phi_ins(iL)%p0 = opt_gz((iL-1)*3+1)
         phi_ins(iL)%p2 = opt_gz((iL-1)*3+2)
         phi_ins(iL)%p1 = opt_gz((iL-1)*3+3)
         R_gz(iL) = GZ_hop(phi_ins(iL))
         write(*,'(8f18.10)') dble(iL),opt_gz((iL-1)*3+1),opt_gz((iL-1)*3+2),opt_gz((iL-1)*3+3),R_gz(iL),GZ_doping(phi_ins(iL))
         write(12,'(8f18.10)') dble(iL),R_gz(iL),GZ_doping(phi_ins(iL)),opt_gz((iL-1)*3+1),opt_gz((iL-1)*3+2),opt_gz((iL-1)*3+3)
      end do
      close(12)
      call free_energy_check(opt_gz,free_ins)      

      !+- compare solutions -+!
      if(free_ins.ge.free_met) then
         write(*,*) 'METALLIC-LIKE SOLUTION',free_ins - free_met
         mu_star = mu_met
         call get_equ_temp_conditions(phi_met)
      else
         mu_star = mu_ins
         call get_equ_temp_conditions(phi_ins)
         write(*,*) 'INSULATING-LIKE SOLUTION',free_met - free_ins
      end if


      write(*,*) 'met',free_met,free_init
      write(*,*) 'ins',free_ins,free_init

      deallocate(opt_gz,phi_met,phi_ins,dop_gz,R_gz)

    END SUBROUTINE EQU_FIXED_DOPING
    







    SUBROUTINE get_equ_temp_conditions(phi_optimized)
      type(gz_projector),dimension(:)  :: phi_optimized
      complex(8),dimension(L,L) :: H_star

      complex(8),dimension(L)   :: R_opt
      real(8),dimension(L)      :: w,dop_opt
      
      integer                   :: iL,jL,iE,ik
      real(8)                   :: theta,gamma,phi,ek
      type(vec2D)               :: k

      eqPhi = phi_optimized
      R_opt = GZ_hop(eqPhi)
      dop_opt = GZ_doping(eqPhi)
      
      open(20,file='out_equ_temp.out')
      do iL=1,L
         write(20,*) iL,abs(R_opt(iL))**2.d0,dop_opt(iL)
      end do
      close(20)

      !+- get uncorrelated density matrix -+!
      do ik=1,Nk_tot !+- k points loop -+!
         k = vec_k(ik)
         ek = square_lattice_disp(k)
         H_star=0.d0
         w = 0.d0
         do iL=1,L
            H_star(iL,iL) =  ek*abs(R_opt(iL))**2 + mu_star(iL)
            if(iL.lt.L) then
               H_star(iL,iL+1) = -1.d0*conjg(R_opt(iL))*R_opt(iL+1)
               H_star(iL+1,iL) = -1.d0*R_opt(iL)*conjg(R_opt(iL+1))
            end if
         end do
         if(pbc) then
            H_star(1,L) = -1.d0*conjg(R_opt(1))*R_opt(L)
            H_star(L,1) = -1.d0*conjg(R_opt(L))*R_opt(1)
         end if
         call  matrix_diagonalize(H_star,w,'V','L')
         do iL=1,L
            do jL=1,L
               Gslab_equ(ik,iL,jL) = 0.d0
               do iE = 1,L
                  Gslab_equ(ik,iL,jL) = Gslab_equ(ik,iL,jL) + Zi*conjg(H_star(jL,iE))*H_star(iL,iE)*fermi(w(iE),beta)
               end do
               if(iL.eq.jL) then
                  Gslab_equ(ik,iL,jL) = Gslab_equ(ik,iL,jL) - Zi
               end if
            end do
         end do
      end do !+- end k points loop -+! 
      
    END SUBROUTINE get_equ_temp_conditions

  END MODULE SLAB_EQU_TEMP_DOP












  subroutine free_energy_check( x, f )
  use global
  !use slab_equ
  implicit none
  !+- routine variables -+!
  real(8), intent( in )           :: x( : )
  real(8), intent( out )          :: f

  !+- inner variables -+!
  integer                         :: iL,jL
  type(gz_projector),dimension(L) :: phi_gz
  complex(8),dimension(L)         :: Rhop
  real(8),dimension(L)            :: Docc,xgz


  !+- slater determinant optimization -+!
  complex(8),dimension(L,L)       :: H_star
  real(8),dimension(L)            :: w
  !character(len=1)        :: jobz
  !character(len=1)        :: uplo

  type(vec2D) :: k
  integer                         :: ik,lx,ly
  real(8)  :: ek
  integer  :: i_theta,i_gm,i_phi

  real(8) :: theta,gm,phi
  real(8) :: e_hop,e_hubb,entropy
  !real(8),dimension(L) :: s_star
  real(8) :: s_star_slab
  real(8) :: p0,p1,p2
  real(8)                                     :: n_temp
  real(8)  :: free_entropy
  real(8),dimension(:),allocatable :: coulomb
  
  do iL=1,L
     phi_gz(iL)%p0 = x((iL-1)*3 + 1)
     phi_gz(iL)%p2 = x((iL-1)*3 + 2)
     phi_gz(iL)%p1 = x((iL-1)*3 + 3)
  end do

  allocate(coulomb(L))
     !+- FREE ENERGY ESTIMATION WITHIN GZ APPROX -+!
  

     !+- GZ parameter reconstruction -+!
     Rhop = GZ_hop(phi_gz)
     Docc = GZ_double(phi_gz)
     xgz  = GZ_doping(phi_gz) 
     !+-------------------------------+!


     !+- minimize F_star with respect to layer dependent chemical potentials -+!
     mu_star = 0.d0  
     !call fzero_broyden(nstar,mu_star)
     call fsolve(nstar,mu_star,tol=1.d-8)
     do iL=1,L
        write(10,'(10(F18.10))') dble(iL),mu_star(iL),tmp_dop(iL),abs(Rhop(iL))
     end do
     !+-----------------------------------------------------------------------+!

     !+- trace over uncorrelated distribution -+!
     hop_intra    = 0.d0
     delta_plus   = 0.d0
     delta_minus  = 0.d0
     free_entropy = 0.d0
     do ik=1,Nk_tot !+- k points loop -+!
        !+- build and diagonalize real space hamiltonian -+!
        k = vec_k(ik)
        ek = square_lattice_disp(k)
        H_star=0.d0
        w = 0.d0
        do iL=1,L
           H_star(iL,iL) =  ek*abs(Rhop(iL))**2 + mu_star(iL)
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
        do iL=1,L
           do jL=1,L
              hop_intra(jL) = hop_intra(jL) + 2.d0*ek*abs(H_star(jL,iL))**2*wt(ik)*fermi(w(iL),beta)
              if(jL.lt.L) then
                 delta_plus(jL)  = delta_plus(jL) + 2.d0*conjg(H_star(jL+1,iL))*H_star(jL,iL)*wt(ik)*fermi(w(iL),beta)
                 delta_minus(jL) = delta_minus(jL) + 2.d0*conjg(H_star(jL,iL))*H_star(jL+1,iL)*wt(ik)*fermi(w(iL),beta)
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
              end if
           end do
        end do
     end do !+- end k points loop -+! 

     coulomb = 0.d0
     if(dabs(alpha_electrostatic).gt.1.d-5) then
        do lx = 1,L
           coulomb(lx) = 0.d0
           do ly = 1,L
              coulomb(lx) = coulomb(lx) + alpha_electrostatic*abs(dble(lx-ly))*x(ly)
           end do
        end do
    end if


     e_hop  = 0.d0
     e_hubb = 0.d0
     entropy = 0.d0
     do iL=1,L
        !+- hopping energy -+!
        e_hop = e_hop + hop_intra(iL)*abs(Rhop(iL))**2
        if (iL.lt.L) then
           e_hop = e_hop - conjg(Rhop(iL+1))*Rhop(iL)*delta_plus(iL)
        end if
        if(iL.gt.1) then
           e_hop = e_hop - conjg(Rhop(iL-1))*Rhop(iL)*delta_minus(iL-1)
        end if
        !+------------------+!

        !+- hubbard energy -+!
        e_hubb = e_hubb + Uz(iL)*Docc(iL)  + local_field(iL)*(1.d0 - xgz(iL)) + coulomb(iL)*xgz(iL)
        !+------------------+!


        p0 = 0.25d0*(1.d0+tmp_dop(iL))**2.d0
        p1 = 0.5d0*(1.d0-tmp_dop(iL)**2.d0)
        p2 = 0.25d0*(1.d0-tmp_dop(iL))**2.d0

        !+- projectors entropy contribution -+!
        if(abs(phi_gz(iL)%p0).gt.1.d-8) then
           entropy = entropy - abs(phi_gz(iL)%p0)**2*log(abs(phi_gz(iL)%p0)**2/p0)
        end if
        if(abs(phi_gz(iL)%p2).gt.1.d-8) then
           entropy = entropy - abs(phi_gz(iL)%p2)**2*log(abs(phi_gz(iL)%p2)**2/p2)
        end if
        if(abs(phi_gz(iL)%p1).gt.1.d-8) then
           entropy = entropy - abs(phi_gz(iL)%p1)**2*log(abs(phi_gz(iL)%p1)**2/p1)
        end if
        !+-----------------------------------+!
     end do

     !+- free energy -+!
     f = e_hop + e_hubb - Temp*(entropy + free_entropy)
     !+---------------+!



contains



  FUNCTION NSTAR(mu)
    USE global
    real(8),dimension(:) :: mu
    real(8),dimension(size(mu)) :: nstar


    real(8) :: fstar_lagrange
    real(8) :: entropy,energy
    real(8) :: n_temp

    hop_intra   = 0.d0
    delta_plus  = 0.d0
    delta_minus = 0.d0
    nstar = 0.d0
    do ik=1,Nk_tot !+- k points loop -+!
       k = vec_k(ik)
       ek = square_lattice_disp(k)
       H_star=0.d0
       w = 0.d0
       do iL=1,L
          H_star(iL,iL) =  ek*abs(Rhop(iL))**2 + mu(iL)
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
       do iL = 1,L
          do jL=1,L
             nstar(jL) = nstar(jL) + abs(H_star(jL,iL))**2*fermi(w(iL),beta)*wt(ik)*2.d0
             hop_intra(jL) = hop_intra(jL) + 2.d0*ek*abs(H_star(jL,iL))**2*wt(ik)*fermi(w(iL),beta)
             if(jL.lt.L) then
                delta_plus(jL)  = delta_plus(jL) + 2.d0*conjg(H_star(jL+1,iL))*H_star(jL,iL)*wt(ik)*fermi(w(iL),beta)
                delta_minus(jL) = delta_minus(jL) + 2.d0*conjg(H_star(jL,iL))*H_star(jL+1,iL)*wt(ik)*fermi(w(iL),beta)
             end if
          end do
       end do
    end do !+- end k points loop -+! 

    nstar = nstar - (1.d0 - tmp_dop)

    return
  END FUNCTION NSTAR


end subroutine free_energy_check
















subroutine free_energy_GZ_parameters( x, f, i )
  use global
  !use slab_equ
  implicit none
  !+- routine variables -+!
  real(8), intent( in )           :: x( : )
  real(8), intent( out )          :: f
  integer, intent( in ), optional :: i

  !+- inner variables -+!
  integer                         :: iL,jL
  type(gz_projector),dimension(L) :: phi_gz
  complex(8),dimension(:),allocatable         :: Rhop
  real(8),dimension(L)            :: Docc,xgz


  !+- slater determinant optimization -+!
  complex(8),dimension(L,L)       :: H_star
  real(8),dimension(L)            :: w
  !character(len=1)        :: jobz
  !character(len=1)        :: uplo

  type(vec2D) :: k
  integer                         :: ik,lx,ly
  real(8)  :: ek
  integer  :: i_theta,i_gm,i_phi

  real(8) :: theta,gm,phi
  real(8) :: e_hop,e_hubb,entropy
  !real(8),dimension(L) :: s_star
  real(8) :: s_star_slab
  real(8) :: p0,p1,p2
  real(8)                                     :: n_temp
  real(8)  :: free_entropy
  real(8),dimension(:),allocatable :: coulomb
  
  do iL=1,L
     phi_gz(iL)%p0 = x((iL-1)*3 + 1)
     phi_gz(iL)%p2 = x((iL-1)*3 + 2)
     phi_gz(iL)%p1 = x((iL-1)*3 + 3)
  end do


  allocate(Rhop(L),coulomb(L))

  if ( .not. present( i ) ) then
     !+- FREE ENERGY ESTIMATION WITHIN GZ APPROX -+!


     !+- GZ parameter reconstruction -+!
     Rhop = GZ_hop(phi_gz)
     Docc = GZ_double(phi_gz)
     xgz  = GZ_doping(phi_gz) 
     !+-------------------------------+!

     !+- minimize F_star with respect to layer dependent chemical potentials -+!
     mu_star = 0.d0  
     call fsolve(nstar,mu_star,tol=1.d-10)
     do iL=1,L
        write(10,'(10(F18.10))') dble(iL),mu_star(iL),tmp_dop(iL),abs(Rhop(iL))
     end do
     !+-----------------------------------------------------------------------+!

     !+- trace over uncorrelated distribution -+!
     hop_intra    = 0.d0
     delta_plus   = 0.d0
     delta_minus  = 0.d0
     free_entropy = 0.d0
     do ik=1,Nk_tot !+- k points loop -+!
        !+- build and diagonalize real space hamiltonian -+!
        k = vec_k(ik)
        ek = square_lattice_disp(k)
        H_star=0.d0
        w = 0.d0
        do iL=1,L
           H_star(iL,iL) =  ek*abs(Rhop(iL))**2 + mu_star(iL) 
           if(iL.lt.L) then
              H_star(iL,iL+1) = -1.d0*conjg(Rhop(iL))*Rhop(iL+1)
              H_star(iL+1,iL) = -1.d0*Rhop(iL)*conjg(Rhop(iL+1))
           end if
        end do
        if(pbc) then
           H_star(1,L) = -1.d0*conjg(Rhop(1))*Rhop(L)
           H_star(L,1) = -1.d0*conjg(Rhop(L))*Rhop(1)
        end if
        !write(*,'(10f18.10)') Rhop(1),Rhop(2),Rhop(3),Rhop(20)
        call  matrix_diagonalize(H_star,w,'V','L')
        !+------------------------------------------------+!
        do iL=1,L
           do jL=1,L
              hop_intra(jL) = hop_intra(jL) + 2.d0*ek*abs(H_star(jL,iL))**2*wt(ik)*fermi(w(iL),beta)
              if(jL.lt.L) then
                 delta_plus(jL)  = delta_plus(jL) + 2.d0*conjg(H_star(jL+1,iL))*H_star(jL,iL)*wt(ik)*fermi(w(iL),beta)
                 delta_minus(jL) = delta_minus(jL) + 2.d0*conjg(H_star(jL,iL))*H_star(jL+1,iL)*wt(ik)*fermi(w(iL),beta)
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
              end if
           end do
        end do
     end do !+- end k points loop -+! 


     coulomb = 0.d0
     if(dabs(alpha_electrostatic).gt.1.d-5) then
        do lx = 1,L
           coulomb(lx) = 0.d0
           do ly = 1,L
              coulomb(lx) = coulomb(lx) + alpha_electrostatic*abs(dble(lx-ly))*x(ly)
           end do
        end do
    end if


     e_hop  = 0.d0
     e_hubb = 0.d0
     entropy = 0.d0
     do iL=1,L
        !+- hopping energy -+!
        e_hop = e_hop + hop_intra(iL)*abs(Rhop(iL))**2
        if (iL.lt.L) then
           e_hop = e_hop - conjg(Rhop(iL+1))*Rhop(iL)*delta_plus(iL)
        end if
        if(iL.gt.1) then
           e_hop = e_hop - conjg(Rhop(iL-1))*Rhop(iL)*delta_minus(iL-1)
        end if
        !+------------------+!

        !+- hubbard energy -+!
        e_hubb = e_hubb + Uz(iL)*Docc(iL)  + local_field(iL)*(1.d0 - xgz(iL)) + coulomb(iL)*xgz(iL)
        !+------------------+!


        p0 = 0.25d0*(1.d0+tmp_dop(iL))**2.d0
        p1 = 0.5d0*(1.d0-tmp_dop(iL)**2.d0)
        p2 = 0.25d0*(1.d0-tmp_dop(iL))**2.d0

        !+- projectors entropy contribution -+!
        if(abs(phi_gz(iL)%p0).gt.1.d-8) then
           entropy = entropy - abs(phi_gz(iL)%p0)**2*log(abs(phi_gz(iL)%p0)**2/p0)
        end if
        if(abs(phi_gz(iL)%p2).gt.1.d-8) then
           entropy = entropy - abs(phi_gz(iL)%p2)**2*log(abs(phi_gz(iL)%p2)**2/p2)
        end if
        if(abs(phi_gz(iL)%p1).gt.1.d-8) then
           entropy = entropy - abs(phi_gz(iL)%p1)**2*log(abs(phi_gz(iL)%p1)**2/p1)
        end if
        !+-----------------------------------+!
     end do

     !+- free energy -+!
     f = e_hop + e_hubb - Temp*(entropy + free_entropy)
     !+---------------+!
    
     write(11,*) f

  else

     !+- CONSTRAINTS ON GUTZWILLER PARAMETERS -+!
     if(i.le.L) then
        f = abs(phi_gz(i)%p0)**2 - abs(phi_gz(i)%p2)**2 - tmp_dop(i)        !+----------------------------------------+!
     else
        f = abs(phi_gz(i-L)%p0)**2 + abs(phi_gz(i-L)%p2)**2 + abs(phi_gz(i-L)%p1)**2 - 1.d0
     end if

  end if


contains



  FUNCTION NSTAR(mu)
    USE global
    real(8),dimension(:) :: mu
    real(8),dimension(size(mu)) :: nstar


    real(8) :: fstar_lagrange
    real(8) :: entropy,energy
    real(8) :: n_temp




    hop_intra   = 0.d0
    delta_plus  = 0.d0
    delta_minus = 0.d0
    nstar = 0.d0
    do ik=1,Nk_tot !+- k points loop -+!
       k = vec_k(ik)
       ek = square_lattice_disp(k)
       H_star=0.d0
       w = 0.d0
       do iL=1,L
          H_star(iL,iL) =  ek*abs(Rhop(iL))**2 + mu(iL)
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
       do iL = 1,L
          do jL=1,L
             nstar(jL) = nstar(jL) + abs(H_star(jL,iL))**2*fermi(w(iL),beta)*wt(ik)*2.d0
             hop_intra(jL) = hop_intra(jL) + 2.d0*ek*abs(H_star(jL,iL))**2*wt(ik)*fermi(w(iL),beta)
             if(jL.lt.L) then
                delta_plus(jL)  = delta_plus(jL) + 2.d0*conjg(H_star(jL+1,iL))*H_star(jL,iL)*wt(ik)*fermi(w(iL),beta)
                delta_minus(jL) = delta_minus(jL) + 2.d0*conjg(H_star(jL,iL))*H_star(jL+1,iL)*wt(ik)*fermi(w(iL),beta)
             end if
          end do
       end do
    end do !+- end k points loop -+! 

    nstar = nstar - (1.d0 - tmp_dop)

    return
  END FUNCTION NSTAR




end subroutine free_energy_GZ_parameters







!#################################################################################################!
!#################################################################################################!
!#################################################################################################!
!#################################################################################################!
!#################################################################################################!
!#################################################################################################!
!#################################################################################################!
!#################################################################################################!
!#################################################################################################!
!#################################################################################################!
!#################################################################################################!
!#################################################################################################!
!#################################################################################################!
!#################################################################################################!
!#################################################################################################!
!#################################################################################################!
!#################################################################################################!
!#################################################################################################!
!#################################################################################################!


subroutine fun ( x, f, i )
  !.............................................................................
  use global
  real(8), intent( in )   :: x( : )
  real(8), intent( out )  :: f
  integer, intent( in ), optional   :: i

  if ( .not. present( i ) ) then
     !       the objective function value (user defined)
     !==============================================================================
     f = 100.0d0*(x(2)-x(1)**2)**2 +(1.0d0-x(1))**2                      !
     !==============================================================================
  else
     select case ( i )
     case ( 1 )
        !               the equality constraint value (user defined)
        !==============================================================================
        f = x(1)+3.0d0*x(2)-3.0d0                                   !
        !==============================================================================
     case ( 2 )
        !               the inequality constraint value (user defined)
        !==============================================================================
        f = x(1)**2+x(2)**2-4.0d0                                    !
        !==============================================================================
     end select
  end if
  return
end subroutine fun
