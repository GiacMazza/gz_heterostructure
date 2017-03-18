  !+-----------------------------------------------------------------------------------------+!
  !+- SLAB AT EQUILIBRIUM USING GUTZWILLER ESTIMATION OF THE FINITE TEMPERATURE FREE ENERGY -+!
  !+- Author: G. Mazza @ SISSA                                                              -+!
!+-----------------------------------------------------------------------------------------+!
MODULE SLAB_EQU_TEMP
  USE MIN_AMOEBA
  USE GLOBAL
  USE LANCELOT_simple_double
  implicit none
  private
  integer :: amoeba_funct_calls
  logical :: print_iter
  real(8),allocatable :: dop_ext(:)
  public  :: GZ_equilibrium_temp
  public  :: GZ_equilibrium_temp_china

CONTAINS


  SUBROUTINE GZ_EQUILIBRIUM_TEMP
    integer :: iL
    allocate(local_field(L))
    select case(inner_field_type)
    case('none')
       local_field = 0.d0
       call GZ_EQUILIBRIUM_TEMP_PH
       write(*,*) 'minimization achived at ph symmetry'
    case('linear')
       do iL=1,L
          !local_field(iL) = -(mu_L + (mu_R - mu_L )*dble(iL-1)/dble(L-1))
          local_field(iL) = chem_equ
       end do
       call GZ_EQUILIBRIUM_TEMP_DENSITY
       write(*,*) 'minimization achived at free layer density'
    end select
    deallocate(local_field)
  END SUBROUTINE GZ_EQUILIBRIUM_TEMP



  !+--------------------------------------------------------------+!
  !   FULL MINIMIZATION OF THE GUTZWILLER FREE ENERGY ESTIMATION   !
  !+--------------------------------------------------------------+!
  SUBROUTINE GZ_EQUILIBRIUM_TEMP_DENSITY
    real (8), allocatable, dimension(:) :: x
    real (8), allocatable, dimension(:) :: x_start_ins,x_start_met

    real(8),allocatable,dimension(:,:)  :: p
    real(8),allocatable,dimension(:)    :: y

    real(8)                             :: f_start_ins,f_start_met,f_test
    real(8)                             :: ftol

    integer                             :: iL,jL,i_theta
    integer                             :: np,mp,n_min
    integer                             :: i,iter,j,ndim,i_vertex,i_dim

    real(8)                             :: d_max,d_min


    allocate(mu_star(L),hop_intra(L),delta_plus(L),delta_minus(L),tmp_dop(L),phi_out(L))
    allocate(mu_met(L),mu_ins(L))

    print_iter=.true.

    d_max = -0.15
    d_min = -0.01

    allocate(x(L),x_start_met(L),x_start_ins(L))

    NP = L
    MP = NP + 1
    allocate(p(MP,NP),y(MP))

    amoeba_funct_calls = 0

    !+-----------------------------------------------------+!
    !+- MINIMIZATION STARTING FROM METALLIC-LIKE SOLUTION -+!
    !+-----------------------------------------------------+!
    do iL=1,L/2
       if(Uz(iL).lt.16.d0) then
          x(iL)     = acos(Uz(iL)/16.d0)
          x(iL+L/2) = pi*0.2  !0.5d0*acos(-( d_max+(d_min-d_max)/(L/2-1)*(iL-1) ) / sin(x(iL))**2 )
       else
          x(iL)     = 0.5d0
          x(iL+L/2) = pi*0.25d0
       end if
    end do
    open(11,file='free_iter_met.out')
    open(10,file='profile_iter_met.out')

    !+- initialize simplex -+!
    do i_dim = 1,NP
       p(:,i_dim) = x(i_dim)
    end do

    do i_vertex=2,MP
       p(i_vertex,i_vertex-1) = 0.1d0*pi
    end do
    do i_vertex=1,MP
       y(i_vertex)=f_theta_gamma(p(i_vertex,:))
       write(*,*) 'vertex',i_vertex,'free energy',y(i_vertex)
    end do

    !+- 1ST minimization attempt -+!
    ftol=1.d-6
    call amoeba(p(1:MP,1:NP),y(1:MP),FTOL,f_theta_gamma,iter)

    !+- reinitialize symplex -+!
    do i_dim = 1,NP
       p(:,i_dim) = p(1,i_dim)
    end do
    do i_vertex=2,MP
       p(i_vertex,i_vertex-1) = 0.25d0*pi
    end do
    do i_vertex=1,MP
       y(i_vertex)=f_theta_gamma(p(i_vertex,:))
       write(*,*) 'vertex',i_vertex,'free energy',y(i_vertex)
    end do
    !+- 2ND attempt -+!
    ftol = 1.d-8
    call amoeba(p(1:MP,1:NP),y(1:MP),FTOL,f_theta_gamma,iter)

    !+- reinitialize symplex -+!
    do i_dim = 1,NP
       p(:,i_dim) = p(1,i_dim)
    end do
    do i_vertex=2,MP
       p(i_vertex,i_vertex-1) = 0.25d0*pi
    end do
    do i_vertex=1,MP
       y(i_vertex)=f_theta_gamma(p(i_vertex,:))
       write(*,*) 'vertex',i_vertex,'free energy',y(i_vertex)
    end do
    !+- 3RD attemp -+!
    ftol = 1.d-9
    call amoeba(p(1:MP,1:NP),y(1:MP),FTOL,f_theta_gamma,iter)


    x_start_met = p(1,:)
    f_start_met = f_theta_gamma(x_start_met)
    mu_met = mu_star

    close(11)
    close(10)


    amoeba_funct_calls = 0

    !+--------------------------------------------------+!
    !+- MINIMIZATION STARTING FROM INSULATING SOLUTION -+!
    !+--------------------------------------------------+!
    open(11,file='free_iter_ins.out')
    open(10,file='profile_iter_ins.out')
    do iL=1,L/2
       x(iL)     = 0.0001d0
       x(iL+L/2) = pi*0.25d0
    end do

    !+- initialize simplex -+!
    do i_dim = 1,NP
       p(:,i_dim) = x(i_dim)
    end do

    do i_vertex=2,MP
       p(i_vertex,i_vertex-1) = 0.15d0*pi
    end do
    do i_vertex=1,MP
       y(i_vertex)=f_theta_gamma(p(i_vertex,:))
       write(*,*) 'vertex',i_vertex,'free energy',y(i_vertex)
    end do

    !+- 1ST minimization attempt -+!
    ftol=1.d-6
    call amoeba(p(1:MP,1:NP),y(1:MP),FTOL,f_theta_gamma,iter)

    !+- reinitialize symplex -+!
    do i_dim = 1,NP
       p(:,i_dim) = p(1,i_dim)
    end do
    do i_vertex=2,MP
       p(i_vertex,i_vertex-1) = 0.15d0*pi
    end do
    do i_vertex=1,MP
       y(i_vertex)=f_theta_gamma(p(i_vertex,:))
       write(*,*) 'vertex',i_vertex,'free energy',y(i_vertex)
    end do

    !+- 2ND minimization attempt -+!
    ftol = 1.d-8
    call amoeba(p(1:MP,1:NP),y(1:MP),FTOL,f_theta_gamma,iter)

    !+- reinitialize symplex -+!
    do i_dim = 1,NP
       p(:,i_dim) = p(1,i_dim)
    end do
    do i_vertex=2,MP
       p(i_vertex,i_vertex-1) = 0.15d0*pi
    end do
    do i_vertex=1,MP
       y(i_vertex)=f_theta_gamma(p(i_vertex,:))
       write(*,*) 'vertex',i_vertex,'free energy',y(i_vertex)
    end do

    !+- 3RD minimization attempt -+!
    ftol = 1.d-12
    call amoeba(p(1:MP,1:NP),y(1:MP),FTOL,f_theta_gamma,iter)

    x_start_ins = p(1,:)
    f_start_ins = f_theta_gamma(x_start_ins)
    mu_ins = mu_star

    if(f_start_ins.ge.f_start_met) then
       write(*,*) 'METALLIC-LIKE SOLUTION',f_start_ins-f_start_met
       mu_star = mu_met
       call get_equ_temp_conditions(x_start_met)
    else
       mu_star = mu_ins
       call get_equ_temp_conditions(x_start_ins)
       write(*,*) 'INSULATING-LIKE SOLUTION',f_start_met-f_start_ins
    end if
    close(11)
    close(10)
    deallocate(y,p)
    deallocate(x,x_start_ins,x_start_met)

    deallocate(mu_star,hop_intra,delta_plus,delta_minus,tmp_dop,phi_out)
    deallocate(mu_met,mu_ins)

  END SUBROUTINE GZ_EQUILIBRIUM_TEMP_DENSITY



  SUBROUTINE GZ_EQUILIBRIUM_TEMP_PH
    real (8), allocatable, dimension(:) :: x
    real (8), allocatable, dimension(:) :: x_start_ins,x_start_met

    real(8),allocatable,dimension(:,:)  :: p
    real(8),allocatable,dimension(:)    :: y

    real(8)                             :: f_start_ins,f_start_met,f_test
    real(8)                             :: ftol

    integer                             :: iL,jL,i_theta
    integer                             :: np,mp,n_min
    integer                             :: i,iter,j,ndim,i_vertex,i_dim

    real(8)                             :: d_max,d_min


    allocate(mu_star(L),hop_intra(L),delta_plus(L),delta_minus(L),tmp_dop(L),phi_out(L))
    allocate(mu_met(L),mu_ins(L))

    print_iter=.true.


    allocate(x(L),x_start_met(L),x_start_ins(L))

    NP = L/2
    MP = NP + 1
    allocate(p(MP,NP),y(MP))

    amoeba_funct_calls = 0

    !+-----------------------------------------------------+!
    !+- MINIMIZATION STARTING FROM METALLIC-LIKE SOLUTION -+!
    !+-----------------------------------------------------+!
    do iL=1,L/2
       if(Uz(iL).lt.16.d0) then
          x(iL)     = acos(Uz(iL)/16.d0)
       else
          x(iL)     = 0.00001d0
       end if
    end do
    open(11,file='free_iter_met_ph.out')
    open(10,file='profile_iter_met_ph.out')

    !+- initialize simplex -+!
    do i_dim = 1,NP
       p(:,i_dim) = x(i_dim)
    end do

    do i_vertex=2,MP
       p(i_vertex,i_vertex-1) = 0.25d0*pi
    end do
    do i_vertex=1,MP
       y(i_vertex)=f_theta(p(i_vertex,:))
       write(*,*) 'vertex',i_vertex,'free energy',y(i_vertex)
    end do

    !+- 1ST minimization attempt -+!
    ftol=1.d-6
    call amoeba(p(1:MP,1:NP),y(1:MP),FTOL,f_theta,iter)

    !+- reinitialize symplex -+!
    do i_dim = 1,NP
       p(:,i_dim) = p(1,i_dim)
    end do
    do i_vertex=2,MP
       p(i_vertex,i_vertex-1) = 0.25d0*pi
    end do
    do i_vertex=1,MP
       y(i_vertex)=f_theta(p(i_vertex,:))
       write(*,*) 'vertex',i_vertex,'free energy',y(i_vertex)
    end do
    !+- 2ND attempt -+!
    ftol = 1.d-8
    call amoeba(p(1:MP,1:NP),y(1:MP),FTOL,f_theta,iter)

    !+- reinitialize symplex -+!
    do i_dim = 1,NP
       p(:,i_dim) = p(1,i_dim)
    end do
    do i_vertex=2,MP
       p(i_vertex,i_vertex-1) = 0.25d0*pi
    end do
    do i_vertex=1,MP
       y(i_vertex)=f_theta(p(i_vertex,:))
       write(*,*) 'vertex',i_vertex,'free energy',y(i_vertex)
    end do
    !+- 3RD attemp -+!
    ftol = 1.d-12
    call amoeba(p(1:MP,1:NP),y(1:MP),FTOL,f_theta,iter)


    x_start_met = p(1,:)
    f_start_met = f_theta(x_start_met)
    mu_met = mu_star

    close(11)
    close(10)


    amoeba_funct_calls = 0

    !+--------------------------------------------------+!
    !+- MINIMIZATION STARTING FROM INSULATING SOLUTION -+!
    !+--------------------------------------------------+!
    open(11,file='free_iter_ins_ph.out')
    open(10,file='profile_iter_ins_ph.out')
    do iL=1,L/2
       x(iL)     = 0.0001d0
    end do

    !+- initialize simplex -+!
    do i_dim = 1,NP
       p(:,i_dim) = x(i_dim)
    end do
    do i_vertex=2,MP
       p(i_vertex,i_vertex-1) = 0.15d0*pi
    end do
    do i_vertex=1,MP
       y(i_vertex)=f_theta(p(i_vertex,:))
       write(*,*) 'vertex',i_vertex,'free energy',y(i_vertex)
    end do
    !+- 1ST minimization attempt -+!
    ftol=1.d-6
    call amoeba(p(1:MP,1:NP),y(1:MP),FTOL,f_theta,iter)
    !+- reinitialize symplex -+!
    do i_dim = 1,NP
       p(:,i_dim) = p(1,i_dim)
    end do
    do i_vertex=2,MP
       p(i_vertex,i_vertex-1) = 0.15d0*pi
    end do
    do i_vertex=1,MP
       y(i_vertex)=f_theta(p(i_vertex,:))
       write(*,*) 'vertex',i_vertex,'free energy',y(i_vertex)
    end do
    !+- 2ND minimization attempt -+!
    ftol = 1.d-8
    call amoeba(p(1:MP,1:NP),y(1:MP),FTOL,f_theta,iter)
    !+- reinitialize symplex -+!
    do i_dim = 1,NP
       p(:,i_dim) = p(1,i_dim)
    end do
    do i_vertex=2,MP
       p(i_vertex,i_vertex-1) = 0.15d0*pi
    end do
    do i_vertex=1,MP
       y(i_vertex)=f_theta(p(i_vertex,:))
       write(*,*) 'vertex',i_vertex,'free energy',y(i_vertex)
    end do
    !+- 3RD minimization attempt -+!
    ftol = 1.d-12
    call amoeba(p(1:MP,1:NP),y(1:MP),FTOL,f_theta,iter)
    x_start_ins = p(1,:)
    f_start_ins = f_theta(x_start_ins)
    mu_ins = mu_star


    !+- CHECK SOLUTIONS -+!
    if(f_start_ins.ge.f_start_met) then
       write(*,*) 'METALLIC-LIKE SOLUTION'
       write(*,*) 'F_ins - F_met = ',f_start_ins-f_start_met
       mu_star = mu_met
       call get_equ_temp_conditions_ph(x_start_met)
    else
       mu_star = mu_ins
       call get_equ_temp_conditions_ph(x_start_ins)
       write(*,*) 'INSULATING-LIKE SOLUTION',f_start_met-f_start_ins
       write(*,*) 'F_met - F_ins = ',f_start_met-f_start_ins
    end if
    !+-------------------+!

    close(11)
    close(10)
    deallocate(y,p)
    deallocate(x,x_start_ins,x_start_met)
    deallocate(mu_star,hop_intra,delta_plus,delta_minus,tmp_dop,phi_out)
    deallocate(mu_met,mu_ins)
  END SUBROUTINE GZ_EQUILIBRIUM_TEMP_PH





  SUBROUTINE GZ_EQUILIBRIUM_TEMP_CHINA
    real (8), allocatable, dimension(:) :: x
    real (8), allocatable, dimension(:) :: x_start_ins,x_start_met

    real(8),allocatable,dimension(:,:)  :: p
    real(8),allocatable,dimension(:)    :: y

    real(8)                             :: f_start_ins,f_start_met,f_test
    real(8)                             :: ftol

    integer                             :: iL,jL,i_theta
    integer                             :: np,mp,n_min
    integer                             :: i,iter,j,ndim,i_vertex,i_dim

    real(8)                             :: d_max,d_min


    allocate(mu_star(L),hop_intra(L),delta_plus(L),delta_minus(L),tmp_dop(L),phi_out(L))
    allocate(mu_met(L),mu_ins(L))

    print_iter=.true.


    allocate(x(L),x_start_met(L),x_start_ins(L))

    NP = L/2
    MP = NP + 1
    allocate(p(MP,NP),y(MP))

    amoeba_funct_calls = 0

    !+-----------------------------------------------------+!
    !+- MINIMIZATION STARTING FROM METALLIC-LIKE SOLUTION -+!
    !+-----------------------------------------------------+!
    do iL=1,L/2
       if(Uz(iL).lt.16.d0) then
          x(iL)     = acos(Uz(iL)/16.d0)
       else
          x(iL)     = 0.00011d0
       end if
    end do
    open(11,file='free_iter_met_ph.out')
    open(10,file='profile_iter_met_ph.out')

    !+- initialize simplex -+!
    do i_dim = 1,NP
       p(:,i_dim) = x(i_dim)
    end do

    do i_vertex=2,MP
       p(i_vertex,i_vertex-1) = 0.25d0*pi
    end do
    do i_vertex=1,MP
       y(i_vertex)=f_theta(p(i_vertex,:))
       write(*,*) 'vertex',i_vertex,'free energy',y(i_vertex)
    end do

    !+- 1ST minimization attempt -+!
    ftol=1.d-6
    call amoeba(p(1:MP,1:NP),y(1:MP),FTOL,f_theta,iter)

    !+- reinitialize symplex -+!
    do i_dim = 1,NP
       p(:,i_dim) = p(1,i_dim)
    end do
    do i_vertex=2,MP
       p(i_vertex,i_vertex-1) = 0.25d0*pi
    end do
    do i_vertex=1,MP
       y(i_vertex)=f_theta(p(i_vertex,:))
       write(*,*) 'vertex',i_vertex,'free energy',y(i_vertex)
    end do
    !+- 2ND attempt -+!
    ftol = 1.d-8
    call amoeba(p(1:MP,1:NP),y(1:MP),FTOL,f_theta,iter)

    !+- reinitialize symplex -+!
    do i_dim = 1,NP
       p(:,i_dim) = p(1,i_dim)
    end do
    do i_vertex=2,MP
       p(i_vertex,i_vertex-1) = 0.25d0*pi
    end do
    do i_vertex=1,MP
       y(i_vertex)=f_theta(p(i_vertex,:))
       write(*,*) 'vertex',i_vertex,'free energy',y(i_vertex)
    end do
    !+- 3RD attemp -+!
    ftol = conv_treshold
    call amoeba(p(1:MP,1:NP),y(1:MP),FTOL,f_theta,iter)


    x_start_met = p(1,:)
    f_start_met = f_theta(x_start_met)
    mu_met = mu_star

    close(11)
    close(10)


    
    mu_star = mu_met
    call get_equ_temp_conditions_ph(x_start_met)

    deallocate(y,p)
    deallocate(x,x_start_ins,x_start_met)
    deallocate(mu_star,hop_intra,delta_plus,delta_minus,tmp_dop,phi_out)
    deallocate(mu_met,mu_ins)
  END SUBROUTINE GZ_EQUILIBRIUM_TEMP_CHINA





  FUNCTION F_THETA_GAMMA(gz_angle)
    real(8),intent(in)                          :: gz_angle(:)
    real(8)                                     :: f_theta_gamma

    !+- inner variables -+!
    type(gz_projector),dimension(:),allocatable :: phi_gz
    real(8),dimension(:),allocatable            :: opt_gz,gz_dop,Docc,n_test
    complex(8),dimension(:),allocatable         :: Rhop
    complex(8),allocatable,dimension(:,:)       :: H_star
    real(8),allocatable,dimension(:)            :: w

    !+- correction term -+!
    !real(8),allocatable,dimension(:)            :: e0_temp
    real(8)                                     :: e0_temp


    real(8)                                     :: theta,gm,phi
    real(8)                                     :: ek
    real(8)                                     :: free_entropy,prj_entropy,energy
    real(8)                                     :: p0,p1,p2
    real(8)                                     :: free
    real(8)                                     :: n_temp

    type(vec2D)                                 :: k

    integer                                     :: iL,jL
    integer                                     :: ik

    allocate(phi_gz(L),H_star(L,L),w(L),Rhop(L),gz_dop(L),Docc(L),n_test(L))

    do iL=1,L/2
       theta = gz_angle(iL)
       gm = gz_angle(iL+L/2)
       phi=0.d0
       ! left 
       phi_gz(iL)%p0 = sin(theta*0.5d0)*sin(gm)*Exp(Zi*phi)
       phi_gz(iL)%p2 = sin(theta*0.5d0)*cos(gm)*Exp(Zi*phi)
       phi_gz(iL)%p1 = cos(theta*0.5d0)
       ! right 
       ! phi_gz(L-(iL-1))%p0 = sin(theta*0.5d0)*sin(pi*0.5d0-gm)*Exp(Zi*phi)
       ! phi_gz(L-(iL-1))%p2 = sin(theta*0.5d0)*cos(pi*0.5d0-gm)*Exp(Zi*phi)
       ! phi_gz(L-(iL-1))%p1 = cos(theta*0.5d0)
       phi_gz(L-(iL-1))%p0 = sin(theta*0.5d0)*sin(gm)*Exp(Zi*phi)
       phi_gz(L-(iL-1))%p2 = sin(theta*0.5d0)*cos(gm)*Exp(Zi*phi)
       phi_gz(L-(iL-1))%p1 = cos(theta*0.5d0)
    end do

    Rhop   = GZ_hop(phi_gz)
    gz_dop = GZ_doping(phi_gz)
    Docc   = GZ_double(phi_gz)


    
    !+- minimize F_star with respect to layer dependent chemical potentials -+!
    mu_star = 0.d0  
    call fsolve(nstar_slab,mu_star,tol=1.d-15)
    !+-----------------------------------------------------------------------+!


    !+------------------------+!
    !+- FREE ENERGY ESTIMATE -+!
    !+------------------------+!

    free = 0.d0

    !+- trace over uncorrelated distribution -+!
    n_test       = 0.d0
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
             n_test(jL) = n_test(jL) + 2.d0*abs(H_star(jL,iL))**2*wt(ik)*fermi(w(iL),beta)
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
    
    if(print_iter) then
       do iL=1,L
          write(10,'(10(F18.10))') dble(iL),mu_star(iL),gz_dop(iL),1.d0-n_test(iL),abs(Rhop(iL))
       end do
       write(10,*)
       write(10,*)
    end if
    
    !+- GZ projectors contribution to the free energy -+!
    prj_entropy = 0.d0
    do iL=1,L

       !+- energy contribution to free energy -+!
       free = free + hop_intra(iL)*ABS(Rhop(iL))**2 
       if(iL.gt.1) then
          free = free - conjg(Rhop(iL-1))*Rhop(iL)*delta_minus(iL-1)
       end if
       if(iL.lt.L) then
          free = free - conjg(Rhop(iL+1))*Rhop(iL)*delta_plus(iL)
       end if
       free = free + Uz(iL)*Docc(iL) + local_field(iL)*(1.d0-gz_dop(iL))
       !+--------------------------------------+!

       !+- projector entropy -+!
       p0 = (1.d0 + gz_dop(iL))**2.d0/4.d0
       p1 = (1.d0 - gz_dop(iL)**2.d0)/2.d0
       p2 = (1.d0 - gz_dop(iL))**2.d0/4.d0

       if(abs(phi_gz(iL)%p0).gt.1.d-10) then
          prj_entropy = prj_entropy - abs(phi_gz(iL)%p0)**2.d0*log(abs(phi_gz(iL)%p0)**2.d0/p0) 
       end if

       if(abs(phi_gz(iL)%p2).gt.1.d-10) then
          prj_entropy = prj_entropy - abs(phi_gz(iL)%p2)**2.d0*log(abs(phi_gz(iL)%p2)**2.d0/p2) 
       end if

       if(abs(phi_gz(iL)%p1).gt.1.d-10) then
          prj_entropy = prj_entropy - abs(phi_gz(iL)%p1)**2.d0*log(abs(phi_gz(iL)%p1)**2.d0/p1) 
       end if
       !+---------------------+!
    end do

    !+- internal energy -+!
    energy = free

    !+- free energy estimation -+!
    !f_theta_gamma = energy - temp*(prj_entropy + free_entropy)
    f_theta_gamma = energy
    
    amoeba_funct_calls = amoeba_funct_calls + 1
    write(11,'(6f18.10)') dble(amoeba_funct_calls),f_theta_gamma,energy,free_entropy,prj_entropy
    write(*,*) amoeba_funct_calls,f_theta_gamma

    deallocate(phi_gz,H_star,w,Rhop,gz_dop,Docc)

  CONTAINS

    FUNCTION NSTAR_SLAB(mu)
      real(8),dimension(:)                  :: mu
      real(8),dimension(size(mu))           :: nstar_slab

      real(8)                               :: fstar_lagrange
      real(8)                               :: entropy,energy
      real(8)                               :: n_temp


      complex(8),dimension(L,L) :: H_star
      real(8),dimension(L)      :: w

      nstar_slab = 0.d0
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
               nstar_slab(jL) = nstar_slab(jL) + abs(H_star(jL,iL))**2*fermi(w(iL),beta)*wt(ik)*2.d0
            end do
         end do
      end do !+- end k points loop -+! 

      nstar_slab = nstar_slab - (1.d0 - gz_dop)
      return
    END FUNCTION NSTAR_SLAB
  END FUNCTION F_THETA_GAMMA




  FUNCTION F_THETA(gz_angle)
    real(8),intent(in)                          :: gz_angle(:)
    real(8)                                     :: f_theta
    !+- inner variables -+!
    type(gz_projector),dimension(:),allocatable :: phi_gz
    real(8),dimension(:),allocatable            :: opt_gz,gz_dop,Docc,n_test
    complex(8),dimension(:),allocatable         :: Rhop
    complex(8),allocatable,dimension(:,:)       :: H_star
    real(8),allocatable,dimension(:)            :: w
    real(8)                                     :: theta,gm,phi
    real(8)                                     :: ek
    real(8)                                     :: free_entropy,prj_entropy,energy
    real(8)                                     :: p0,p1,p2
    real(8)                                     :: free
    real(8)                                     :: n_temp
    type(vec2D)                                 :: k
    integer                                     :: iL,jL
    integer                                     :: ik
    real(8)                                     :: e0_temp
    real(8)                                     :: ehop_temp,D_ave
    real(8),dimension(:),allocatable            :: e0_intra,e0_plus,e0_minus
    allocate(phi_gz(L),H_star(L,L),w(L),Rhop(L),gz_dop(L),Docc(L),n_test(L))
    allocate(e0_intra(L),e0_plus(L),e0_minus(L))

    e0_intra=0.d0
    e0_plus=0.d0
    e0_minus=0.d0
    do ik=1,Nk_tot !+- k points loop -+!
       !+- build and diagonalize real space hamiltonian -+!
       k = vec_k(ik)
       ek = square_lattice_disp(k)
       H_star=0.d0
       w = 0.d0
       do iL=1,L
          H_star(iL,iL)   = ek
          if(iL.lt.L) then
             H_star(iL,iL+1) = -1.d0
             H_star(iL+1,iL) = -1.d0
          end if
       end do
       call  matrix_diagonalize(H_star,w,'V','L')
       !+------------------------------------------------+!
       do iL=1,L
          n_temp = fermi(w(iL),beta)
          do jL=1,L
             n_test(jL) = n_test(jL) + 2.d0*abs(H_star(jL,iL))**2*wt(ik)*n_temp
             e0_intra(jL) = e0_intra(jL) + 2.d0*ek*abs(H_star(jL,iL))**2*wt(ik)*n_temp
             if(jL.lt.L) then
                e0_plus(jL)  = delta_plus(jL) + 2.d0*conjg(H_star(jL+1,iL))*H_star(jL,iL)*wt(ik)*n_temp
                e0_minus(jL) = delta_minus(jL) + 2.d0*conjg(H_star(jL,iL))*H_star(jL+1,iL)*wt(ik)*n_temp
             endif
          enddo
       enddo
    enddo !+- end k points loop -+! 

    !+--+!
    ! n_temp = fermi(ene_dos(ieps)*(R**2 + Gamma*(1.d0 - 4.d0*double)**2.d0*n0),beta)
    ! energy_free = energy_free + (R**2+ Gamma*(1.d0 - 4.d0*double)**2.d0*n0)*n_temp*ene_dos(ieps)*flat_dos(ieps)*2.d0
    !+--+!


    do iL=1,L/2
       theta = gz_angle(iL)
       gm = 0.25d0*pi
       phi=0.d0
       phi_gz(iL)%p0 = sin(theta*0.5d0)*sin(gm)*Exp(Zi*phi)
       phi_gz(iL)%p2 = sin(theta*0.5d0)*cos(gm)*Exp(Zi*phi)
       phi_gz(iL)%p1 = cos(theta*0.5d0)
       !
       phi_gz(L-(iL-1))%p0 = sin(theta*0.5d0)*sin(gm)*Exp(Zi*phi)
       phi_gz(L-(iL-1))%p2 = sin(theta*0.5d0)*cos(gm)*Exp(Zi*phi)
       phi_gz(L-(iL-1))%p1 = cos(theta*0.5d0)
       !
    end do
    Rhop   = GZ_hop(phi_gz)
    gz_dop = GZ_doping(phi_gz)
    Docc   = GZ_double(phi_gz)
    D_ave = sum(Docc)/dble(L)
    !+- minimize F_star with respect to layer dependent chemical potentials -+!
    mu_star = 0.d0  
    !+------------------------+!
    !+- FREE ENERGY ESTIMATE -+!
    !+------------------------+!
    free = 0.d0
    !+- trace over uncorrelated distribution -+!
    n_test       = 0.d0
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
          H_star(iL,iL) =  ek*abs(Rhop(iL))**2 + mu_star(iL) !- Gamma*(1.d0-4.d0*Docc(iL))**2.d0*e0_intra(iL)
          if(iL.lt.L) then
             H_star(iL,iL+1) = -1.d0*conjg(Rhop(iL))*Rhop(iL+1) !- Gamma*(1.d0-4.d0*Docc(iL))*(1.d0-4.d0*Docc(iL+1))*e0_plus(iL)
             H_star(iL+1,iL) = -1.d0*Rhop(iL)*conjg(Rhop(iL+1)) !- Gamma*(1.d0-4.d0*Docc(iL))*(1.d0-4.d0*Docc(iL+1))*e0_minus(iL)
          end if
       end do
       if(pbc) then
          H_star(1,L) = -1.d0*conjg(Rhop(1))*Rhop(L)
          H_star(L,1) = -1.d0*conjg(Rhop(L))*Rhop(1)
       end if
       call  matrix_diagonalize(H_star,w,'V','L')
       !+------------------------------------------------+!
       do iL=1,L
          n_temp = fermi(w(iL),beta)
          !n_temp = fermi(w(iL)- Gamma*(1.d0-4.d0*Docc(iL))**2.d0*e0_intra(iL),beta)
          do jL=1,L
             n_test(jL) = n_test(jL) + 2.d0*abs(H_star(jL,iL))**2*wt(ik)*n_temp
             hop_intra(jL) = hop_intra(jL) + 2.d0*ek*abs(H_star(jL,iL))**2*wt(ik)*n_temp
             if(jL.lt.L) then
                delta_plus(jL)  = delta_plus(jL) + 2.d0*conjg(H_star(jL+1,iL))*H_star(jL,iL)*wt(ik)*n_temp
                delta_minus(jL) = delta_minus(jL) + 2.d0*conjg(H_star(jL,iL))*H_star(jL+1,iL)*wt(ik)*n_temp
                !+- free electrons entropy -+!
                if(jL.eq.1) then
                   if(n_temp.gt.1.d-10) then
                      free_entropy = free_entropy - ( n_temp*log(n_temp))*wt(ik)*2.d0 
                   endif
                   if(abs(1.d0-n_temp).gt.1.d-10) then
                      free_entropy = free_entropy - (1.d0-n_temp)*log(1.d0-n_temp )*wt(ik)*2.d0
                   endif
                endif
                !+--------------------------+!
             endif
          enddo
       enddo
    enddo !+- end k points loop -+! 
    if(print_iter) then
       do iL=1,L
          write(10,'(10(F18.10))') dble(iL),mu_star(iL),gz_dop(iL),1.d0-n_test(iL),abs(Rhop(iL))
       end do
       write(10,*)
       write(10,*)
    end if
    !+- GZ projectors contribution to the free energy -+!
    prj_entropy = 0.d0
    do iL=1,L
       !+- energy contribution to free energy -+!
       write(*,*) Gamma,D_ave
       free = free + hop_intra(iL)*(ABS(Rhop(iL))**2 - Gamma*(1.d0-4.d0*Docc(iL))**2.d0*e0_intra(iL))
       if(iL.gt.1) then
          free = free - (conjg(Rhop(iL-1))*Rhop(iL) - Gamma*(1.d0-4.d0*Docc(iL))*(1.d0-4.d0*Docc(iL-1))*e0_intra(iL))*delta_minus(iL-1)
       end if
       if(iL.lt.L) then
          free = free - (conjg(Rhop(iL+1))*Rhop(iL) - Gamma*(1.d0-4.d0*Docc(iL))*(1.d0-4.d0*Docc(iL+1))*e0_intra(iL))*delta_plus(iL)
       end if
       free = free + Uz(iL)*Docc(iL) + local_field(iL)*(1.d0-gz_dop(iL))
       !+--------------------------------------+!
       !+- projector entropy -+!
       p0 = (1.d0 + gz_dop(iL))**2.d0/4.d0
       p1 = (1.d0 - gz_dop(iL)**2.d0)/2.d0
       p2 = (1.d0 - gz_dop(iL))**2.d0/4.d0
       if(abs(phi_gz(iL)%p0).gt.1.d-10) then
          prj_entropy = prj_entropy - abs(phi_gz(iL)%p0)**2.d0*log(abs(phi_gz(iL)%p0)**2.d0/p0) 
       end if
       if(abs(phi_gz(iL)%p2).gt.1.d-10) then
          prj_entropy = prj_entropy - abs(phi_gz(iL)%p2)**2.d0*log(abs(phi_gz(iL)%p2)**2.d0/p2) 
       end if
       if(abs(phi_gz(iL)%p1).gt.1.d-10) then
          prj_entropy = prj_entropy - abs(phi_gz(iL)%p1)**2.d0*log(abs(phi_gz(iL)%p1)**2.d0/p1) 
       end if
       !+---------------------+!
    end do
    !+- internal energy -+!
    energy = free
    !+- free energy estimation -+!
    f_theta = energy - temp*(prj_entropy + free_entropy)
    amoeba_funct_calls = amoeba_funct_calls + 1
    write(11,'(6f18.10)') dble(amoeba_funct_calls),f_theta,energy,free_entropy,prj_entropy
    write(*,*) amoeba_funct_calls,f_theta
    deallocate(phi_gz,H_star,w,Rhop,gz_dop,Docc)
  END FUNCTION F_THETA









  SUBROUTINE get_equ_temp_conditions(x_opt)
    real(8),dimension(:)      :: x_opt
    complex(8),dimension(L,L) :: H_star

    complex(8),dimension(L)   :: R_opt
    real(8),dimension(L)      :: w,dop_opt

    integer                   :: iL,jL,iE,ik
    real(8)                   :: theta,gamma,phi,ek
    type(vec2D)               :: k
    real(8)                   :: Rave


    !+- get GZ parameters -+!
    do iL=1,L/2
       theta = x_opt(iL)
       gamma = x_opt(iL+L/2)
       phi=0.d0
       ! left 
       eqPhi(iL)%p0 = sin(theta*0.5d0)*sin(gamma)*Exp(Zi*phi)
       eqPhi(iL)%p2 = sin(theta*0.5d0)*cos(gamma)*Exp(Zi*phi)
       eqPhi(iL)%p1 = cos(theta*0.5d0)
       ! right 
       eqPhi(L-(iL-1))%p0 = sin(theta*0.5d0)*sin(pi*0.5d0-gamma)*Exp(Zi*phi)
       eqPhi(L-(iL-1))%p2 = sin(theta*0.5d0)*cos(pi*0.5d0-gamma)*Exp(Zi*phi)
       eqPhi(L-(iL-1))%p1 = cos(theta*0.5d0)
    end do
    R_opt = GZ_hop(eqPhi)
    dop_opt = GZ_doping(eqPhi)

    open(20,file='out_equ_temp.out')
    do iL=1,L
       write(20,*) iL,abs(R_opt(iL))**2.d0,dop_opt(iL)
    end do
    close(20)

    Rave = sum(R_opt)/dble(L)

    if(Rave.gt.1.d-3) then

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

    else

       do ik=1,Nk_tot !+- k points loop -+!
          k = vec_k(ik)
          ek = square_lattice_disp(k)
          H_star=0.d0
          w = 0.d0
          do iL=1,L
             H_star(iL,iL) =  ek + mu_star(iL)
             if(iL.lt.L) then
                H_star(iL,iL+1) = -1.d0
                H_star(iL+1,iL) = -1.d0
             end if
          end do
          if(pbc) then
             H_star(1,L) = -1.d0
             H_star(L,1) = -1.d0
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

    end if

  END SUBROUTINE get_equ_temp_conditions




  SUBROUTINE get_equ_temp_conditions_ph(x_opt)
    real(8),dimension(:)      :: x_opt
    complex(8),dimension(L,L) :: H_star

    complex(8),dimension(L)   :: R_opt
    real(8),dimension(L)      :: w,dop_opt

    integer                   :: iL,jL,iE,ik
    real(8)                   :: theta,gamma,phi,ek,Rave
    type(vec2D)               :: k



    !+- get GZ parameters -+!
    do iL=1,L/2
       theta = x_opt(iL)
       gamma = 0.25d0*pi
       phi=0.d0
       !
       eqPhi(iL)%p0 = sin(theta*0.5d0)*sin(gamma)*Exp(Zi*phi)
       eqPhi(iL)%p2 = sin(theta*0.5d0)*cos(gamma)*Exp(Zi*phi)
       eqPhi(iL)%p1 = cos(theta*0.5d0)
       !
       eqPhi(L-(iL-1))%p0 = sin(theta*0.5d0)*sin(gamma)*Exp(Zi*phi)
       eqPhi(L-(iL-1))%p2 = sin(theta*0.5d0)*cos(gamma)*Exp(Zi*phi)
       eqPhi(L-(iL-1))%p1 = cos(theta*0.5d0)
       !
    end do
    R_opt = GZ_hop(eqPhi)
    dop_opt = GZ_doping(eqPhi)

    Rave = sum(R_opt)/dble(L)

    open(20,file='out_equ_temp.out')
    do iL=1,L
       write(20,*) iL,abs(R_opt(iL))**2.d0,dop_opt(iL)
    end do
    close(20)

    if(abs(Rave).gt.1.d-4) then

       write(*,*) 'RAVE gtr',Rave


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

    else

       write(*,*) 'RAVE less',Rave


       !+- get uncorrelated density matrix -+!
       do ik=1,Nk_tot !+- k points loop -+!
          k = vec_k(ik)
          ek = square_lattice_disp(k)
          H_star=0.d0
          w = 0.d0
          do iL=1,L
             H_star(iL,iL) =  ek + mu_star(iL)
             if(iL.lt.L) then
                H_star(iL,iL+1) = -1.d0
                H_star(iL+1,iL) = -1.d0
             end if
          end do
          if(pbc) then
             H_star(1,L) = -1.d0
             H_star(L,1) = -1.d0
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

    end if

  END SUBROUTINE get_equ_temp_conditions_ph


  END MODULE SLAB_EQU_TEMP


