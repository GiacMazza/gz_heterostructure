!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BRILLOUIN ZONE STUFF
!Author G Mazza
!Usefull help A Amaricci
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE BZ_POINTS
  USE SF_CONSTANTS
  USE VECTORS
  USE SF_SPECIAL
  USE SF_IOTOOLS
  USE SF_ARRAYS
!  USE FUNCTIONS
  implicit none
  private

  real(8), parameter,public                     :: alat=1.0

  real(8),dimension(:),allocatable,public       :: k_orth   !grid of k_orth points
  type(vec2D),dimension(:,:),allocatable,public :: kgrid    !grid of k vectors
  type(vec2D),dimension(:),allocatable,public   :: vec_k    !one d aray of k vectors
  type(vec2D),dimension(:,:),allocatable,public :: mp_kgrid !mp grid of k vectors
  integer,dimension(:),allocatable,public       :: x_index,y_index
  real(8),dimension(:),allocatable,public       :: wt       !weight of the grid
  real(8),dimension(:),allocatable,public       :: wt_orth  !weight of orth k-points

  real(8),dimension(:),allocatable,public       :: wt_bath
  real(8),dimension(:),allocatable,public       :: wt_bath_log
  real(8),dimension(:),allocatable,public       :: ene_bath
  real(8),dimension(:),allocatable,public       :: ene_bath_log
  real(8),dimension(:),allocatable,public       :: dos_bath
  real(8),dimension(:),allocatable,public       :: dos_bath_log
  real(8),dimension(:),allocatable,public       :: k_bath_log
  real(8),dimension(:),allocatable,public       :: k_bath

  real(8),dimension(:),allocatable,public       :: ene_dos
  real(8),dimension(:),allocatable,public       :: wt_dos


  integer,dimension(:,:),allocatable,public     :: kindex   !k grid index

  public :: reduced_BZ  
  public :: square_lattice_disp
  public :: chain_disp
  public :: fermi_zero
  public :: bath_dos
  public :: bath_coupling
  public :: build_layer_dos
  public :: build_orth_BZ
  !    public :: Npoints_BZ

CONTAINS


  !+------------------------------------------+!
  ! routine to build up the k_grid in the 
  ! Reduced Brillouin Zone
  !+------------------------------------------+!
  SUBROUTINE build_layer_dos(nx,wband,layer_dos_)
    integer               :: nx,i
    real(8)               :: wband
    real(8),dimension(nx) :: dos
    real(8)               :: dw,x,kint,eint,test
    character(len=16),optional     :: layer_dos_
    character(len=16)     :: layer_dos

    layer_dos='bethe'
    if(present(layer_dos_)) layer_dos=layer_dos_

    allocate(wt_dos(nx),ene_dos(nx))

    dos=0.d0
    test=0.d0
    !wband=4.d0
    ene_dos = linspace(-wband*0.5,wband*0.5,nx,mesh=dw)
    do i=1,nx
       x=0.5d0*(ene_dos(i)/wband)**2-1.d0
       call comelp(x,kint,eint)
       !+- select layer density of states -+!
       select case(layer_dos)
       case('2d_square')
          dos(i)=2.d0/wband/pi**2*kint*heaviside(wband-abs(ene_dos(i)))
       case('flat')
          dos(i)=0.5d0/4.d0*heaviside(4.d0-abs(ene_dos(i)))
       case('bethe')
          if(ene_dos(i).ge.-wband/2.0.and.ene_dos(i).le.wband/2.0) then
             dos(i) = 4.d0/Wband/pi*sqrt(1.d0-(2.d0*ene_dos(i)/Wband)**2.d0)
          else
             dos(i) = 0.d0
          end if
          !insert here the bethe lattice dos
       end select
       !
       if(i.eq.nx.or.i.eq.1) then
          wt_dos(i) = dw*0.5d0*dos(i)
       else
          wt_dos(i) = dw*dos(i)
       end if
       !
       test=test+wt_dos(i)
    enddo
    dos=dos/sum(dos)/dw
    call splot("dos2d.out",ene_dos,dos,wt_dos)
!    write(*,*) test;stop
  END SUBROUTINE build_layer_dos



  SUBROUTINE build_orth_BZ(Nk)
    integer   ::  Nk
    real(8)   ::  kz
    real(8)   ::  dk
    real(8)   ::  test
    integer   ::  ik

    allocate(k_orth(Nk),wt_orth(Nk))

    open(unit=10,file='k_orth.out')

    !+- kz = pi/(N+1) ... pi*N/(N+1) -+!
    test = 0.d0

    dk = pi/dble(Nk+1)
    kz = 0.d0
    do ik = 1,Nk
       kz = kz + dk
       k_orth(ik) = kz 

       wt_orth(ik) = 1.d0/dble(Nk)

       write(10,*) k_orth(ik),wt_orth(ik),-2.d0*cos(k_orth(ik))
       test = test + wt_orth(ik)

    end do

    write(10,*) test

    close(10)

  END SUBROUTINE build_orth_BZ



  SUBROUTINE reduced_BZ(nx,Nkorth,off_set)
    integer     :: nx,i
    integer     :: Nk
    integer     :: Nkorth
    type(vec2D) :: ai,aj !real lattice vectors
    type(vec2D) :: bi,bj !reciprocal space vectors
    logical, intent(in) :: off_set
    real(8) :: norm

    ai=Xver*alat
    aj=Yver*alat

    bi=Xver*(2*pi/alat)
    bj=Yver*(2*pi/alat)

    Nk=Npoints_BZ(nx)
    write(*,*) 'number of k points',Nk
    
    call build_reduced_BZ(off_set)
    !call build_orth_BZ(Nkorth)
    call print_reduced_BZ

    norm = 0.d0
    do i=1,Nk
       norm = norm +wt(i)
    end do

    if(dabs(norm-1.0).gt.1e-5) then
       write(*,*) norm,'error in BZ mesh'
    end if


  CONTAINS  !#########################################

    !BZ GRID
    SUBROUTINE build_reduced_BZ(off_set)
      integer            :: ix,iy,ik
      integer            :: weight
      integer            :: nk_bz
      logical,intent(in) :: off_set
      real(8)            :: x,y,dx,dy
      type(vec2D)        :: k
      real(8)            :: test

      nk_bz = 4*(nx-1)*(nx-1)

      allocate(kgrid(nx,nx),kindex(nx,nx),mp_kgrid(nx,nx))
      allocate(wt(Nk),vec_k(Nk),x_index(Nk),y_index(Nk))


      if(.not.off_set) then

         !write(*,*) 'no off set'
         if(nx.gt.1) then
            dx=1./dble(nx-1)
         else
            dx=0.d0
         end if
         x=-dx
         ik=0
         test = 0.d0
         do ix=1,nx
            x=x+dx
            y=-dx
            do iy=1,ix
               y=y+dx
               kgrid(ix,iy)=.5*x*bi+.5*y*bj
               ik=ik+1
               vec_k(ik) = .5*x*bi+.5*y*bj
               kindex(ix,iy)=ik
               x_index(ik) = ix
               y_index(ik) = iy
               !define the weight of each k point
               if(ix.eq.1) then
                  weight = 1.0            !center of BZ  (Gamma)
               elseif(ix.eq.nx) then
                  if(iy.eq.ix) then
                     weight=1.0           !corner of BZ (X)
                  elseif(iy.eq.1) then
                     weight=2.0           !M point (0,pi)
                  else
                     weight=4.0           !edge  (M->X)
                  end if
               else
                  if(iy.eq.ix) then
                     weight = 4.0         !diagonal (Gamma->x)
                  elseif(iy.eq.1) then
                     weight = 4.0         !horizontal cut (Gamma->M)
                  else
                     weight = 8.0         !BZ bulk
                  end if
               end if
               wt(ik) = dble(weight)/dble(nk_bz)
               test = test + wt(ik)
            end do
         end do

         !           write(*,*) ik

      else
         !           write(*,*) 'zero off_set'

         dx=1./dble(nx)
         x=-dx*.5
         ik=0
         test = 0.d0
         do ix=1,nx
            x=x+dx
            y=-dx*.5

            do iy=1,ix
               y=y+dx
               kgrid(ix,iy)=.5*x*bi+.5*y*bj

               ik=ik+1
               kindex(ix,iy)=ik
               vec_k(ik) = .5*x*bi+.5*y*bj              
               
               !define the weight of each k point
               if(ix.eq.iy) then
                  weight = 4.0
               else
                  weight = 8.0
               end if
               wt(ik) = dble(weight)/dble(nx)/dble(nx)/4!*0.5d0
               test = test + wt(ik)
            end do
         end do

         !           write(*,*) ik

      end if


      write(*,*) 'TEST 2d BZ',test

    END SUBROUTINE build_reduced_BZ


    !+----------------------------------------+!
    ! alternative choise : Monkhorst Pack grid !
    !+----------------------------------------+!
    SUBROUTINE get_monk_pack_grid
      integer :: i,j
      real(8) :: mp_x,mp_y

      do i=1,nx

         mp_x = (2.d0*i-nx-1.d0)/(2.d0*nx)

         do j=1,nx
            mp_y = (2.d0*j-nx-1.d0)/(2.d0*nx)
            mp_kgrid(i,j) = bi*mp_x + bj*mp_y
         end do

      end do

    END SUBROUTINE get_monk_pack_grid




    !nr of points in the reduced bz
    FUNCTION Npoints_BZ(nx) RESULT(NBZ)
      integer :: nx
      integer :: NBZ
      NBZ = nx*nx+nx
      NBZ=NBZ/2
      return
    END FUNCTION Npoints_BZ


    !print BZ function
    SUBROUTINE print_reduced_BZ(pm3d)
      integer, optional :: pm3d
      integer :: ix,iy,ik

      open(10,file='BZ_mesh.out')


      do ix=1,nx
         do iy=1,ix
            
            ik=kindex(ix,iy)
            write(10,"(3(F18.10),I4)") kgrid(ix,iy)%x,kgrid(ix,iy)%y,wt(ik),ik
            
         end do

         if(present(pm3d)) then
            write(10,*)
         end if
         
      end do
      

      close(10)
      close(20)


    END SUBROUTINE print_reduced_BZ




  END SUBROUTINE reduced_BZ


  !+--------------------------------+!
  !+- BUILD BATH DENSITY OF STATES -+!
  !+--------------------------------+!

  SUBROUTINE bath_dos(W_bath,N_bath,ibath)
    implicit none
    real(8)  :: W_bath
    real(8)  :: d_ene
    integer  :: N_bath
    integer  :: N_star
    integer  :: iw,ibath
    real(8)  :: test
    real(8)  :: beta

    beta = 10

    N_star = (N_bath+1)/2
    open(unit=50,file='bath_dos.out')

    allocate(wt_bath(N_bath),dos_bath(N_bath),ene_bath_log(N_bath),dos_bath_log(N_bath),wt_bath_log(N_bath),k_bath_log(N_bath),k_bath(N_bath))

    !+-------------------------------------------------+!
    !+- build up log-mesh for the one dimensional DOS -+!
    !+-------------------------------------------------+!
    write(*,*)"WARNING: there is an error here. log-mesh is wrong somewhere, we`re too lazy to fix it now..."
    allocate(ene_bath(N_star))
    ene_bath = logspace(W_bath,0.d0,N_star,10.d0)

    do iw = 1,N_star
       ene_bath_log(iw) = ene_bath(N_star+1-iw) - 2.d0
    end do

    do iw = 1,N_star-1
       ene_bath_log(iw+N_star) = ene_bath(iw+1) + (1.d0-ene_bath(iw+1))*2.d0
    end do

    test = 0.d0
    do iw = 1,N_bath
       dos_bath_log(iw) = 1.d0/sqrt(1-(ene_bath_log(iw)/2)**2)/2.d0/pi
       if(iw.gt.1.and.iw.lt.N_bath) then
          wt_bath_log(iw) = (ene_bath_log(iw+1) - ene_bath_log(iw-1))*0.5d0
       else
          if(iw.eq.1) then
             wt_bath_log(iw) = (ene_bath_log(iw+1) - ene_bath_log(iw))*0.5d0
          else
             wt_bath_log(iw) = (ene_bath_log(iw) - ene_bath_log(iw-1))*0.5d0
          end if

       end if

!!$         if(ene_bath_log(iw).le.0.d0) then
!!$            test = test + dos_bath_log(iw)*wt_bath_log(iw)*ene_bath_log(iw)
!!$         end if

       test = test + dos_bath_log(iw)*wt_bath_log(iw)

       k_bath_log(iw) = acos(-ene_bath_log(iw)/2.d0)

       write(50,"(4(F18.10))") ene_bath_log(iw),dos_bath_log(iw),wt_bath_log(iw),k_bath_log(iw)
    end do

    write(50,*)
    write(50,*)



    deallocate(ene_bath)


    !+----------------------------------------------------------------+!
    !+- build up the linear mesh for other DOS (flat,semicircular..) -+!
    !+----------------------------------------------------------------+!
    allocate(ene_bath(N_bath))

    ene_bath = linspace(-W_bath,W_bath,N_bath,.true.,mesh=d_ene)

    write(*,*) 'de_bath',d_ene

    test = 0.d0
    do iw = 1,N_bath

       if(iw.gt.1.and.iw.lt.N_bath) then
          wt_bath(iw) = (ene_bath(iw+1) - ene_bath(iw-1))*0.5d0
       else

          if(iw.eq.1) then
             wt_bath(iw) = (ene_bath(iw+1) - ene_bath(iw))*0.5d0
          else
             wt_bath(iw) = (ene_bath(iw) - ene_bath(iw-1))*0.5d0
          end if

       end if

       select case(ibath)

       case(0)
          dos_bath(iw) = 1.d0/2.d0/W_bath
       case(1)
          dos_bath(iw) = 2.d0/pi/W_bath/W_bath*sqrt(W_bath**2-ene_bath(iw)**2)
          !+- insert rounded density of states -+!
       end select


       test = test + dos_bath(iw)*wt_bath(iw)!*ene_bath_log(iw)


       k_bath(iw) = acos(-ene_bath(iw)/2.d0)

       write(50,*) ene_bath(iw),dos_bath(iw),wt_bath(iw)
    end do

    write(*,*) 'test',test,1./pi      

    close(50)

  END SUBROUTINE bath_dos




  FUNCTION bath_coupling(k,iw) result(vk)
    implicit none 
    real(8) :: k,vk
    integer :: iw

    !vk = sqrt(dos_bath(iw)*wt_bath(iw))
    vk = sqrt(dos_bath(iw)*wt_bath(iw))

    return
  END FUNCTION bath_coupling







  !+-------------------------------------------+!
  !--------DISPERSION IN SQUARE LATTICEs--------!
  !+-------------------------------------------+!

  FUNCTION chain_disp(k)
    implicit none
    real(8) :: chain_disp
    real(8) :: k
    chain_disp = -2.0d0*cos(k)
    return
  END FUNCTION chain_disp

  FUNCTION square_lattice_disp(k,shift)
    implicit none
    real(8) :: square_lattice_disp
    type(vec2D) :: k
    real(8),dimension(2),optional :: shift
    real(8),dimension(2)    :: shift_
    shift_=0.d0
    if(present(shift)) shift_=shift
    square_lattice_disp = -2*(cos(k%x-shift_(1))+cos(k%y-shift_(2)))
    return
  END FUNCTION square_lattice_disp



  FUNCTION cubic_lattice_disp(k)
    implicit none
    real(8) :: cubic_lattice_disp
    type(vec3D) :: k
    cubic_lattice_disp = -2*(dcos(k%x)+dcos(k%y)+dcos(k%z))
    return
  END FUNCTION cubic_lattice_disp


  !+- FERMI DISTRIBUTION FUNCTION T=0 -+!
  FUNCTION fermi_zero(ek,chem_pot) result(f)
    implicit none 
    real(8)   :: chem_pot
    real(8)   :: ek
    real(8)   :: f
    real(8)   :: beta

    beta = 1.d+3
    if(ek.gt.0.1/beta) then
       f = 0.d0
    else
       if(ek.lt.-0.1/beta) then
          f=1.d0
       else
          f = 1.d0/(exp(beta*ek)+1.d0)
       end if
    end if
    !
    if(ek.gt.0.d0) then
       f=0.d0
    else
       if(ek.eq.0.d0) then
          f=0.5d0
       else
          f=1.d0
       end if
    end if
    ! f = 1.d0/(exp(beta*ek)+1.d0)

    return
  END FUNCTION fermi_zero

END MODULE BZ_POINTS
