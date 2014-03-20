!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Compute stuff on the GUTZWILLER wfc in  the slab geometry
!Author: G Mazza
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE GUTZWILLER
  implicit none
  private

!#############################################################################################!
!(...work in prog...NEW MODULE)

  
  !+--------------------------------------------------------------+!
  !          define type containing GZ PROJECTORS                  !
  !+--------------------------------------------------------------+!
  type,public :: gz_projector
     complex(8) :: p0,p1,p2   !projectors
     complex(8) :: H0,H1,H2   !hamiltonian matrix elements
  end type gz_projector
  
  public :: GZ_init0
  public :: GZ_init_half
  public :: GZ_doping
  public :: GZ_hop
  public :: GZ_double
  public :: GZ_normalization
  public :: slab_hubb_ene
  

CONTAINS

  !+----------------------------------------------------+!
  ! allocate GZ projectors and initialize to bulk values !
  !+----------------------------------------------------+!
!  SUBROUTINE allocate_GZ(Phi,L)
!    implicit none
!    integer :: L
!    integer :: start_gz
!    type(gz_projector),dimension(:),allocatable :: Phi
!    allocate(Phi(L))
!    call Gz_init_half(Phi)
!  END SUBROUTINE allocate_GZ

  !+---------------------------------------------------+!
  !          INITIALIZE GUTZWILLER PROJECTORS           !
  !+---------------------------------------------------+!
  !--> not interacting values
  elemental function GZ_init0(dop)
    implicit none
    type(gz_projector) :: GZ_init0
    real(8),intent(in) :: dop
    GZ_init0%p0 = (1+dop)/2
    GZ_init0%p1 = (1-dop*dop)/2
    GZ_init0%p2 = (1-dop)/2
    return
  end function GZ_init0

  !--> bulk interacting values
  elemental function GZ_init_half(Ulayer)
    implicit none
    type(gz_projector) :: GZ_init_half
    real(8),intent(in) :: Ulayer
    real(8)            :: Uc
    Uc=16.d0
    if(Ulayer.lt.Uc) then
       GZ_init_half%p0 = sqrt(1-Ulayer/Uc)/2.d0
       GZ_init_half%p2 = sqrt(1-Ulayer/Uc)/2.d0
       GZ_init_half%p1 = sqrt(1+Ulayer/Uc)/sqrt(2.d0)
    else
       GZ_init_half%p0 = sqrt(1-.95)/2.d0
       GZ_init_half%p2 = sqrt(1-.95)/2.d0
       GZ_init_half%p1 = sqrt(1+.95)/sqrt(2.d0)
    end if
    return
  end function GZ_init_half


  !+--------------------------------------------------------------+!
  !------------- COMPUTE STUFF USING GZ PROJECTORS ----------------!
  !+--------------------------------------------------------------+!
    
  !---> number of particles
  elemental function GZ_doping(Phi)
        type(gz_projector),intent(in) :: Phi
    real(8) :: GZ_doping
    GZ_doping = ABS(Phi%p0)*ABS(Phi%p0)
    GZ_doping = GZ_doping - ABS(Phi%p2)*ABS(Phi%p2)
    return
  end function GZ_doping
    
  !---> hopping renormalization factors 
  elemental function GZ_hop(Phi)
    type(gz_projector),intent(in) :: Phi
    complex(8) :: GZ_hop
    real(8) :: x,n0
    x=GZ_doping(Phi)
    n0=(1-x)*.5d0
    
    GZ_hop = conjg(Phi%p0)*Phi%p1
    GZ_hop = GZ_hop + conjg(Phi%p1)*Phi%p2

    GZ_hop = GZ_hop*sqrt(2.d0/(1-x*x))

    
    return
  end function GZ_hop


  !---> double occupancy  
  elemental function GZ_double(Phi)
    type(gz_projector),intent(in) :: Phi
    real(8) :: GZ_double
    GZ_double = ABS(Phi%p0)*ABS(Phi%p0)
    GZ_double = GZ_double + ABS(Phi%p2)*ABS(Phi%p2)
    GZ_double = GZ_double*.5
    return
  end function GZ_double
  
  !---> GZ normalization
  elemental function GZ_normalization(Phi)
    type(gz_projector), intent(in) :: Phi
    real(8)  :: GZ_normalization
    GZ_normalization = ABS(Phi%p0)*ABS(Phi%p0)
    GZ_normalization = GZ_normalization + ABS(Phi%p1)*ABS(Phi%p1)
    GZ_normalization = GZ_normalization + ABS(Phi%p2)*ABS(Phi%p2)
    return
  end function GZ_normalization

  !---> Hubbard energy
  function slab_hubb_ene(Phi,L,Ulayer)
    implicit none
    type(gz_projector),dimension(L),intent(in) :: Phi
    real(8),dimension(L),intent(in) :: Ulayer
    real(8),dimension(L) :: d_occ
    real(8) :: slab_hubb_ene
    integer :: i,L

    d_occ = GZ_double(Phi)
    slab_hubb_ene = 0.d0
    do i=1,L
       slab_hubb_ene = slab_hubb_ene + d_occ(i)*Ulayer(i)
    end do
    return
  end function slab_hubb_ene
  
  
END MODULE GUTZWILLER
