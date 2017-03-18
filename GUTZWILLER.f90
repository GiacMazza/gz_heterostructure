!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Compute stuff on the GUTZWILLER wfc in  the slab geometry
!Author: G Mazza
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE GUTZWILLER
  implicit none
  private

!#############################################################################################!
!(...work in prog...NEW MODULE)

  complex(8),parameter                                   :: Z1 = (1.d0,0.d0)
  complex(8),parameter                                   :: Z0 = (0.d0,0.d0)
  complex(8),parameter                                   :: Zi = (0.d0,1.d0)
  

  real(8),parameter  :: tiny=0.d0
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

  public :: GZ_hop_hd
  public :: GZ_doping_hd
  public :: deriveR_doublon
  public :: deriveR_holon
  public :: deriveR_ph_doublon
  public :: deriveR_ph_holon

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

  !##################################################################!
  !           COMPUTATION WITH HOLONS AND DOUBLONS                   !
  !##################################################################!

  elemental function GZ_hop_hd(doublon,holon,ph_doublon,ph_holon)
    real(8),intent(in) :: doublon,holon,ph_doublon,ph_holon 
    complex(8) :: GZ_hop_hd
    real(8) :: x,p1
    x=GZ_doping_hd(doublon,holon,ph_doublon,ph_holon)
    p1 = 1.d0-holon-doublon
    GZ_hop_hd = sqrt(2.d0*p1/(1-x*x))*(sqrt(abs(holon))*Exp(Zi*ph_holon)+sqrt(abs(doublon))*Exp(-Zi*ph_doublon))
    return
  end function GZ_hop_hd


  elemental function GZ_doping_hd(doublon,holon,ph_doublon,ph_holon)
    real(8),intent(in) :: doublon,holon,ph_doublon,ph_holon 
    complex(8) :: GZ_doping_hd
    real(8) :: x,n0
    GZ_doping_hd = (holon - doublon)
    return
  end function GZ_doping_hd


  elemental function deriveR_doublon(doublon,holon,ph_doublon,ph_holon) result(dr_dd)
    real(8),intent(in) :: doublon,holon,ph_doublon,ph_holon 
    complex(8) :: dr_dd
    real(8) :: x,p1
    complex(8) :: R
    x  = GZ_doping_hd(doublon,holon,ph_doublon,ph_holon)
    p1 = 1.d0 - holon - doublon
    R = GZ_hop_hd(doublon,holon,ph_doublon,ph_holon)
    !
    dr_dd = -x/(1.d0-x*x)*R + sqrt(2.d0/(1.d0-x*x))* &
         ( -dsqrt(dabs(holon))*Exp(Zi*ph_holon)/2.d0/dsqrt(dabs(p1)) + ( 1.d0-2.d0*doublon-holon )*Exp(-Zi*ph_doublon)/2.d0/dsqrt(dabs(doublon*p1)) )
    
    ! dr_dd = -x/(1.d0-x*x)*R + sqrt(2.d0/(1.d0-x*x))* &
    !      ( -holon*Exp(Zi*ph_holon)/2.d0/(sqrt(holon*p1)) + ( 1.d0-2.d0*doublon-holon )*Exp(-Zi*ph_doublon)/2.d0/(sqrt(doublon*p1)) )
    !
  end function deriveR_doublon


  elemental function deriveR_holon(doublon,holon,ph_doublon,ph_holon) result(dr_dh)
    real(8),intent(in) :: doublon,holon,ph_doublon,ph_holon 
    complex(8) :: dr_dh
    real(8) :: x,p1
    complex(8) :: R
    x  = GZ_doping_hd(doublon,holon,ph_doublon,ph_holon)
    p1 = 1.d0 - holon - doublon
    R = GZ_hop_hd(doublon,holon,ph_doublon,ph_holon)
    ! dr_dh = x/(1.d0 - x*x)*R + sqrt(2.d0/(1.d0 - x*x))* &
    !      ( (1.d0 - 2.d0*holon-doublon)*Exp(Zi*ph_holon)/2.d0/(sqrt(holon*p1)) - doublon*Exp(-Zi*ph_doublon)/2.d0/(sqrt(doublon*p1)) )
    dr_dh = x/(1.d0 - x*x)*R + sqrt(2.d0/(1.d0 - x*x))* &
         ( (1.d0 - 2.d0*holon-doublon)*Exp(Zi*ph_holon)/2.d0/dsqrt(dabs(holon*p1)) - dsqrt(dabs(doublon))*Exp(-Zi*ph_doublon)/2.d0/dsqrt(p1) )
  end function deriveR_holon


  elemental function deriveR_ph_doublon(doublon,holon,ph_doublon,ph_holon) result(dr_dpd)
    real(8),intent(in) :: doublon,holon,ph_doublon,ph_holon 
    complex(8) :: dr_dpd
    real(8) :: x,p1
    complex(8) :: R
    x  = GZ_doping_hd(doublon,holon,ph_doublon,ph_holon)
    p1 = 1.d0 - holon - doublon
    R = GZ_hop_hd(doublon,holon,ph_doublon,ph_holon)
    dr_dpd =  -Zi*sqrt(2.d0*abs(doublon*p1)/(1.d0-x**2))*Exp(-Zi*ph_doublon)
    !write(*,*) dr_dpd
  end function deriveR_ph_doublon


  elemental function deriveR_ph_holon(doublon,holon,ph_doublon,ph_holon) result(dr_dph)
    real(8),intent(in) :: doublon,holon,ph_doublon,ph_holon 
    complex(8) :: dr_dph
    real(8) :: x,p1
    complex(8) :: R
    x  = GZ_doping_hd(doublon,holon,ph_doublon,ph_holon)
    p1 = 1.d0 - holon - doublon
    R = GZ_hop_hd(doublon,holon,ph_doublon,ph_holon)
    dr_dph = Zi*sqrt(2.d0*abs(holon*p1)/(1.d0-x**2))*Exp(Zi*ph_holon)    
  end function deriveR_ph_holon


  
  
END MODULE GUTZWILLER
