!----------------------------------------------!
!DEFINE VECTORS STRUCTURES
!Author G Mazza
!Crucial hints A Amaricci
!----------------------------------------------!
MODULE VECTORS
  implicit none
  private

  type,public :: vec2D
     real(8) :: x
     real(8) :: y
  end type vec2D

  type,public :: vec3D
     real(8) :: x
     real(8) :: y
     real(8) :: z
  end type vec3D

  interface operator(+)
     module procedure  add
  end interface

  interface operator(-)
     module procedure subtract
  end interface

  interface operator(*)
     module procedure prodL, prodR
  end interface
  
  interface operator(.dot.)
     module procedure dot
  end interface
  
  interface assignment(=)
     module procedure vector_eq, scalar_eq
  end interface

  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(.dot.)
  public :: assignment(=)
  public :: modulo,add,subtract,prodL,prodR,dot,scalar_eq,vector_eq

  !some special vectors

  type(vec2D),parameter,public    ::  Xver = vec2D(1.d0,0.0)
  type(vec2D),parameter,public    ::  Yver = vec2D(0.0,1.d0)
  type(vec2D),parameter,public    ::  Uvec = vec2D(1.d0,1.d0)
  type(vec2D),parameter,public    ::  Vzero = vec2D(0.d0,0.d0)


CONTAINS
  
  !+---------------------------------------------+!
  !             define some algebra               !
  !+---------------------------------------------+!

  elemental function add(V,W)
    type(vec2D),intent(in)  :: V,W
    type(vec2D) :: add
    add%x = V%x + W%x
    add%y = V%y + W%y
  end function add
  !-------------------------------------

  !+-----------------------------------------------------------------+
  elemental function modulo(V)
    type(vec2D),intent(in)  :: V
    real(8) :: modulo
    modulo = sqrt(V%x**2 + V%y**2)
  end function modulo


  !-------------------------------------
  elemental function subtract(V,W)
    type(vec2D),intent(in)  :: V,W
    type(vec2D) :: subtract
    subtract%x = V%x - W%x
    subtract%y = V%y - W%y
  end function subtract
  !-------------------------------------

  !-------------------------------------
  elemental function prodL(C,V)
    real(8),intent(in) :: C
    type(vec2D),intent(in)  :: V
    type(vec2D) :: prodL
    prodL%x = C*V%x 
    prodL%y = C*V%y
  end function prodL
  !-------------------------------------

  !-------------------------------------
  elemental function prodR(V,C)
    real(8),intent(in) :: C
    type(vec2D),intent(in)  :: V
    type(vec2D) :: prodR
    prodR%x = V%x*C 
    prodR%y = V%y*C
  end function prodR
  !-------------------------------------

  !-------------------------------------
  elemental function dot(V,W)
    type(vec2D),intent(in)  :: V,W
    real(8) :: dot
    dot = V%x*W%x + V%y*W%y
  end function dot
  !-------------------------------------

  !-------------------------------------
  elemental subroutine scalar_eq(V,C)
    real(8),intent(in)  :: C
    type(vec2D),intent(inout) :: V
    V%x = C
    V%y = C
  end subroutine  scalar_eq
  !-------------------------------------

  !-------------------------------------
  elemental subroutine vector_eq(V,W)
    type(vec2D),intent(in)  :: W
    type(vec2D),intent(inout) :: V
    V%x = W%x
    V%y = W%y
  end subroutine  vector_eq


END MODULE VECTORS

