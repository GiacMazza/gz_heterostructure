!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!DEFINE GLOBAL VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE global
  USE SCIFOR
  USE TIMER
  USE COMMON_VARS
  USE MATRIX
  USE PARSE_INPUT
  USE TOOLS
  !local modules
  USE VECTORS
  USE BZ_POINTS
  USE GUTZWILLER                !derived type for GZ, auxiliary function in:projector out: observables 

  !+------------------------------------------------+!
  ! slab details
  !+------------------------------------------------+!
  integer,public                              :: L    !number of slab layers
  integer,public                              :: Slab !..I have no idea which explain the choise of this variable name..
  integer                                     :: N    !number of leads layers
  integer                                     :: tot_size !total number of layers
  logical                                     :: pbc
  complex(8),dimension(:,:,:),allocatable     :: Gslab_equ


  !+----- EQUILIBRIUM TEMPERATURE -----+!
  real(8) :: Temp,beta
  real(8),dimension(:),allocatable  :: hop_intra
  real(8),dimension(:),allocatable  :: tmp_dop
  complex(8),dimension(:),allocatable  :: delta_plus,delta_minus
  real(8)  :: test_entropy
  type(gz_projector),allocatable,dimension(:)  ::  phi_out
  real(8),allocatable,dimension(:) :: mu_star
  real(8),allocatable,dimension(:) :: local_field
  real(8),dimension(:),allocatable :: mu_met,mu_ins  


  !hoppings (tslab is the unit of energy)
  real(8)                                     :: t_lead !
  logical                                     :: Vhyb   !Hybridization potential
  real(8)                                     :: vL,vR
  real(8),dimension(:,:),allocatable          :: vk_L,vk_R
  real(8)                                     :: W_bath,Gamma
  !Hubbard repulsion
  real(8)                                     :: U
  real(8),dimension(:),allocatable            :: Uz         
  real(8),dimension(:,:),allocatable          :: Uz_time         


  !leads chemical potential!
  real(8)                                     :: mu_L,mu_R
  real(8),dimension(:),allocatable            :: muL,muR
  real(8),dimension(:,:),allocatable          :: e_loc
  real(8)                                     :: chem_shift
  real(8)                                     :: switch_time
  logical                                     :: shift_slab
  logical                                     :: left,right
  logical                                     :: kleads
  logical                                     :: real_space
  real(8)                                     :: w_hyb,alpha
  
  real(8)                                     :: dop_layer

  real(8)                                     :: ramp_time,time_bias,ramp_bias
  real(8)                                     :: cut_time 


  !+------------------------------------------------+!
  ! k-grid details
  !+------------------------------------------------+!
  integer                                     :: nx_grid !nr of points in the x direction
  integer                                     :: Nk_tot  !total nr of k-points
  logical                                     :: off_set !off_set at gamma point
  integer                                     :: Nk_orth !nr of k-points in the z-direction
  integer                                     :: ik0
  

  !+------------------------+!
  !  minimization algorithm  !
  !+------------------------+!
  integer                                     :: n_max         !maximum number of step
  integer                                     :: start_gz      !gz starting projectors
  real(8)                                     :: conv_treshold !convergence treshold


  !+-----------------------+!
  !  GUTZWILLER PROJECTORS  !
  !+-----------------------+!
  type(gz_projector),dimension(:),allocatable :: eqPhi  !gz projectors at equilibrium
  
  !+----------+!
  !  DYNAMICS  !
  !+----------+!
  real(8)                                     ::  dt     ! time step
  real(8),dimension(:),allocatable            ::  t_grid,t_finer
  integer                                     ::  Nt     ! Number of time steps
  integer                                     ::  mrk    ! Rk order
  
  real(8)                                     ::  adiabatic
  real(8)                                     ::  init
  
  real(8)                                     ::  alpha_electrostatic
  logical                                     ::  linear_field,lead_chem           
  character(len=16)                           ::  inner_field_type
  character(len=16)                           ::  i_bath
  character(len=16)                           ::  lead_type

  !+-------------------------------+!
  !+- OTHER VERY USEFUL VARIABLES -+!
  !+-------------------------------+!
  complex(8),parameter                                   :: Z1 = (1.d0,0.d0)
  complex(8),parameter                                   :: Z0 = (0.d0,0.d0)
  complex(8),parameter                                   :: Zi = (0.d0,1.d0)


  !+-----------+!
  !  NAMELISTS  !
  !+-----------+!
  namelist/variables/        &
       Slab                , &
       Temp                , &
       t_lead              , &
       Vhyb                , &
       U                   , &
       nx_grid             , &
       Nk_orth             , &
       off_set             , &
       pbc                 , &
       start_gz            , &
       n_max               , &
       conv_treshold       , & 
       dt                  , &
       Nt                  , &
       mrk                 , &
       chem_shift          , &
       left                , &
       right               , &
       init                , &
       vL                  , &
       vR                  , &
       W_bath              , &
       Gamma               , &
       i_bath              , &
       lead_type           , &
       kleads              , &
       real_space          , &
       alpha_electrostatic , &
       inner_field_type    , &
       lead_chem           , &
       ramp_time           , &
       ramp_bias           , &
       time_bias           , &
       dop_layer
       

CONTAINS

  !+---------------------+!
  !    READ INPUT DATA    !
  !+---------------------+!
  subroutine read_input
    implicit none
    logical :: control
    integer :: i,iU


    !DEFAULT INPUT FILE:
    SLab                = 4
    U                   = 5.d0
    Temp                = 1.d-3
    Nx_grid             = 8
    Nk_orth             = 150
    Vhyb                = .false.
    t_lead              = 1.d0
    off_set             = .false.
    pbc                 = .false.
    start_gz            = 0
    n_max               = 5000
    conv_treshold       = 1.d-6
    dt                  = 0.1d0
    Nt                  = 50
    mrk                 = 4
    vL                  = 0.5d0
    vR                  = 0.5d0
    chem_shift          = 0.5d0
    left                = .true.
    right               = .true.
    init                = 1.d-3
    W_bath              = 10.d0
    kleads              = .false.
    real_space          = .false.
    Gamma               = 1.d-2
    alpha_electrostatic = 0.0
    linear_field        = .false.
    lead_chem           = .true.
    ramp_time           = 0.d0
    time_bias           = ramp_time
    ramp_bias           = ramp_time
    cut_time            = 1000
    i_bath              = 'flat'
    lead_type           = 'generic_bath'
    inner_field_type    = 'none'
    shift_slab          = .true.
    dop_layer           = 0.d0


    inquire(file="input_slab.in",exist=control)    
    if(control)then
       open(unit=10,file="input_slab.in")
       read(10,nml=variables)
       close(10)
    else
       print*,"Can not find INPUT file"
       print*,"Printing a default version in default.input_slab.in"
       open(50,file="default.input_slab.in")
       write(50,nml=variables)
       write(50,*)""
       stop
    endif

    call parse_cmd_variable(Slab                ,"L")
    call parse_cmd_variable(U                   ,"U")         
    call parse_cmd_variable(Temp                ,"TEMP")         
    call parse_cmd_variable(Nx_grid             ,"NX_GRID")
    call parse_cmd_variable(Nk_orth             ,"NK_ORTH")       
    call parse_cmd_variable(off_set             ,"OFF_SET")
    call parse_cmd_variable(pbc                 ,"PBC")
    call parse_cmd_variable(start_gz            ,"START_GZ")
    call parse_cmd_variable(n_max               ,"N_MAX")
    call parse_cmd_variable(conv_treshold       ,"CONV_TRESHOLD")
    call parse_cmd_variable(dt                  ,"DT")
    call parse_cmd_variable(Nt                  ,"NT")    
    call parse_cmd_variable(mrk                 ,"MRK")
    call parse_cmd_variable(vL                  ,"VL")       
    call parse_cmd_variable(vR                  ,"VR")        
    call parse_cmd_variable(t_lead              ,"T_LEAD")        
    call parse_cmd_variable(Vhyb                ,"VHYB")        
    call parse_cmd_variable(chem_shift          ,"CHEM_SHIFT")
    call parse_cmd_variable(left                ,"LEFT")
    call parse_cmd_variable(right               ,"RIGHT")
    call parse_cmd_variable(init                ,"INIT")
    call parse_cmd_variable(W_bath              ,"W_BATH")
    call parse_cmd_variable(i_bath              ,"BATH_TYPE")
    call parse_cmd_variable(kleads              ,"kLEADS")
    call parse_cmd_variable(real_space          ,"REAL_SPACE")
    call parse_cmd_variable(Gamma               ,"GAMMA")
    call parse_cmd_variable(alpha_electrostatic ,"ALPHA_ELECTROSTATIC")  
    call parse_cmd_variable(ramp_time           ,"RAMP")  
    call parse_cmd_variable(time_bias           ,"TIME_BIAS")  
    call parse_cmd_variable(linear_field        ,"LINEAR_FIELD")  
    call parse_cmd_variable(inner_field_type    ,"INNER_FIELD")  
    call parse_cmd_variable(lead_chem           ,"LEAD_CHEM")  
    call parse_cmd_variable(lead_type           ,"LEAD_TYPE")  
    call parse_cmd_variable(dop_layer           ,"DOP_LAYER")  

    open(unit=10,file="used_input.out")
    write(10,nml=variables)
    close(10)

    beta = 1.d0/Temp

    tot_size = 2*N+L

    mu_L = 0.d0
    mu_R = 0.d0

    cut_time = 100000.d0

    if(left)  mu_L= chem_shift*0.5d0
    if(right) mu_R=-chem_shift*0.5d0

    select case(lead_type)
    case('3d_tb')
       !t_lead = 1.d0
    case('generic_bath')
       t_lead = W_bath*0.5d0
    end select

    Nk_tot = nx_grid*(nx_grid+1)/2

    if(.not.Vhyb) then
       L = Slab
       allocate(Uz(L),Gslab_equ(Nk_tot,L,L),eqPhi(L))
       Uz = U
       Uz(1) = U*1.d0
    else
       if(.not.real_space) Nt=0
       L = Slab + 2*Nk_orth
       allocate(Uz(L),Gslab_equ(Nk_tot,L,L),eqPhi(L))
       Uz(1:Nk_orth) = 0.d0
       Uz(Nk_orth+1:Nk_orth+Slab) = U
       Uz(Nk_orth+Slab+1:L) = 0.d0
       do ik0=1,L
          write(*,*) Uz(ik0),U
       end do
    end if


  END SUBROUTINE read_input

END MODULE global
