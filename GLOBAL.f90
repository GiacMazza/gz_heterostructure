!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!DEFINE GLOBAL VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE global
  USE SCIFOR
!  USE TIMER
!  USE CONSTANTS
!  USE MATRIX
!  USE PARSE_INPUT
!  USE TOOLS
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

  !# Electric field #!
  character(len=16)    :: field_type    !choose the profile of the electric field
  real(8)              :: Dpulse
  real(8)              :: Efield        !Electric field strength
  real(8),dimension(3) :: Evect         !Electric field vectors as input
  real(8)              :: Ton,Toff      !turn on/off time, t0 also center of the pulse
  integer              :: Ncycles       !Number of cycles in pulsed light packet
  real(8)              :: omega0        !parameter for the Oscilatting field and Pulsed light
  real(8)              :: E1            !Electric field strenght for the AC+DC case (tune to resonate)


  !hoppings (tslab is the unit of energy)
  real(8)                                     :: t_lead !
  real(8)                                     :: t_perp !
  real(8)                                     :: k_diss,beta_diss
  real(8),dimension(:,:),allocatable                                     :: Hslab !
  logical                                     :: Vhyb   !Hybridization potential
  real(8)                                     :: vL,vR
  real(8),dimension(:,:),allocatable          :: vk_L,vk_R
  real(8)                                     :: W_bath,Gamma
  !Hubbard repulsion
  real(8)                                     :: U
  real(8),dimension(:),allocatable            :: Uz         
  real(8),dimension(:,:),allocatable          :: Uz_time
  real(8),dimension(:,:,:),allocatable          ::   ek_time         


  !leads chemical potential!
  real(8)                                     :: mu_L,mu_R
  real(8)                                     :: beta_left,beta_right
  real(8)                                     :: mu_pulseL,mu_pulseR
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
  real(8)                                     :: cut_time,pulse_time

  real(8)                                     :: eta_bath
  real(8)                                     :: chem_equ


  !+------------------------------------------------+!
  ! k-grid details
  !+------------------------------------------------+!
  integer                                     :: nx_grid !nr of points in the x direction
  integer                                     :: Nprint
  integer                                     :: Nprint_nk
  integer                                     :: Nk_tot  !total nr of k-points
  logical                                     :: off_set !off_set at gamma point
  integer                                     :: Nk_orth !nr of k-points in the z-direction
  integer                                     :: ik0
  logical :: dos_plane  

  !+------------------------+!
  !  minimization algorithm  !
  !+------------------------+!
  integer                                     :: n_max         !maximum number of step
  integer                                     :: start_gz      !gz starting projectors
  real(8)                                     :: conv_treshold !convergence treshold
  real(8),dimension(:),allocatable                                     :: epsik
  integer(8),dimension(:),allocatable         :: ik_ord


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
       t_perp              , &
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
  subroutine read_input(INPUTunit)
    implicit none
    character(len=*) :: INPUTunit
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
    t_perp              = 0.5d0
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


    ! inquire(file="input_slab.in",exist=control)    
    ! if(control)then
    !    open(unit=10,file="input_slab.in")
    !    read(10,nml=variables)
    !    close(10)
    ! else
    !    print*,"Can not find INPUT file"
    !    print*,"Printing a default version in default.input_slab.in"
    !    open(50,file="default.input_slab.in")
    !    write(50,nml=variables)
    !    write(50,*)""
    !    stop
    ! endif




    call parse_input_variable(Slab                ,"L",INPUTunit,default=Slab)
    call parse_input_variable(U                   ,"U",INPUTunit,default=U)         
    call parse_input_variable(Temp                ,"TEMP",INPUTunit,default=Temp)         
    call parse_input_variable(Nx_grid             ,"NX_GRID",INPUTunit,default=Nx_grid)
    call parse_input_variable(Nk_orth             ,"NK_ORTH",INPUTunit,default=Nk_orth)       
    call parse_input_variable(off_set             ,"OFF_SET",INPUTunit,default=off_set)
    call parse_input_variable(pbc                 ,"PBC",INPUTunit,default=pbc)
    call parse_input_variable(start_gz            ,"START_GZ",INPUTunit,default=start_gz)
    call parse_input_variable(n_max               ,"N_MAX",INPUTunit,default=n_max)
    call parse_input_variable(conv_treshold       ,"CONV_TRESHOLD",INPUTunit,default=conv_treshold)
    call parse_input_variable(dt                  ,"DT",INPUTunit,default=dt)
    call parse_input_variable(Nt                  ,"NT",INPUTunit,default=Nt)    
    call parse_input_variable(mrk                 ,"MRK",INPUTunit,default=mrk)
    call parse_input_variable(vL                  ,"VL",INPUTunit,default=vL)       
    call parse_input_variable(vR                  ,"VR",INPUTunit,default=vR)        
    call parse_input_variable(t_lead              ,"T_LEAD",INPUTunit,default=t_lead)        
    call parse_input_variable(t_perp              ,"T_PERP",INPUTunit,default=t_perp)
    call parse_input_variable(eta_bath                  ,"eta_bath",INPUTunit,default=0.d0)               
    call parse_input_variable(Vhyb                ,"VHYB",INPUTunit,default=Vhyb)        
    call parse_input_variable(chem_shift          ,"CHEM_SHIFT",INPUTunit,default=chem_shift)
    call parse_input_variable(left                ,"LEFT",INPUTunit,default=left)
    call parse_input_variable(right               ,"RIGHT",INPUTunit,default=right)
    call parse_input_variable(init                ,"INIT",INPUTunit,default=init)
    call parse_input_variable(W_bath              ,"W_BATH",INPUTunit,default=W_bath)
    call parse_input_variable(i_bath              ,"BATH_TYPE",INPUTunit,default=i_bath)
    call parse_input_variable(kleads              ,"kLEADS",INPUTunit,default=kleads)
    call parse_input_variable(real_space          ,"REAL_SPACE",INPUTunit,default=real_space)
    call parse_input_variable(Gamma               ,"GAMMA",INPUTunit,default=Gamma)
    call parse_input_variable(alpha_electrostatic ,"ALPHA_ELECTROSTATIC",INPUTunit,default=alpha_electrostatic)  
    call parse_input_variable(ramp_time           ,"RAMP",INPUTunit,default=ramp_time)  
    call parse_input_variable(time_bias           ,"TIME_BIAS",INPUTunit,default=time_bias)  
    call parse_input_variable(linear_field        ,"LINEAR_FIELD",INPUTunit,default=linear_field)  
    call parse_input_variable(inner_field_type    ,"INNER_FIELD",INPUTunit,default=inner_field_type)  
    call parse_input_variable(lead_chem           ,"LEAD_CHEM",INPUTunit,default=lead_chem)  
    call parse_input_variable(lead_type           ,"LEAD_TYPE",INPUTunit,default=lead_type)  
    call parse_input_variable(dop_layer           ,"DOP_LAYER",INPUTunit,default=dop_layer)  
    call parse_input_variable(beta_left           ,"BETA_L",INPUTunit,default=1000.d0)  
    call parse_input_variable(beta_right           ,"BETA_R",INPUTunit,default=1000.d0)  
    call parse_input_variable(Nprint           ,"Nprint",INPUTunit,default=10)  
    call parse_input_variable(Nprint_Nk           ,"Nprint_NK",INPUTunit,default=100)  

    !ELECTRIC FIELD VARIABLES
    call parse_input_variable(field_type,"FIELD_TYPE",INPUTunit,default ='pulse',comment="profile type of the electric field ")
    call parse_input_variable(Efield,"EFIELD",INPUTunit,default=0d0,comment="electric field strength")
    call parse_input_variable(Evect,"EVECT",INPUTunit,default=[1d0,0d0,0d0],comment="electric field direction (normalized)")
    call parse_input_variable(ton,"TON",INPUTunit,default=0d0,comment="turn on time or center of the pulse")
    call parse_input_variable(toff,"TOFF",INPUTunit,default=10000d0,comment="turn off time")
    call parse_input_variable(Dpulse,"DPULSE",INPUTunit,default=10.d0,comment="time width of the field pulse")
    call parse_input_variable(omega0,"OMEGA0",INPUTunit,default=acos(-1d0) , comment="parameter for the Oscilatting field and Pulsed light")


    ! call parse_cmd_variable(Slab                ,"L")
    ! call parse_cmd_variable(U                   ,"U")         
    ! call parse_cmd_variable(Temp                ,"TEMP")         
    ! call parse_cmd_variable(Nx_grid             ,"NX_GRID")
    ! call parse_cmd_variable(Nk_orth             ,"NK_ORTH")       
    ! call parse_cmd_variable(off_set             ,"OFF_SET")
    ! call parse_cmd_variable(pbc                 ,"PBC")
    ! call parse_cmd_variable(start_gz            ,"START_GZ")
    ! call parse_cmd_variable(n_max               ,"N_MAX")
    ! call parse_cmd_variable(conv_treshold       ,"CONV_TRESHOLD")
    ! call parse_cmd_variable(dt                  ,"DT")
    ! call parse_cmd_variable(Nt                  ,"NT")    
    ! call parse_cmd_variable(mrk                 ,"MRK")
    ! call parse_cmd_variable(vL                  ,"VL")       
    ! call parse_cmd_variable(vR                  ,"VR")        
    ! call parse_cmd_variable(t_lead              ,"T_LEAD")        
    ! call parse_cmd_variable(t_perp              ,"T_PERP")        
    ! call parse_cmd_variable(Vhyb                ,"VHYB")        
    ! call parse_cmd_variable(chem_shift          ,"CHEM_SHIFT")
    ! call parse_cmd_variable(left                ,"LEFT")
    ! call parse_cmd_variable(right               ,"RIGHT")
    ! call parse_cmd_variable(init                ,"INIT")
    ! call parse_cmd_variable(W_bath              ,"W_BATH")
    ! call parse_cmd_variable(i_bath              ,"BATH_TYPE")
    ! call parse_cmd_variable(kleads              ,"kLEADS")
    ! call parse_cmd_variable(real_space          ,"REAL_SPACE")
    ! call parse_cmd_variable(Gamma               ,"GAMMA")
    ! call parse_cmd_variable(alpha_electrostatic ,"ALPHA_ELECTROSTATIC")  
    ! call parse_cmd_variable(ramp_time           ,"RAMP")  
    ! call parse_cmd_variable(time_bias           ,"TIME_BIAS")  
    ! call parse_cmd_variable(linear_field        ,"LINEAR_FIELD")  
    ! call parse_cmd_variable(inner_field_type    ,"INNER_FIELD")  
    ! call parse_cmd_variable(lead_chem           ,"LEAD_CHEM")  
    ! call parse_cmd_variable(lead_type           ,"LEAD_TYPE")  
    ! call parse_cmd_variable(dop_layer           ,"DOP_LAYER")  

    ! call parse_cmd_variable(beta_left           ,"BETA_L")  
    ! call parse_cmd_variable(beta_right           ,"BETA_R")  

    ! open(unit=10,file="used_input.out")
    ! write(10,nml=variables)
    ! close(10)

    beta = 1.d0/Temp

    tot_size = 2*N+L


    select case(lead_type)
    case('3d_tb')
       !t_lead = 1.d0
    case('generic_bath')
       t_lead = W_bath*0.5d0
    end select

    !Nk_tot = nx_grid*(nx_grid+1)/2

    ! if(.not.Vhyb) then
    !    L = Slab
    !    !allocate(Uz(L),Gslab_equ(Nk_tot,L,L),eqPhi(L))
    !    allocate(Uz(L))
    !    Uz = U
    !    Uz(1) = U*1.d0
    ! else
    !    if(.not.real_space) Nt=0
    !    L = Slab + 2*Nk_orth
    !    !allocate(Uz(L),Gslab_equ(Nk_tot,L,L),eqPhi(L))
    !    allocate(Uz(L))
    !    Uz(1:Nk_orth) = 0.d0
    !    Uz(Nk_orth+1:Nk_orth+Slab) = U
    !    Uz(Nk_orth+Slab+1:L) = 0.d0
    !    do ik0=1,L
    !       write(*,*) Uz(ik0),U
    !    end do
    ! end if

    
    eta_bath=0.d0

  END SUBROUTINE read_input

END MODULE global
