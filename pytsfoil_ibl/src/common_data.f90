! common_data.f90
! Global data:
!   - Constants
!   - User-input parameters
!   - Mesh and geometry parameters and arrays
!   - File unit numbers

module common_data
    implicit none
    
    ! ------------------------------------------------
    ! Constants
    ! ------------------------------------------------
    
    integer, parameter :: N_MESH_POINTS = 1000            ! Mesh size parameter - change this to adjust mesh dimensions
    integer, parameter :: NMP_plus2 = N_MESH_POINTS + 2   ! Number of mesh points + 2
    integer, parameter :: NMP_plus1 = N_MESH_POINTS + 1   ! Number of mesh points + 1

    integer, parameter :: IMIN = 1
    integer, parameter :: JMIN = 1

    real, parameter :: GAM = 1.4            ! Specific heat ratio
    real, parameter :: GAM1 = GAM + 1.0     ! gamma + 1

    ! ------------------------------------------------
    ! User-input parameters
    ! ------------------------------------------------

    real :: CLSET = 0.0     ! Lift coefficient set-point

    integer :: BCTYPE = 1   ! Boundary condition identifiers (1 = free air, 2 = tunnel)

    real :: XL(N_MESH_POINTS) = 0.0, XU(N_MESH_POINTS) = 0.0 ! Airfoil X-coordinates
    real :: YL(N_MESH_POINTS) = 0.0, YU(N_MESH_POINTS) = 0.0 ! Airfoil Y-coordinates
    real :: DELTA = 0.0       ! Maximum thickness of airfoil (set to zero to raise error)

    logical :: KUTTA = .true. ! Whether Kutta condition is enforced
    logical :: FCR = .true.   ! Whether difference equations are fully conservative

    integer :: SIMDEF = 3   ! Similarity scaling (1 = Cole, 2 = Spreiter, 3 = Krupp)

    real :: AK = 0.0        ! Free stream similarity parameter
    real :: RIGF = 0.0      ! Rigidity factor for transonic effects
    real :: POR = 0.0       ! Porosity

    real :: WCIRC = 1.0         ! Weight for circulation jump at trailing edge (0.0-1.0)

    real :: EPS = 0.2           ! Convergence tolerance 
    real :: CVERGE = 0.00001    ! Error criterion for convergence
    real :: DVERGE = 10.0       ! Error criterion for divergence

    ! Wall/tunnel constants (Optional)
    real :: F = 0.0
    real :: H = 0.0

    ! Flap parameters (Optional)
    integer :: IFLAP = 0      ! Flap flag
    real :: DELFLP = 0.0      ! Flap deflection angle  
    real :: FLPLOC = 0.77     ! Flap location


    ! ------------------------------------------------
    ! Mesh and geometry parameters and arrays
    ! ------------------------------------------------
    ! Mesh indices
    integer :: IMAX, JMAX   ! maximum number of grid points in X/Y-direction used in code
    integer :: IUP, IDOWN   ! upstream/downstream indices
    integer :: ILE, ITE     ! leading/trailing edge i-indices
    integer :: JUP          ! upper surface j-indices, index of first point where Y > 0.0 (calculated by JSLIT)
    integer :: JLOW         ! lower surface j-indices, JLOW = JUP - 1 (calculated by JSLIT)
    integer :: JTOP, JBOT   ! far-field top/bottom j-indices

    ! Mesh coordinate arrays
    real :: X(NMP_plus2) = 0.0, Y(NMP_plus2) = 0.0  ! Mesh coordinate arrays
    real :: XDIFF(N_MESH_POINTS) = 0.0, YDIFF(N_MESH_POINTS) = 0.0 ! mesh derivative arrays

    ! Airfoil arrays
    real :: VOL = 0.0
    real :: FU(N_MESH_POINTS) = 0.0 ! Upper surface interpolation
    real :: FL(N_MESH_POINTS) = 0.0 ! Lower surface interpolation
    real :: FXU(N_MESH_POINTS) = 0.0 ! Derivative of upper surface to X-coordinate
    real :: FXL(N_MESH_POINTS) = 0.0 ! Derivative of lower surface to X-coordinate
    
contains

    ! Initialize common data arrays and parameters
    subroutine initialize_common()
        implicit none

        IMAX = N_MESH_POINTS
        JMAX = N_MESH_POINTS
        
        ! Initialize mesh indices to safe defaults (will be recalculated later)
        IUP = 2
        IDOWN = IMAX - 1
        ILE = IMIN + 5  ! Safe default
        ITE = IMAX - 5  ! Safe default
        JUP = (JMAX + JMIN) / 2 + 1   ! Safe default above center
        JLOW = (JMAX + JMIN) / 2 - 1  ! Safe default below center
        JTOP = JMAX - 1
        JBOT = JMIN + 1

        ! ------------------------------------------------
        ! Initialize all common data
        ! ------------------------------------------------
        
        CLSET = 0.0
    
        BCTYPE = 1
    
        XL = 0.0
        XU = 0.0
        YL = 0.0
        YU = 0.0
        DELTA = 0.0
            
        KUTTA = .true.
        FCR = .true.
    
        SIMDEF = 3
    
        AK = 0.0
        RIGF = 0.0
        POR = 0.0
    
        WCIRC = 1.0
    
        EPS = 0.2
        CVERGE = 0.00001
        DVERGE = 10.0
        F = 0.0
        H = 0.0
    
        IFLAP = 0
        DELFLP = 0.0
        FLPLOC = 0.77

        X = 0.0
        Y = 0.0
        XDIFF = 0.0
        YDIFF = 0.0
    
        VOL = 0.0

    end subroutine initialize_common

end module common_data
