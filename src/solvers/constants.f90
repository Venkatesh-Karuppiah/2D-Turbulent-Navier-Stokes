!===============================================================================
! constants.f90 - Physical and numerical constants for Navier-Stokes solver
!===============================================================================
MODULE CONSTANTS
    IMPLICIT NONE

    ! Maximum problem size parameters
    INTEGER, PARAMETER :: NN = 60000   ! Maximum number of nodes/elements
    INTEGER, PARAMETER :: NB = 6000    ! Maximum boundary faces

    ! Gas dynamics constants
    DOUBLE PRECISION, PARAMETER :: RUNIV_VAL = 8314.0D0   ! Universal gas constant [J/(kmol*K)]

    ! Turbulence model constants (Baldwin-Lomax style)
    DOUBLE PRECISION, PARAMETER :: CMU_VAL    = 0.09D0
    DOUBLE PRECISION, PARAMETER :: ALPHA_TURB = 0.00052D0
    DOUBLE PRECISION, PARAMETER :: BETA_TURB  = 0.2D0
    DOUBLE PRECISION, PARAMETER :: CONST1_VAL = 21.479D0
    DOUBLE PRECISION, PARAMETER :: CONST2_VAL = 21.394D0
    DOUBLE PRECISION, PARAMETER :: CONST3_VAL = 1.5D0

    ! Sutherland's law coefficient
    DOUBLE PRECISION, PARAMETER :: SUTH_COEF = 0.00000011848D0

    ! Dissipation scheme constants (Jameson's scheme)
    DOUBLE PRECISION, PARAMETER :: C02_VAL = 0.25D0
    DOUBLE PRECISION, PARAMETER :: C04_VAL = 0.00390625D0   ! 1/256

    ! RK4 coefficients (Jameson's 4-stage scheme)
    DOUBLE PRECISION, PARAMETER :: ALPHA1_VAL = 0.25D0
    DOUBLE PRECISION, PARAMETER :: ALPHA2_VAL = 1.0D0/3.0D0
    DOUBLE PRECISION, PARAMETER :: ALPHA3_VAL = 0.5D0
    DOUBLE PRECISION, PARAMETER :: ALPHA4_VAL = 1.0D0

    ! Turbulent transport
    DOUBLE PRECISION, PARAMETER :: CFLTURB_VAL = 0.5D0

    ! Divisor constants
    DOUBLE PRECISION, PARAMETER :: DIV23_VAL = 2.0D0/3.0D0
    DOUBLE PRECISION, PARAMETER :: QUARTER   = 0.25D0
    DOUBLE PRECISION, PARAMETER :: HALF      = 0.5D0

    ! Convergence monitoring
    DOUBLE PRECISION, PARAMETER :: CDIFF1_INIT = 1.0E-10
    DOUBLE PRECISION, PARAMETER :: CDIFF5_INIT = 1.0E-20

    ! Boundary layer profile exponent
    DOUBLE PRECISION, PARAMETER :: BL_EXPONENT = 1.0D0/7.0D0

    ! Unit conversion
    DOUBLE PRECISION, PARAMETER :: MM_TO_M   = 0.001D0     ! Millimeters to meters
    DOUBLE PRECISION, PARAMETER :: PA_TO_BAR = 0.00001D0   ! Pascal to Bar

    ! Combustion reference
    DOUBLE PRECISION, PARAMETER :: CONVKSC_VAL = 1.0D0/98066.5D0
    DOUBLE PRECISION, PARAMETER :: GTH_VAL     = 35.0D0

    ! Threshold values
    DOUBLE PRECISION, PARAMETER :: YAXIS_THRESH = 0.001D0  ! Axisymmetric Y-axis threshold
    DOUBLE PRECISION, PARAMETER :: AMU_THRESH   = 0.001D0  ! Viscosity threshold

END MODULE CONSTANTS
