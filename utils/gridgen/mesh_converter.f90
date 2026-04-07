!===============================================================================
! mesh_converter.f90 - GAMBIT Neutral to Solver Mesh Format Converter
!
! Reads:    tape7.neu  (GAMBIT neutral file format)
! Writes:   grid.in    (Solver mesh: nodes + elements with neighbor info)
!           bc.in      (Solver boundary: ghost cell definitions)
!
! Auto-detects boundary types based on coordinate positions:
!   Type 2  →  Inlet/outlet (at domain extents)
!   Type 3  →  Farfield (at xmax)
!   Type 4  →  Symmetry (at y=0)
!   Type 5  →  Wall (default for exposed faces)
!===============================================================================
PROGRAM MESH_CONVERT
    IMPLICIT NONE

    INTEGER, PARAMETER :: NN = 90000
    INTEGER, PARAMETER :: TOL_MATCH = 1.0D-4   ! Coordinate tolerance for BC detection

    ! Mesh data
    INTEGER :: NNODES, NELES, N, I, J, NEL
    INTEGER :: N1, N2, N3, N4, M1, M2, M3, M4
    INTEGER :: ITTEST
    DOUBLE PRECISION :: X(NN), Y(NN)
    INTEGER :: LNODE(NN, 4), NC1(NN), NC2(NN), NC3(NN), NC4(NN)

    ! Boundary detection
    INTEGER :: NGHOST, NTYPE, NSIDE
    DOUBLE PRECISION :: XMIN, XMAX, YMIN, YMAX, Z

    ! I/O unit numbers
    INTEGER, PARAMETER :: IUNIT_NEU  = 21
    INTEGER, PARAMETER :: IUNIT_GRID = 22
    INTEGER, PARAMETER :: IUNIT_BCTMP= 23
    INTEGER, PARAMETER :: IUNIT_BC   = 25
    INTEGER, PARAMETER :: IUNIT_DBG  = 27

    !---------------------------------------------------------------------------
    ! Open files
    !---------------------------------------------------------------------------
    OPEN(IUNIT_NEU,   FILE='tape7.neu',    STATUS='OLD')
    OPEN(IUNIT_GRID,  FILE='grid.dat',     STATUS='REPLACE')
    OPEN(IUNIT_BCTMP, FILE='bctemp.dat',   STATUS='REPLACE')
    OPEN(IUNIT_DBG,   FILE='grid_debug.dat', STATUS='REPLACE')

    !---------------------------------------------------------------------------
    ! Read header
    !---------------------------------------------------------------------------
    READ(IUNIT_NEU, *) NNODES, NELES
    WRITE(IUNIT_GRID, *) NNODES, NELES

    !---------------------------------------------------------------------------
    ! Read nodes and compute domain extents
    !---------------------------------------------------------------------------
    XMAX = -1.0D10; XMIN =  1.0D10
    YMAX = -1.0D10; YMIN =  1.0D10
    Z = 0.0D0

    DO I = 1, NNODES
        READ(IUNIT_NEU, *) N, X(N), Y(N)
        WRITE(IUNIT_GRID, 100) N, X(N), Y(N), Z
        WRITE(IUNIT_DBG,  100) N, X(N), Y(N), Z
        XMAX = MAX(XMAX, X(N)); XMIN = MIN(XMIN, X(N))
        YMAX = MAX(YMAX, Y(N)); YMIN = MIN(YMIN, Y(N))
    END DO

100 FORMAT(I8, 3F18.10)

    PRINT *, 'Domain extents: X = [', XMIN, ', ', XMAX, ']'
    PRINT *, '                Y = [', YMIN, ', ', YMAX, ']'
    PRINT *, 'Nodes =', NNODES, '  Elements =', NELES

    !---------------------------------------------------------------------------
    ! Read element connectivity
    !---------------------------------------------------------------------------
    DO I = 1, NELES
        READ(IUNIT_NEU, *) NEL, N1, N2, N3, N4
        LNODE(NEL, 1) = N1
        LNODE(NEL, 2) = N2
        LNODE(NEL, 3) = N3
        LNODE(NEL, 4) = N4
    END DO

    !---------------------------------------------------------------------------
    ! Build neighbor connectivity: O(N²) brute-force face matching
    !---------------------------------------------------------------------------
    NGHOST = NELES
    NC1 = 0; NC2 = 0; NC3 = 0; NC4 = 0

    DO I = 1, NELES
        N1 = LNODE(I, 1); N2 = LNODE(I, 2)
        N3 = LNODE(I, 3); N4 = LNODE(I, 4)
        ITTEST = 0

        DO J = 1, NELES
            IF (I == J) CYCLE
            M1 = LNODE(J, 1); M2 = LNODE(J, 2)
            M3 = LNODE(J, 3); M4 = LNODE(J, 4)

            IF (FACES_MATCH(N1, N2, M1, M2, M3, M4)) THEN
                NC1(I) = J; ITTEST = ITTEST + 1
            ELSE IF (FACES_MATCH(N2, N3, M1, M2, M3, M4)) THEN
                NC2(I) = J; ITTEST = ITTEST + 1
            ELSE IF (FACES_MATCH(N3, N4, M1, M2, M3, M4)) THEN
                NC3(I) = J; ITTEST = ITTEST + 1
            ELSE IF (FACES_MATCH(N4, N1, M1, M2, M3, M4)) THEN
                NC4(I) = J; ITTEST = ITTEST + 1
            END IF

            IF (ITTEST >= 4) EXIT
        END DO

        ! Identify boundary faces and auto-detect type
        CALL CLASSIFY_BOUNDARY(I, N1, N2, 1, NC1(I), NGHOST, X, Y, XMIN, XMAX, YMIN, YMAX)
        CALL CLASSIFY_BOUNDARY(I, N2, N3, 2, NC2(I), NGHOST, X, Y, XMIN, XMAX, YMIN, YMAX)
        CALL CLASSIFY_BOUNDARY(I, N3, N4, 3, NC3(I), NGHOST, X, Y, XMIN, XMAX, YMIN, YMAX)
        CALL CLASSIFY_BOUNDARY(I, N4, N1, 4, NC4(I), NGHOST, X, Y, XMIN, XMAX, YMIN, YMAX)

        ! Write element line to grid output
        WRITE(IUNIT_GRID, 200) I, N1, N2, N3, N4, NC1(I), NC2(I), NC3(I), NC4(I)
    END DO

200 FORMAT(9I8)

    !---------------------------------------------------------------------------
    ! Write boundary condition file
    !---------------------------------------------------------------------------
    CLOSE(IUNIT_NEU)
    CLOSE(IUNIT_GRID)
    CLOSE(IUNIT_BCTMP)

    OPEN(IUNIT_BCTMP, FILE='bctemp.dat', STATUS='OLD')
    OPEN(IUNIT_BC,    FILE='bc.in',      STATUS='REPLACE')

    WRITE(IUNIT_BC, *) NGHOST - NELES
    PRINT *, 'Total boundary faces:', NGHOST - NELES

    DO I = 1, NGHOST - NELES
        READ(IUNIT_BCTMP, *) N, NGHOST, NTYPE, NSIDE
        WRITE(IUNIT_BC, *) N, NGHOST, NTYPE, NSIDE
    END DO

    CLOSE(IUNIT_BCTMP)
    CLOSE(IUNIT_BC)
    CLOSE(IUNIT_DBG)

    ! Clean up temporary file
    CALL REMOVE_FILE('bctemp.dat')

    PRINT *, 'Mesh conversion complete. Output: grid.in, bc.in'

CONTAINS

    !===========================================================================
    ! Check if two faces share the same pair of nodes (in any order)
    !===========================================================================
    LOGICAL FUNCTION FACES_MATCH(NA, NB, MA, MB, MC, MD)
        INTEGER, INTENT(IN) :: NA, NB, MA, MB, MC, MD
        FACES_MATCH = &
            (NA == MB .AND. NB == MA) .OR. &
            (NA == MA .AND. NB == MB) .OR. &
            (NA == MC .AND. NB == MB) .OR. &
            (NA == MB .AND. NB == MC) .OR. &
            (NA == MD .AND. NB == MC) .OR. &
            (NA == MC .AND. NB == MD) .OR. &
            (NA == MA .AND. NB == MD) .OR. &
            (NA == MD .AND. NB == MA)
    END FUNCTION FACES_MATCH

    !===========================================================================
    ! Classify a face: if no neighbor, assign ghost ID and auto-detect BC type
    !===========================================================================
    SUBROUTINE CLASSIFY_BOUNDARY(IELEM, NA, NB, FACE, NC_VAL, NGCOUNT, XARR, YARR, &
                                 XLO, XHI, YLO, YHI)
        INTEGER, INTENT(IN)    :: IELEM, NA, NB, FACE
        INTEGER, INTENT(INOUT) :: NC_VAL, NGCOUNT
        DOUBLE PRECISION, INTENT(IN) :: XARR(NN), YARR(NN)
        DOUBLE PRECISION, INTENT(IN) :: XLO, XHI, YLO, YHI
        INTEGER :: NT

        IF (NC_VAL /= 0) RETURN   ! Has a neighbor → interior face

        NGCOUNT = NGCOUNT + 1
        NC_VAL = NGCOUNT
        NT = 5   ! Default: wall

        ! Auto-detect boundary type from coordinate position
        IF (ABS(XARR(NA) - XLO) < TOL_MATCH .AND. ABS(XARR(NB) - XLO) < TOL_MATCH) THEN
            NT = 2   ! Inlet at xmin
        END IF
        IF (ABS(XARR(NA) - XHI) < TOL_MATCH .AND. ABS(XARR(NB) - XHI) < TOL_MATCH) THEN
            NT = 3   ! Farfield at xmax
        END IF
        IF (ABS(YARR(NA) - YHI) < TOL_MATCH .AND. ABS(YARR(NB) - YHI) < TOL_MATCH) THEN
            NT = 2   ! Outlet at ymax
        END IF
        IF (ABS(YARR(NA)) < TOL_MATCH .AND. ABS(YARR(NB)) < TOL_MATCH) THEN
            NT = 4   ! Symmetry at y=0
        END IF

        WRITE(IUNIT_BCTMP, *) IELEM, NGCOUNT, NT, FACE
    END SUBROUTINE CLASSIFY_BOUNDARY

    !===========================================================================
    ! Portable file removal
    !===========================================================================
    SUBROUTINE REMOVE_FILE(FNAME)
        CHARACTER(*), INTENT(IN) :: FNAME
        INTEGER :: IERR
        INQUIRE(FILE=FNAME, EXIST=IERR)
        IF (IERR) THEN
            OPEN(NEWUNIT=IERR, FILE=FNAME, STATUS='SCRATCH')
            CLOSE(IERR, STATUS='DELETE')
        END IF
    END SUBROUTINE REMOVE_FILE

END PROGRAM MESH_CONVERT
