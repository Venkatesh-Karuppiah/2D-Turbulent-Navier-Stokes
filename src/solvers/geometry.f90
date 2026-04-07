!===============================================================================
! geometry.f90 - Mesh reading and geometric property computation
!                Supports both 2D Cartesian (IAXIS=0) and Axisymmetric (IAXIS=1)
!===============================================================================
MODULE GEOMETRY
    USE SOLVER_STATE, ONLY: NN, NODS, NELES, NGHOSTS, FNAMEGRID, FNAMEBC,     &
        NOD, NC1, NC2, NC3, NC4, X, Y, XCC, YCC, VOL, DL, SC1X, SC2X, SC3X, &
        SC4X, SC1Y, SC2Y, SC3Y, SC4Y, NPARENT, NGHOST, NTYPE, NSIDE, ITEST,  &
        N1, N2, N3, N4, NG, NP, NS, X1, X2, X3, X4, Y1, Y2, Y3, Y4,         &
        DL12, DL23, DL34, DL41, XCCFACE1, XCCFACE2, XCCFACE3, XCCFACE4,     &
        YCCFACE1, YCCFACE2, YCCFACE3, YCCFACE4, I
    USE CONSTANTS,  ONLY: MM_TO_M, HALF
    IMPLICIT NONE

    INTEGER :: IAXIS  ! 0 = 2D Cartesian, 1 = Axisymmetric

CONTAINS

    !===========================================================================
    ! GEOMETRY - Read mesh and compute geometric properties
    !===========================================================================
    SUBROUTINE COMPUTE_GEOMETRY
        IMPLICIT NONE

        NC1 = -1
        NC2 = -1
        NC3 = -1
        NC4 = -1

        ! Read grid file
        OPEN(21, FILE=FNAMEGRID)
        READ(21, *) NODS, NELES
        DO I = 1, NODS
            READ(21, *) N, X(N), Y(N)
            X(N) = X(N) * MM_TO_M
            Y(N) = Y(N) * MM_TO_M
        END DO

        ITEST = 0
        DO I = 1, NELES
            READ(21, *) N, NOD(I,1), NOD(I,2), NOD(I,3), NOD(I,4), &
                        NC1(I), NC2(I), NC3(I), NC4(I)
        END DO
        CLOSE(21)

        ! Read boundary condition file
        OPEN(22, FILE=FNAMEBC)
        READ(22, *) NGHOSTS
        DO I = 1, NGHOSTS
            READ(22, *) NPARENT(I), NGHOST(I), NTYPE(I), NSIDE(I)
        END DO
        CLOSE(22)

        ! Compute element geometric properties
        DO I = 1, NELES
            N1 = NOD(I,1)
            N2 = NOD(I,2)
            N3 = NOD(I,3)
            N4 = NOD(I,4)
            ITEST(N1) = ITEST(N1) + 1
            ITEST(N2) = ITEST(N2) + 1
            ITEST(N3) = ITEST(N3) + 1
            ITEST(N4) = ITEST(N4) + 1
            X1 = X(N1); X2 = X(N2); X3 = X(N3); X4 = X(N4)
            Y1 = Y(N1); Y2 = Y(N2); Y3 = Y(N3); Y4 = Y(N4)

            XCC(I) = (X1 + X2 + X3 + X4) * 0.25D0
            YCC(I) = (Y1 + Y2 + Y3 + Y4) * 0.25D0

            DL12 = DSQRT((X1-X2)*(X1-X2) + (Y1-Y2)*(Y1-Y2))
            DL23 = DSQRT((X2-X3)*(X2-X3) + (Y2-Y3)*(Y2-Y3))
            DL34 = DSQRT((X3-X4)*(X3-X4) + (Y3-Y4)*(Y3-Y4))
            DL41 = DSQRT((X4-X1)*(X4-X1) + (Y4-Y1)*(Y4-Y1))
            DL(I) = DMIN1(DL12, DL23, DL34, DL41)

            ! Surface vectors and volume differ between 2D and axisymmetric
            IF (IAXIS .EQ. 1) THEN
                ! Axisymmetric: include Y-weighting
                SC1X(I) = (Y2-Y1) * (Y1+Y2) * HALF
                SC1Y(I) = -(X2-X1) * (Y1+Y2) * HALF
                SC2X(I) = (Y3-Y2) * (Y2+Y3) * HALF
                SC2Y(I) = -(X3-X2) * (Y2+Y3) * HALF
                SC3X(I) = (Y4-Y3) * (Y3+Y4) * HALF
                SC3Y(I) = -(X4-X3) * (Y3+Y4) * HALF
                SC4X(I) = (Y1-Y4) * (Y4+Y1) * HALF
                SC4Y(I) = -(X1-X4) * (Y4+Y1) * HALF
                VOL(I) = HALF * ((X3-X1)*(Y4-Y2) - (X4-X2)*(Y3-Y1)) * YCC(I)
            ELSE
                ! 2D Cartesian
                SC1X(I) =  (Y2 - Y1)
                SC1Y(I) = -(X2 - X1)
                SC2X(I) =  (Y3 - Y2)
                SC2Y(I) = -(X3 - X2)
                SC3X(I) =  (Y4 - Y3)
                SC3Y(I) = -(X4 - X3)
                SC4X(I) =  (Y1 - Y4)
                SC4Y(I) = -(X1 - X4)
                VOL(I) = HALF * ((X3-X1)*(Y4-Y2) - (X4-X2)*(Y3-Y1))
            END IF

            IF (VOL(I) .LE. 0.0D0) VOL(I) = -VOL(I)
        END DO

        ! Compute ghost cell properties
        DO I = 1, NGHOSTS
            NP = NPARENT(I)
            NG = NGHOST(I)
            NS = NSIDE(I)
            N1 = NOD(NP,1); N2 = NOD(NP,2)
            N3 = NOD(NP,3); N4 = NOD(NP,4)
            X1 = X(N1); X2 = X(N2); X3 = X(N3); X4 = X(N4)
            Y1 = Y(N1); Y2 = Y(N2); Y3 = Y(N3); Y4 = Y(N4)

            IF (NS .EQ. 1) THEN
                XCCFACE1 = HALF*(X1 + X2)
                YCCFACE1 = HALF*(Y1 + Y2)
                XCC(NG) = 2.0D0*XCCFACE1 - XCC(NP)
                YCC(NG) = 2.0D0*YCCFACE1 - YCC(NP)
                ITEST(N1) = ITEST(N1) + 1
                ITEST(N2) = ITEST(N2) + 1
            ELSE IF (NS .EQ. 2) THEN
                XCCFACE2 = HALF*(X2 + X3)
                YCCFACE2 = HALF*(Y2 + Y3)
                XCC(NG) = 2.0D0*XCCFACE2 - XCC(NP)
                YCC(NG) = 2.0D0*YCCFACE2 - YCC(NP)
                ITEST(N2) = ITEST(N2) + 1
                ITEST(N3) = ITEST(N3) + 1
            ELSE IF (NS .EQ. 3) THEN
                XCCFACE3 = HALF*(X3 + X4)
                YCCFACE3 = HALF*(Y3 + Y4)
                XCC(NG) = 2.0D0*XCCFACE3 - XCC(NP)
                YCC(NG) = 2.0D0*YCCFACE3 - YCC(NP)
                ITEST(N3) = ITEST(N3) + 1
                ITEST(N4) = ITEST(N4) + 1
            ELSE IF (NS .EQ. 4) THEN
                XCCFACE4 = HALF*(X4 + X1)
                YCCFACE4 = HALF*(Y4 + Y1)
                XCC(NG) = 2.0D0*XCCFACE4 - XCC(NP)
                YCC(NG) = 2.0D0*YCCFACE4 - YCC(NP)
                ITEST(N1) = ITEST(N1) + 1
                ITEST(N4) = ITEST(N4) + 1
            END IF

            VOL(NG) = VOL(NP)
            NC1(NG) = NG
            NC2(NG) = NG
            NC3(NG) = NG
            NC4(NG) = NG
        END DO

    END SUBROUTINE COMPUTE_GEOMETRY

END MODULE GEOMETRY
