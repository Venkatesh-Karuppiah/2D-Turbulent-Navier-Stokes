!===============================================================================
! solver_state.f90 - Shared state variables and common subroutines
!                    Used by both 2D Cartesian and Axisymmetric solvers
!===============================================================================
MODULE SOLVER_STATE
    USE CONSTANTS, ONLY: NN, C02_VAL, C04_VAL, HALF
    IMPLICIT NONE

    !---------------------------------------------------------------------------
    ! File names and control strings
    !---------------------------------------------------------------------------
    CHARACTER*50 :: FNAMEGRID, FNAMEBC, CHOICE, RESTARTTEMPFILE

    !---------------------------------------------------------------------------
    ! Integer scalars
    !---------------------------------------------------------------------------
    INTEGER :: NODS, NELES, NGHOSTS, LOCAL, NITER, ISTART, IAXIS, IVISCOUS, &
        NCC1, NCC2, NCC3, NCC4, N1, N2, N3, N4, NG, NP, NT, NS, NCOPP,    &
        I, IE, ITER, K, J, L, N, NITEROLD, NOLD, nel, NPSEUDO, num, na, nbb

    !---------------------------------------------------------------------------
    ! Double precision scalars
    !---------------------------------------------------------------------------
    DOUBLE PRECISION :: ALPHA1, ALPHA2, ALPHA3, ALPHA4, ALPHA, BETA, GAMA,   &
        AMOL, PRANDTL, PRANDTLT, CFL, P0, T0, PAMB, AMACHINIT, CP, CV, RR,   &
        S0, DIV23, C02, C04, RTOLERANCE, ANR, DNR, DUX, DUY, DVX, DVY,      &
        DIVVEL, CMU, TXX, TXY, TYY, QXX, QYY, DQX, DQY, DRX, DRY, AMU, AK,  &
        AMUINIT, RHOINIT, VINIT, UINIT, TINIT, EINIT, RINIT, AMUT, AKT,     &
        AMUEFF, AKEFF, SX, SY, XA, XB, XC, XD, YA, YB, YC, YD, UA, UB, UC,  &
        UD, VA, VB, VC, VD, QA, QB, QC, QD, RA, RB, RC, RD, CHI, TTT, RHO0, &
        DYNRATIO, RHO, E, DTMIN, VEL, DTBYVOL, PHI, E2NC1, E2NC2, E2NC3,    &
        E2NC4, E4NC1, E4NC2, E4NC3, E4NC4, D2NC1, D2NC2, D2NC3, D2NC4,     &
        D24NC1, D24NC2, D24NC3, D24NC4, RI, FR2, FMU, PK, DEE, UAVG, VAVG,  &
        TAVG, DL12, DL23, DL34, DL41, XCCFACE1, XCCFACE2, XCCFACE3,         &
        XCCFACE4, YCCFACE1, YCCFACE2, YCCFACE3, YCCFACE4, RT, ENX, ENY,     &
        TFLAME, ANCOMB, ACOMB, CONVKSC, RHOPROP, D21, D22, D23, D2F, D25,   &
        D41, D42, D43, D4F, D45, ANUN, ANUD, YAVG, VBYY, CDIFF1, CDIFF2,    &
        CDIFF3, CDIFF4, CDIFF5, PNN, TNN, AMN, UNN, VNN, RON, DTX, DTY, TA, &
        TB, TC, TD, C5I, PIE, D4NC1, D4NC2, D4NC3, D4NC4, X1, X2, X3, X4,  &
        Y1, Y2, Y3, Y4, PINIT, THRUST, UI, VI, TI, RHOI, UNG, RNG, TNG,     &
        VNG, DFLOATITEST, RAVG, RHOAVG, RHOB, RHOD, RHONG, VOLBYDT, ANUI,   &
        ANURTI, DELTABL, PINF, TINF, AMINF, TWALL, F1I, F2I, F3I, F4I, F5I, &
        G1I, G2I, G3I, G4I, G5I, SX1, SX2, SX3, SX4, SY1, SY2, SY3, SY4,   &
        CFLTURB, RESNN, CSOUND, RHOA, RHOC, RTRED, BURNRATE, C1I, C2I, C3I, &
        C4I, D2Q1I, D2Q2I, D2Q3I, D2Q4I, D2Q5I, rdot, GTH, g, rdoteb,      &
        Gport, dhyd, DELY, RUNIV, CONST1, CONST2, CONST3

    !---------------------------------------------------------------------------
    ! Double precision arrays
    !---------------------------------------------------------------------------
    DOUBLE PRECISION, DIMENSION(NN) :: X, Y, XCC, YCC, VOL, DL, DT, SC1X,    &
        SC2X, SC3X, SC4X, SC1Y, SC2Y, SC3Y, SC4Y, C1, C2, C3, C4, C5, F1,    &
        F2, F3, F4, F5, G1, G2, G3, G4, G5, R1, R2, R3, R4, R5, COLD1,      &
        COLD2, COLD3, COLD4, COLD5, D2Q1, D2Q2, D2Q3, D2Q4, D2Q5, D241,     &
        D242, D243, D244, D245, U, V, P, T, VISU, VISV, VIST, VISR, UN, VN,  &
        TN, RN, RHON, SOURCE3, SOURCE5, ANU, ANURT

    !---------------------------------------------------------------------------
    ! Integer arrays
    !---------------------------------------------------------------------------
    INTEGER, DIMENSION(NN, 4) :: NOD
    INTEGER, DIMENSION(NN) :: NC1, NC2, NC3, NC4, ITEST, NPARENT, NGHOST,    &
        NTYPE, NSIDE

CONTAINS

    !===========================================================================
    ! DISSIPATION - Jameson's artificial dissipation scheme (identical in both)
    !===========================================================================
    SUBROUTINE DISSIPATION
        IMPLICIT NONE
        DOUBLE PRECISION :: E2N, E4N, D24N, D2N, D4N

        C02 = C02_VAL
        C04 = C04_VAL

        ! Copy ghost cells from parent
        DO I = 1, NGHOSTS
            NG = NGHOST(I)
            NP = NPARENT(I)
            P(NG) = P(NP)
            DT(NG) = DT(NP)
            C1(NG) = C1(NP)
            C2(NG) = C2(NP)
            C3(NG) = C3(NP)
            C4(NG) = C4(NP)
            C5(NG) = C5(NP)
        END DO

        ! Compute artificial viscosity coefficients
        DO I = 1, NELES
            NCC1 = NC1(I); NCC2 = NC2(I); NCC3 = NC3(I); NCC4 = NC4(I)
            PIE = P(I); C1I = C1(I); C5I = C5(I)
            ANUN = DABS(P(NCC1)-PIE)+DABS(P(NCC2)-PIE)+&
                   DABS(P(NCC3)-PIE)+DABS(P(NCC4)-PIE)
            ANUD = DABS(P(NCC1)+PIE)+DABS(P(NCC2)+PIE)+&
                   DABS(P(NCC3)+PIE)+DABS(P(NCC4)+PIE)
            ANU(I) = ANUN / ANUD

            ANUN = DABS(C5(NCC1)-C5I)+DABS(C5(NCC2)-C5I)+&
                   DABS(C5(NCC3)-C5I)+DABS(C5(NCC4)-C5I)
            ANUD = DABS(C5(NCC1)+C5I)+DABS(C5(NCC2)+C5I)+&
                   DABS(C5(NCC3)+C5I)+DABS(C5(NCC4)+C5I)
            ANURT(I) = ANUN / DMAX1(ANUD, 8.0D0*C1I*RTOLERANCE)

            D2Q1(I)=C1(NCC1)+C1(NCC2)+C1(NCC3)+C1(NCC4)-4.0D0*C1(I)
            D2Q2(I)=C2(NCC1)+C2(NCC2)+C2(NCC3)+C2(NCC4)-4.0D0*C2(I)
            D2Q3(I)=C3(NCC1)+C3(NCC2)+C3(NCC3)+C3(NCC4)-4.0D0*C3(I)
            D2Q4(I)=C4(NCC1)+C4(NCC2)+C4(NCC3)+C4(NCC4)-4.0D0*C4(I)
            D2Q5(I)=C5(NCC1)+C5(NCC2)+C5(NCC3)+C5(NCC4)-4.0D0*C5(I)
        END DO

        ! Copy dissipation to ghost cells
        DO I = 1, NGHOSTS
            NG = NGHOST(I); NP = NPARENT(I)
            ANU(NG)=ANU(NP); D2Q1(NG)=D2Q1(NP); D2Q2(NG)=D2Q2(NP)
            D2Q3(NG)=D2Q3(NP); D2Q4(NG)=D2Q4(NP); D2Q5(NG)=D2Q5(NP)
        END DO

        ! Compute dissipation contributions
        DO I = 1, NELES
            VOLBYDT = VOL(I)/DT(I)
            ANUI=ANU(I); ANURTI=ANURT(I)
            C1I=C1(I); C2I=C2(I); C3I=C3(I); C4I=C4(I); C5I=C5(I)
            D2Q1I=D2Q1(I); D2Q2I=D2Q2(I); D2Q3I=D2Q3(I)
            D2Q4I=D2Q4(I); D2Q5I=D2Q5(I)
            NCC1=NC1(I); NCC2=NC2(I); NCC3=NC3(I); NCC4=NC4(I)

            ! --- Momentum/Energy: Equation 1 ---
            E2N=C02*DMAX1(ANUI,ANU(NCC1)); E4N=DMAX1(0.0D0,C04-E2N)
            D24N=HALF*(VOLBYDT+VOL(NCC1)/DT(NCC1))
            D2N=E2N*D24N; D4N=E4N*D24N
            D21=D2N*(C1(NCC1)-C1I); D41=D4N*(D2Q1(NCC1)-D2Q1I)
            E2N=C02*DMAX1(ANUI,ANU(NCC2)); E4N=DMAX1(0.0D0,C04-E2N)
            D24N=HALF*(VOLBYDT+VOL(NCC2)/DT(NCC2))
            D2N=E2N*D24N; D4N=E4N*D24N
            D21=D21+D2N*(C1(NCC2)-C1I); D41=D41+D4N*(D2Q1(NCC2)-D2Q1I)
            E2N=C02*DMAX1(ANUI,ANU(NCC3)); E4N=DMAX1(0.0D0,C04-E2N)
            D24N=HALF*(VOLBYDT+VOL(NCC3)/DT(NCC3))
            D2N=E2N*D24N; D4N=E4N*D24N
            D21=D21+D2N*(C1(NCC3)-C1I); D41=D41+D4N*(D2Q1(NCC3)-D2Q1I)
            E2N=C02*DMAX1(ANUI,ANU(NCC4)); E4N=DMAX1(0.0D0,C04-E2N)
            D24N=HALF*(VOLBYDT+VOL(NCC4)/DT(NCC4))
            D2N=E2N*D24N; D4N=E4N*D24N
            D21=D21+D2N*(C1(NCC4)-C1I); D41=D41+D4N*(D2Q1(NCC4)-D2Q1I)
            D241(I)=D21-D41

            ! --- Equation 2 ---
            E2N=C02*DMAX1(ANUI,ANU(NCC1)); E4N=DMAX1(0.0D0,C04-E2N)
            D24N=HALF*(VOLBYDT+VOL(NCC1)/DT(NCC1))
            D2N=E2N*D24N; D4N=E4N*D24N
            D22=D2N*(C2(NCC1)-C2I); D42=D4N*(D2Q2(NCC1)-D2Q2I)
            E2N=C02*DMAX1(ANUI,ANU(NCC2)); E4N=DMAX1(0.0D0,C04-E2N)
            D24N=HALF*(VOLBYDT+VOL(NCC2)/DT(NCC2))
            D2N=E2N*D24N; D4N=E4N*D24N
            D22=D22+D2N*(C2(NCC2)-C2I); D42=D42+D4N*(D2Q2(NCC2)-D2Q2I)
            E2N=C02*DMAX1(ANUI,ANU(NCC3)); E4N=DMAX1(0.0D0,C04-E2N)
            D24N=HALF*(VOLBYDT+VOL(NCC3)/DT(NCC3))
            D2N=E2N*D24N; D4N=E4N*D24N
            D22=D22+D2N*(C2(NCC3)-C2I); D42=D42+D4N*(D2Q2(NCC3)-D2Q2I)
            E2N=C02*DMAX1(ANUI,ANU(NCC4)); E4N=DMAX1(0.0D0,C04-E2N)
            D24N=HALF*(VOLBYDT+VOL(NCC4)/DT(NCC4))
            D2N=E2N*D24N; D4N=E4N*D24N
            D22=D22+D2N*(C2(NCC4)-C2I); D42=D42+D4N*(D2Q2(NCC4)-D2Q2I)
            D242(I)=D22-D42

            ! --- Equation 3 ---
            E2N=C02*DMAX1(ANUI,ANU(NCC1)); E4N=DMAX1(0.0D0,C04-E2N)
            D24N=HALF*(VOLBYDT+VOL(NCC1)/DT(NCC1))
            D2N=E2N*D24N; D4N=E4N*D24N
            D23=D2N*(C3(NCC1)-C3I); D43=D4N*(D2Q3(NCC1)-D2Q3I)
            E2N=C02*DMAX1(ANUI,ANU(NCC2)); E4N=DMAX1(0.0D0,C04-E2N)
            D24N=HALF*(VOLBYDT+VOL(NCC2)/DT(NCC2))
            D2N=E2N*D24N; D4N=E4N*D24N
            D23=D23+D2N*(C3(NCC2)-C3I); D43=D43+D4N*(D2Q3(NCC2)-D2Q3I)
            E2N=C02*DMAX1(ANUI,ANU(NCC3)); E4N=DMAX1(0.0D0,C04-E2N)
            D24N=HALF*(VOLBYDT+VOL(NCC3)/DT(NCC3))
            D2N=E2N*D24N; D4N=E4N*D24N
            D23=D23+D2N*(C3(NCC3)-C3I); D43=D43+D4N*(D2Q3(NCC3)-D2Q3I)
            E2N=C02*DMAX1(ANUI,ANU(NCC4)); E4N=DMAX1(0.0D0,C04-E2N)
            D24N=HALF*(VOLBYDT+VOL(NCC4)/DT(NCC4))
            D2N=E2N*D24N; D4N=E4N*D24N
            D23=D23+D2N*(C3(NCC4)-C3I); D43=D43+D4N*(D2Q3(NCC4)-D2Q3I)
            D243(I)=D23-D43

            ! --- Equation 4 ---
            E2N=C02*DMAX1(ANUI,ANU(NCC1)); E4N=DMAX1(0.0D0,C04-E2N)
            D24N=HALF*(VOLBYDT+VOL(NCC1)/DT(NCC1))
            D2N=E2N*D24N; D4N=E4N*D24N
            D2F=D2N*(C4(NCC1)-C4I); D4F=D4N*(D2Q4(NCC1)-D2Q4I)
            E2N=C02*DMAX1(ANUI,ANU(NCC2)); E4N=DMAX1(0.0D0,C04-E2N)
            D24N=HALF*(VOLBYDT+VOL(NCC2)/DT(NCC2))
            D2N=E2N*D24N; D4N=E4N*D24N
            D2F=D2F+D2N*(C4(NCC2)-C4I); D4F=D4F+D4N*(D2Q4(NCC2)-D2Q4I)
            E2N=C02*DMAX1(ANUI,ANU(NCC3)); E4N=DMAX1(0.0D0,C04-E2N)
            D24N=HALF*(VOLBYDT+VOL(NCC3)/DT(NCC3))
            D2N=E2N*D24N; D4N=E4N*D24N
            D2F=D2F+D2N*(C4(NCC3)-C4I); D4F=D4F+D4N*(D2Q4(NCC3)-D2Q4I)
            E2N=C02*DMAX1(ANUI,ANU(NCC4)); E4N=DMAX1(0.0D0,C04-E2N)
            D24N=HALF*(VOLBYDT+VOL(NCC4)/DT(NCC4))
            D2N=E2N*D24N; D4N=E4N*D24N
            D2F=D2F+D2N*(C4(NCC4)-C4I); D4F=D4F+D4N*(D2Q4(NCC4)-D2Q4I)
            D244(I)=D2F-D4F

            ! --- Turbulence: Equation 5 ---
            E2N=C02*DMAX1(ANURTI,ANURT(NCC1)); E4N=DMAX1(0.0D0,C04-E2N)
            D24N=HALF*(VOLBYDT+VOL(NCC1)/DT(NCC1))
            D2N=E2N*D24N; D4N=E4N*D24N
            D25=D2N*(C5(NCC1)-C5I); D45=D4N*(D2Q5(NCC1)-D2Q5I)
            E2N=C02*DMAX1(ANURTI,ANURT(NCC2)); E4N=DMAX1(0.0D0,C04-E2N)
            D24N=HALF*(VOLBYDT+VOL(NCC2)/DT(NCC2))
            D2N=E2N*D24N; D4N=E4N*D24N
            D25=D25+D2N*(C5(NCC2)-C5I); D45=D45+D4N*(D2Q5(NCC2)-D2Q5I)
            E2N=C02*DMAX1(ANURTI,ANURT(NCC3)); E4N=DMAX1(0.0D0,C04-E2N)
            D24N=HALF*(VOLBYDT+VOL(NCC3)/DT(NCC3))
            D2N=E2N*D24N; D4N=E4N*D24N
            D25=D25+D2N*(C5(NCC3)-C5I); D45=D45+D4N*(D2Q5(NCC3)-D2Q5I)
            E2N=C02*DMAX1(ANURTI,ANURT(NCC4)); E4N=DMAX1(0.0D0,C04-E2N)
            D24N=HALF*(VOLBYDT+VOL(NCC4)/DT(NCC4))
            D2N=E2N*D24N; D4N=E4N*D24N
            D25=D25+D2N*(C5(NCC4)-C5I); D45=D45+D4N*(D2Q5(NCC4)-D2Q5I)
            D245(I)=D25-D45
        END DO

    END SUBROUTINE DISSIPATION

    !===========================================================================
    ! POSTPROCESS - Tecplot output generation (identical in both)
    !===========================================================================
    SUBROUTINE POSTPROCESS
        USE CONSTANTS, ONLY: PA_TO_BAR
        IMPLICIT NONE

        PRINT*, 'Postprocessing.'
        PRINT*, ' '

        C1=0.0D0; C2=0.0D0; C3=0.0D0; C4=0.0D0; C5=0.0D0
        ITEST=0

        DO I=1,NELES
            DO J=1,4
                K=NOD(I,J)
                ITEST(K)=ITEST(K)+1
                C1(K)=C1(K)+P(I); C2(K)=C2(K)+U(I)
                C3(K)=C3(K)+V(I); C4(K)=C4(K)+T(I)
            END DO
        END DO

        DO I=1,NGHOSTS
            NG=NGHOST(I); NP=NPARENT(I); NS=NSIDE(I)
            IF(NS.EQ.1)THEN; N1=NOD(NP,1); N2=NOD(NP,2)
            ELSE IF(NS.EQ.2)THEN; N1=NOD(NP,2); N2=NOD(NP,3)
            ELSE IF(NS.EQ.3)THEN; N1=NOD(NP,3); N2=NOD(NP,4)
            ELSE IF(NS.EQ.4)THEN; N1=NOD(NP,4); N2=NOD(NP,1)
            END IF
            C1(N1)=C1(N1)+P(NG); C2(N1)=C2(N1)+U(NG)
            C3(N1)=C3(N1)+V(NG); C4(N1)=C4(N1)+T(NG)
            C1(N2)=C1(N2)+P(NG); C2(N2)=C2(N2)+U(NG)
            C3(N2)=C3(N2)+V(NG); C4(N2)=C4(N2)+T(NG)
            ITEST(N1)=ITEST(N1)+1; ITEST(N2)=ITEST(N2)+1
        END DO

        DO I=1,NODS
            DFLOATITEST=DFLOAT(ITEST(I))
            C1(I)=C1(I)/DFLOATITEST; C2(I)=C2(I)/DFLOATITEST
            C3(I)=C3(I)/DFLOATITEST; C4(I)=C4(I)/DFLOATITEST
        END DO

        OPEN(UNIT=26,FILE='tec_plot.dat')
        WRITE(26,*)'VARIABLES = "X", "Y", "U", "V", "P", "RHO", "MACH", "T"'
        WRITE(26,*)'ZONE F = FEPOINT, ET = QUADRILATERAL, N = ',NODS,', E = ',NELES
        DO I=1,NODS
            PNN=C1(I); UNN=C2(I); VNN=C3(I); TNN=C4(I)
            AMN=DSQRT(UNN*UNN+VNN*VNN)/DSQRT(GAMA*RR*TNN)
            RON=PNN/(RR*TNN)
            WRITE(26,*)X(I),Y(I),UNN,VNN,PNN*PA_TO_BAR,RON,AMN,TNN
        END DO
        DO IE=1,NELES
            WRITE(26,*)NOD(IE,1),NOD(IE,2),NOD(IE,3),NOD(IE,4)
        END DO
        CLOSE(26)

    END SUBROUTINE POSTPROCESS

END MODULE SOLVER_STATE
