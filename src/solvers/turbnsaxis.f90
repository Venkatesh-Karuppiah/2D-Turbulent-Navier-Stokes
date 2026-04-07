!===============================================================================
! turbnsaxis.f90 - Axisymmetric Turbulent Navier-Stokes Solver
!===============================================================================
PROGRAM TURBNSAXIS
    USE CONSTANTS
    USE SOLVER_STATE
    USE GEOMETRY
    USE VISCOUS
    USE SOURCE_TERMS
    USE BOUNDARY
    USE RESIDUAL
    IMPLICIT NONE

    ! Set Axisymmetric mode
    DATA IAXIS / 1 /

    !---------------------------------------------------------------------------
    ! Read input file
    !---------------------------------------------------------------------------
    OPEN(25, FILE='flow.in')
    READ(25, *) FNAMEGRID, FNAMEBC, GAMA, AMOL, LOCAL, NITER, ISTART, CFL, &
                PRANDTL, PRANDTLT, P0, T0, PAMB, AMACHINIT, RTRED, NPSEUDO, IAXIS
    CLOSE(25)

    RUNIV = RUNIV_VAL
    GTH = GTH_VAL
    CONVKSC = CONVKSC_VAL

    !---------------------------------------------------------------------------
    ! Compute geometry
    !---------------------------------------------------------------------------
    CALL COMPUTE_GEOMETRY

    !---------------------------------------------------------------------------
    ! Thermodynamic properties
    !---------------------------------------------------------------------------
    RR = RUNIV / AMOL
    CV = RR / (GAMA - 1.0D0)
    CP = CV * GAMA
    RHO0 = P0 / (RR*T0)
    S0 = P0 / (RHO0**GAMA)
    TINF = T0 / (1.0D0 + 0.5D0*(GAMA-1.0D0)*AMINF*AMINF)

    VISU = 0.0D0; VISV = 0.0D0; VIST = 0.0D0; VISR = 0.0D0
    SOURCE5 = 0.0D0
    ITEST = 0
    CFLTURB = CFLTURB_VAL
    NITEROLD = 0; NOLD = 0

    !---------------------------------------------------------------------------
    ! Initialize solution
    !---------------------------------------------------------------------------
    IF (ISTART .EQ. 0) THEN
        DYNRATIO = 1.0D0 + (GAMA-1.0D0)*0.5D0*AMACHINIT**2.0D0
        TINIT = T0 / DYNRATIO
        PINIT = P0 / (DYNRATIO**(GAMA/(GAMA-1.0D0)))
        UINIT = AMACHINIT * DSQRT(GAMA*RR*TINIT)
        VINIT = 0.0D0
        EINIT = CV*TINIT + 0.5D0*(UINIT**2.0D0 + VINIT**2.0D0)
        RHOINIT = PINIT / (RR*TINIT)
        AMUINIT = (0.00000011848D0)*DSQRT(AMOL)*(TINIT**0.6D0)
        RINIT = AMUINIT / RHOINIT
        RTOLERANCE = RINIT * RTRED

        DO I = 1, NELES
            C1(I) = RHOINIT
            C2(I) = RHOINIT*UINIT
            C3(I) = RHOINIT*VINIT
            C4(I) = RHOINIT*EINIT
            C5(I) = RHOINIT*RINIT
        END DO

    ELSE IF (ISTART .EQ. 1) THEN
        NOLD = 1
        OPEN(55, FILE='restart.in')
        DO I = 1, NELES
            READ(55, *) C1(I), C2(I), C3(I), C4(I), C5(I)
        END DO
        READ(55, *) RTRED
        READ(55, *) NPSEUDO
        READ(55, *) RINIT
        READ(55, *) NITEROLD
        CLOSE(55)
        RTOLERANCE = RINIT * RTRED
    END IF

    DO I = 1, NELES
        COLD1(I)=C1(I); COLD2(I)=C2(I); COLD3(I)=C3(I)
        COLD4(I)=C4(I); COLD5(I)=C5(I)
    END DO

    !---------------------------------------------------------------------------
    ! RK4 time integration
    !---------------------------------------------------------------------------
    DO ITER = 1, NITER
        !--- Compute time step ---
        DTMIN = 100000.0D0
        DO I = 1, NELES
            RHO = C1(I)
            U(I) = C2(I)/RHO
            V(I) = C3(I)/RHO
            VEL = DSQRT(U(I)*U(I) + V(I)*V(I))
            E = C4(I)/RHO
            T(I) = (E - 0.5D0*(U(I)*U(I)+V(I)*V(I)))/CV
            P(I) = RHO*RR*T(I)
            DT(I) = CFL*DL(I) / (DABS(VEL) + DSQRT(GAMA*RR*T(I)))
            DTMIN = DMIN1(DTMIN, DT(I))
        END DO

        IF (LOCAL .EQ. 0) THEN
            DO I = 1, NELES
                DT(I) = DTMIN
            END DO
        END IF

        CALL DISSIPATION

        !--- Stage 1 ---
        CALL COMPUTE_PRIMITIVE_VARS
        CALL APPLY_BOUNDARY
        CALL COMPUTE_VISCOUS
        CALL COMPUTE_SOURCE_TERMS
        CALL COMPUTE_RESIDUAL
        CALL UPDATE_SOLUTION(ALPHA1_VAL)

        !--- Stage 2 ---
        CALL COMPUTE_PRIMITIVE_VARS
        CALL APPLY_BOUNDARY
        CALL COMPUTE_SOURCE_TERMS
        CALL COMPUTE_RESIDUAL
        CALL UPDATE_SOLUTION(ALPHA2_VAL)

        !--- Stage 3 ---
        CALL COMPUTE_PRIMITIVE_VARS
        CALL APPLY_BOUNDARY
        CALL COMPUTE_SOURCE_TERMS
        CALL COMPUTE_RESIDUAL
        CALL UPDATE_SOLUTION(ALPHA3_VAL)

        !--- Stage 4 ---
        CALL COMPUTE_PRIMITIVE_VARS
        CALL APPLY_BOUNDARY
        CALL COMPUTE_SOURCE_TERMS
        CALL COMPUTE_RESIDUAL
        CALL UPDATE_SOLUTION(ALPHA4_VAL)

        !--- Convergence monitoring ---
        CDIFF1 = CDIFF1_INIT; CDIFF2 = CDIFF1_INIT; CDIFF3 = CDIFF1_INIT
        CDIFF4 = CDIFF1_INIT; CDIFF5 = CDIFF5_INIT
        IF (MOD(ITER,1000) .EQ. 0) THEN
            DO I = 1, NELES
                CDIFF1 = DMAX1(DABS(COLD1(I)-C1(I)), CDIFF1)
                CDIFF2 = DMAX1(DABS(COLD2(I)-C2(I)), CDIFF2)
                CDIFF3 = DMAX1(DABS(COLD3(I)-C3(I)), CDIFF3)
                CDIFF4 = DMAX1(DABS(COLD4(I)-C4(I)), CDIFF4)
                CDIFF5 = DMAX1(DABS(COLD5(I)-C5(I)), CDIFF5)
            END DO
            WRITE(32, 11) CDIFF1, CDIFF2, CDIFF3, CDIFF4, CDIFF5
11          FORMAT(5E15.5)
            PRINT 10, ITER+NITEROLD, ' ', CDIFF1, ' ', CDIFF2, ' ', CDIFF3, &
                   ' ', CDIFF4, ' ', CDIFF5
10          FORMAT(I10, A5, E12.6, A5, E12.6, A5, E12.6, A5, E12.6, A5, E12.6)
        END IF

        DO I = 1, NELES
            COLD1(I)=C1(I); COLD2(I)=C2(I); COLD3(I)=C3(I)
            COLD4(I)=C4(I); COLD5(I)=C5(I)
        END DO

        !--- Periodic restart file ---
        IF (MOD(NITEROLD+ITER, 20000) .EQ. 0) THEN
            WRITE(RESTARTTEMPFILE, "(A7, I10, A3)") "restart", NITEROLD+ITER, ".in"
            OPEN(31, FILE=RESTARTTEMPFILE)
            DO I = 1, NELES
                WRITE(31, *) C1(I), C2(I), C3(I), C4(I), C5(I)
            END DO
            WRITE(31, *) RTRED
            WRITE(31, *) NPSEUDO
            WRITE(31, *) RINIT
            WRITE(31, *) ITER+NITEROLD
            CLOSE(31)
        END IF
    END DO
    ! END OF JAMESON'S RK4

    NITEROLD = NITEROLD + NITER
    OPEN(31, FILE='restart.in')
    DO I = 1, NELES
        WRITE(31, *) C1(I), C2(I), C3(I), C4(I), C5(I)
    END DO
    WRITE(31, *) RTRED
    WRITE(31, *) NPSEUDO
    WRITE(31, *) RINIT
    WRITE(31, *) NITEROLD
    CLOSE(31)

    CALL POSTPROCESS

END PROGRAM TURBNSAXIS

!===============================================================================
! Helper: compute primitive variables from conserved variables
!===============================================================================
SUBROUTINE COMPUTE_PRIMITIVE_VARS
    USE SOLVER_STATE, ONLY: NELES, C1, C2, C3, C4, C5, U, V, T, P, F1, F2, &
        F3, F4, F5, G1, G2, G3, G4, G5, RHO, E, RR, CV, GAMA, RT, I
    IMPLICIT NONE
    DO I = 1, NELES
        RHO = C1(I)
        U(I) = C2(I)/RHO
        V(I) = C3(I)/RHO
        E = C4(I)/RHO
        T(I) = (E - 0.5D0*(U(I)*U(I)+V(I)*V(I)))/CV
        P(I) = RHO*RR*T(I)
        F1(I) = C2(I)
        F2(I) = RHO*U(I)*U(I) + P(I)
        F3(I) = RHO*U(I)*V(I)
        F4(I) = (C4(I)+P(I))*U(I)
        G1(I) = C3(I)
        G2(I) = F3(I)
        G3(I) = RHO*V(I)*V(I) + P(I)
        G4(I) = (C4(I)+P(I))*V(I)
        RT = C5(I)/RHO
        F5(I) = RHO*RT*U(I)
        G5(I) = RHO*RT*V(I)
    END DO
END SUBROUTINE COMPUTE_PRIMITIVE_VARS
