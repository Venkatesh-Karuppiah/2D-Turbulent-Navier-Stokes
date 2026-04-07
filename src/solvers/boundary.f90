!===============================================================================
! boundary.f90 - Ghost cell boundary condition computation
!                Supports both 2D Cartesian (IAXIS=0) and Axisymmetric (IAXIS=1)
!===============================================================================
MODULE BOUNDARY
    USE SOLVER_STATE, ONLY: NN, NGHOSTS, NPARENT, NGHOST, NTYPE, NSIDE, NOD,  &
        NC1, NC2, NC3, NC4, YCC, SC1X, SC1Y, SC2X, SC2Y, SC3X, SC3Y, SC4X,   &
        SC4Y, U, V, P, T, C1, C2, C3, C4, C5, F1, F2, F3, F4, F5, G1, G2,   &
        G3, G4, G5, CP, CV, RR, GAMA, S0, PAMB, AMOL, T0, P0, AMINF, TINF,   &
        PINF, DELTABL, TFLAME, ANCOMB, ACOMB, CONVKSC, RHOPROP, GTH,         &
        NP, NG, NT, NS, NCOPP, N1, N2, N3, N4, DHYD, DNR, ENX, ENY, DELY,    &
        Gport, RDOT, RDOTEB, VEL, AMU, RT, RHO, E, RTOLERANCE, RTRED, Y,     &
        AMUINIT, RINIT, N, I
    USE CONSTANTS, ONLY: PA_TO_BAR, BL_EXPONENT, GTH_VAL
    IMPLICIT NONE

    INTEGER :: IAXIS  ! 0 = 2D Cartesian, 1 = Axisymmetric

CONTAINS

    !===========================================================================
    ! APPLY_BOUNDARY - Set ghost cell values for all boundary types
    !===========================================================================
    SUBROUTINE APPLY_BOUNDARY
        IMPLICIT NONE

        DO I = 1, NGHOSTS
            NP = NPARENT(I)
            NG = NGHOST(I)
            NT = NTYPE(I)
            NS = NSIDE(I)

            IF (NT .EQ. 2) THEN
                CALL BC_SUBSONIC_INLET
            ELSE IF (NT .EQ. 3) THEN
                CALL BC_FARFIELD
            ELSE IF (NT .EQ. 4) THEN
                CALL BC_SYMMETRY
            ELSE IF (NT .EQ. 5) THEN
                CALL BC_WALL
            ELSE IF (NT .EQ. 6) THEN
                CALL BC_COMBUSTION
            ELSE IF (NT .EQ. 7) THEN
                CALL BC_BOUNDLAYER
            END IF
        END DO

    END SUBROUTINE APPLY_BOUNDARY

    !---------------------------------------------------------------------------
    ! BC type 2: Subsonic inlet (isentropic)
    !---------------------------------------------------------------------------
    SUBROUTINE BC_SUBSONIC_INLET
        IMPLICIT NONE
        U(NG) = DABS(U(NP))
        V(NG) = -V(NP)
        T(NG) = T0 - 0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))/CP
        RHO = (S0/(RR*T(NG)))**(1.0D0/(1.0D0-GAMA))
        P(NG) = RHO*RR*T(NG)
        AMU = (0.00000011848D0)*DSQRT(AMOL)*(T(NG)**0.6D0)
        RT = AMU/RHO
        C1(NG) = RHO
        C2(NG) = RHO*U(NG)
        C3(NG) = RHO*V(NG)
        E = CV*T(NG) + 0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))
        C4(NG) = RHO*E
        F1(NG)=C2(NG); F2(NG)=RHO*U(NG)*U(NG)+P(NG); F3(NG)=RHO*U(NG)*V(NG)
        F4(NG)=(C4(NG)+P(NG))*U(NG); G1(NG)=C3(NG); G2(NG)=F3(NG)
        G3(NG)=RHO*V(NG)*V(NG)+P(NG); G4(NG)=(C4(NG)+P(NG))*V(NG)
        C5(NG)=RHO*RT; F5(NG)=C5(NG)*U(NG); G5(NG)=C5(NG)*V(NG)
    END SUBROUTINE BC_SUBSONIC_INLET

    !---------------------------------------------------------------------------
    ! BC type 3: Farfield (Riemann extrapolation)
    !---------------------------------------------------------------------------
    SUBROUTINE BC_FARFIELD
        IMPLICIT NONE
        IF (NS .EQ. 1) THEN;      NCOPP=NC3(NP)
        ELSE IF (NS .EQ. 2) THEN; NCOPP=NC4(NP)
        ELSE IF (NS .EQ. 3) THEN; NCOPP=NC1(NP)
        ELSE IF (NS .EQ. 4) THEN; NCOPP=NC2(NP)
        END IF
        U(NG) = DABS(2.0D0*U(NP)-U(NCOPP))
        V(NG) = 2.0D0*V(NP)-V(NCOPP)
        T(NG) = 2.0D0*T(NP)-T(NCOPP)
        P(NG) = 2.0D0*P(NP)-P(NCOPP)
        IF (PAMB .GE. 0.0D0) P(NG) = PAMB
        RHO = P(NG)/(RR*T(NG))
        C1(NG)=RHO; C2(NG)=RHO*U(NG); C3(NG)=RHO*V(NG)
        E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG)); C4(NG)=RHO*E
        F1(NG)=C2(NG); F2(NG)=RHO*U(NG)*U(NG)+P(NG); F3(NG)=RHO*U(NG)*V(NG)
        F4(NG)=(C4(NG)+P(NG))*U(NG); G1(NG)=C3(NG); G2(NG)=F3(NG)
        G3(NG)=RHO*V(NG)*V(NG)+P(NG); G4(NG)=(C4(NG)+P(NG))*V(NG)
        C5(NG)=C5(NP); F5(NG)=F5(NP); G5(NG)=G5(NP)
    END SUBROUTINE BC_FARFIELD

    !---------------------------------------------------------------------------
    ! BC type 4: Symmetry / slip wall
    !---------------------------------------------------------------------------
    SUBROUTINE BC_SYMMETRY
        IMPLICIT NONE
        U(NG)=U(NP); V(NG)=-V(NP); T(NG)=T(NP); RHO=C1(NP); P(NG)=P(NP)
        C1(NG)=RHO; C2(NG)=RHO*U(NG); C3(NG)=RHO*V(NG)
        E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG)); C4(NG)=RHO*E
        F1(NG)=C2(NG); F2(NG)=RHO*U(NG)*U(NG)+P(NG); F3(NG)=RHO*U(NG)*V(NG)
        F4(NG)=(C4(NG)+P(NG))*U(NG); G1(NG)=C3(NG); G2(NG)=F3(NG)
        G3(NG)=RHO*V(NG)*V(NG)+P(NG); G4(NG)=(C4(NG)+P(NG))*V(NG)
        C5(NG)=C5(NP); F5(NG)=F5(NP); G5(NG)=-G5(NP)
    END SUBROUTINE BC_SYMMETRY

    !---------------------------------------------------------------------------
    ! BC type 5: No-slip wall (2D) or axis (axisymmetric)
    !---------------------------------------------------------------------------
    SUBROUTINE BC_WALL
        IMPLICIT NONE
        U(NG)=-U(NP); V(NG)=-V(NP); T(NG)=T(NP); RHO=C1(NP); P(NG)=P(NP)
        C1(NG)=RHO; C2(NG)=RHO*U(NG); C3(NG)=RHO*V(NG)
        E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG)); C4(NG)=RHO*E
        F1(NG)=-F1(NP); F2(NG)=2.0D0*P(NP)-F2(NP); F3(NG)=-F3(NP); F4(NG)=-F4(NP)
        G1(NG)=-G1(NP); G2(NG)=-G2(NP); G3(NG)=2.0D0*P(NP)-G3(NP); G4(NG)=-G4(NP)
        C5(NG)=-C5(NP); F5(NG)=-F5(NP); G5(NG)=-G5(NP)
    END SUBROUTINE BC_WALL

    !---------------------------------------------------------------------------
    ! BC type 6: Combustion / mass injection (with G-threshold correction)
    ! Differs between 2D and axisymmetric in GPORT computation
    !---------------------------------------------------------------------------
    SUBROUTINE BC_COMBUSTION
        IMPLICIT NONE
        DOUBLE PRECISION :: GVAL

        DHYD = DABS(YCC(NG) + YCC(NP))

        ! Compute mass flow rate GPORT along the boundary
        CALL COMPUTE_GPORT

        RHO = C1(NP)
        RDOT = ACOMB * ((P(NP)*CONVKSC)**(ANCOMB))
        RDOTEB = RDOT
        VEL = RHOPROP/RHO * RDOTEB
        P(NG) = P(NP)
        U(NG) = VEL*ENX + VEL*ENX - U(NP)
        V(NG) = VEL*ENY + VEL*ENY - V(NP)
        T(NG) = TFLAME - 0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))/CP
        AMU = (0.00000011848D0)*DSQRT(AMOL)*(T(NG)**0.6D0)
        GVAL = (GPORT/(RHOPROP*RDOT)) * ((RHOPROP*RDOT*DHYD/AMU)**(-0.125D0))
        IF (GVAL .GT. GTH) THEN
            RDOTEB = RDOTEB*(1.0D0 + 0.023D0*((GVAL**0.8D0)-(GTH**0.8D0)))
        END IF
        VEL = RHOPROP/RHO * RDOTEB
        U(NG) = VEL*ENX + VEL*ENX - U(NP)
        V(NG) = VEL*ENY + VEL*ENY - V(NP)
        T(NG) = TFLAME - 0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))/CP
        RHO = P(NG)/(RR*T(NG))
        AMU = (0.00000011848D0)*DSQRT(AMOL)*(T(NG)**0.6D0)
        RT = AMU/RHO
        C5(NG)=RHO*RT; C1(NG)=RHO; C2(NG)=RHO*U(NG); C3(NG)=RHO*V(NG)
        E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG)); C4(NG)=RHO*E
        F1(NG)=C2(NG); F2(NG)=RHO*U(NG)*U(NG)+P(NG); F3(NG)=RHO*U(NG)*V(NG)
        F4(NG)=(C4(NG)+P(NG))*U(NG); G1(NG)=C3(NG); G2(NG)=F3(NG)
        G3(NG)=RHO*V(NG)*V(NG)+P(NG); G4(NG)=(C4(NG)+P(NG))*V(NG)
        F5(NG)=C5(NG)*U(NG); G5(NG)=C5(NG)*V(NG)
    END SUBROUTINE BC_COMBUSTION

    !---------------------------------------------------------------------------
    ! GPORT computation: differs between 2D and axisymmetric
    !---------------------------------------------------------------------------
    SUBROUTINE COMPUTE_GPORT
        IMPLICIT NONE
        DOUBLE PRECISION :: SCX, SCY

        ! Get surface vector magnitude and normal
        IF (NS .EQ. 1) THEN; SCX=SC1X(NP); SCY=SC1Y(NP)
        ELSE IF (NS .EQ. 2) THEN; SCX=SC2X(NP); SCY=SC2Y(NP)
        ELSE IF (NS .EQ. 3) THEN; SCX=SC3X(NP); SCY=SC3Y(NP)
        ELSE IF (NS .EQ. 4) THEN; SCX=SC4X(NP); SCY=SC4Y(NP)
        END IF
        DNR = DSQRT(SCX*SCX + SCY*SCY)
        ENX = -SCX/DNR; ENY = -SCY/DNR

        ! Compute DELY and initial GPORT for parent element
        N1 = NOD(NP,1); N2 = NOD(NP,2); N3 = NOD(NP,3); N4 = NOD(NP,4)
        IF (NS .EQ. 1) THEN
            DELY = HALF*((Y(N1)+Y(N2))-(Y(N3)+Y(N4)))
            IF (IAXIS .EQ. 1) THEN
                GPORT = C1(NP)*U(NP)*DELY*YCC(NP)
            ELSE
                GPORT = C1(NP)*U(NP)*DELY
            END IF
            CALL ACCUMULATE_GPORT_ALONG_BOUNDARY(NC3, 3)
        ELSE IF (NS .EQ. 2) THEN
            DELY = HALF*((Y(N2)+Y(N3))-(Y(N1)+Y(N4)))
            IF (IAXIS .EQ. 1) THEN
                GPORT = C1(NP)*U(NP)*DELY*YCC(NP)
            ELSE
                GPORT = C1(NP)*U(NP)*DELY
            END IF
            CALL ACCUMULATE_GPORT_ALONG_BOUNDARY(NC4, 4)
        ELSE IF (NS .EQ. 3) THEN
            DELY = HALF*((Y(N3)+Y(N4))-(Y(N1)+Y(N2)))
            IF (IAXIS .EQ. 1) THEN
                GPORT = C1(NP)*U(NP)*DELY*YCC(NP)
            ELSE
                GPORT = C1(NP)*U(NP)*DELY
            END IF
            CALL ACCUMULATE_GPORT_ALONG_BOUNDARY(NC1, 1)
        ELSE IF (NS .EQ. 4) THEN
            DELY = HALF*((Y(N4)+Y(N1))-(Y(N2)+Y(N3)))
            IF (IAXIS .EQ. 1) THEN
                GPORT = C1(NP)*U(NP)*DELY*YCC(NP)
            ELSE
                GPORT = C1(NP)*U(NP)*DELY
            END IF
            CALL ACCUMULATE_GPORT_ALONG_BOUNDARY(NC2, 2)
        END IF

        ! Finalize GPORT (different normalization for 2D vs axisymmetric)
        IF (IAXIS .EQ. 1) THEN
            GPORT = 2.0D0 * GPORT / (DHYD**2)
        ELSE
            GPORT = GPORT / DHYD
        END IF
    END SUBROUTINE COMPUTE_GPORT

    !---------------------------------------------------------------------------
    ! Walk along the boundary accumulating GPORT
    ! NCARR: array accessor function (as integer index into NC1..NC4)
    ! SIDE: face number for DELY computation
    !---------------------------------------------------------------------------
    SUBROUTINE ACCUMULATE_GPORT_ALONG_BOUNDARY(NCARR, SIDE)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NCARR, SIDE
        INTEGER :: NCVAL

        N = NCARR
        N1 = NOD(N,1); N2 = NOD(N,2); N3 = NOD(N,3); N4 = NOD(N,4)

        IF (SIDE .EQ. 1) THEN
            DELY = HALF*((Y(N1)+Y(N2))-(Y(N3)+Y(N4)))
        ELSE IF (SIDE .EQ. 2) THEN
            DELY = HALF*((Y(N2)+Y(N3))-(Y(N1)+Y(N4)))
        ELSE IF (SIDE .EQ. 3) THEN
            DELY = HALF*((Y(N3)+Y(N4))-(Y(N1)+Y(N2)))
        ELSE IF (SIDE .EQ. 4) THEN
            DELY = HALF*((Y(N4)+Y(N1))-(Y(N2)+Y(N3)))
        END IF

        IF (IAXIS .EQ. 1) THEN
            GPORT = GPORT + C1(N)*U(N)*DELY*YCC(N)
        ELSE
            GPORT = GPORT + C1(N)*U(N)*DELY
        END IF
        NCVAL = NCARR
        N = NCVAL
        ! Walk along boundary until we loop back
        IF (SIDE .EQ. 1) THEN
            N = NC1(N)
            DO WHILE (NC1(N) .NE. N)
                N1 = NOD(N,1); N2 = NOD(N,2); N3 = NOD(N,3); N4 = NOD(N,4)
                DELY = HALF*((Y(N3)+Y(N4))-(Y(N1)+Y(N2)))
                IF (IAXIS .EQ. 1) THEN
                    GPORT = GPORT + C1(N)*U(N)*DELY*YCC(N)
                ELSE
                    GPORT = GPORT + C1(N)*U(N)*DELY
                END IF
                N = NC1(N)
            END DO
        ELSE IF (SIDE .EQ. 2) THEN
            N = NC2(N)
            DO WHILE (NC2(N) .NE. N)
                N1 = NOD(N,1); N2 = NOD(N,2); N3 = NOD(N,3); N4 = NOD(N,4)
                DELY = HALF*((Y(N4)+Y(N1))-(Y(N2)+Y(N3)))
                IF (IAXIS .EQ. 1) THEN
                    GPORT = GPORT + C1(N)*U(N)*DELY*YCC(N)
                ELSE
                    GPORT = GPORT + C1(N)*U(N)*DELY
                END IF
                N = NC2(N)
            END DO
        ELSE IF (SIDE .EQ. 3) THEN
            N = NC3(N)
            DO WHILE (NC3(N) .NE. N)
                N1 = NOD(N,1); N2 = NOD(N,2); N3 = NOD(N,3); N4 = NOD(N,4)
                DELY = HALF*((Y(N1)+Y(N2))-(Y(N3)+Y(N4)))
                IF (IAXIS .EQ. 1) THEN
                    GPORT = GPORT + C1(N)*U(N)*DELY*YCC(N)
                ELSE
                    GPORT = GPORT + C1(N)*U(N)*DELY
                END IF
                N = NC3(N)
            END DO
        ELSE IF (SIDE .EQ. 4) THEN
            N = NC4(N)
            DO WHILE (NC4(N) .NE. N)
                N1 = NOD(N,1); N2 = NOD(N,2); N3 = NOD(N,3); N4 = NOD(N,4)
                DELY = HALF*((Y(N2)+Y(N3))-(Y(N1)+Y(N4)))
                IF (IAXIS .EQ. 1) THEN
                    GPORT = GPORT + C1(N)*U(N)*DELY*YCC(N)
                ELSE
                    GPORT = GPORT + C1(N)*U(N)*DELY
                END IF
                N = NC4(N)
            END DO
        END IF
    END SUBROUTINE ACCUMULATE_GPORT_ALONG_BOUNDARY

    !---------------------------------------------------------------------------
    ! BC type 7: Boundary layer inlet profile
    !---------------------------------------------------------------------------
    SUBROUTINE BC_BOUNDLAYER
        IMPLICIT NONE
        P(NG) = PINF; T(NG) = TINF; V(NG) = -V(NP)
        RHO = P(NG)/(RR*T(NG))
        IF (YCC(NG) .LE. DELTABL) THEN
            U(NG) = 2.0D0*AMINF*DSQRT(GAMA*RR*TINF)*(YCC(NG)/DELTABL)**BL_EXPONENT - U(NP)
        ELSE
            U(NG) = 2.0D0*AMINF*DSQRT(GAMA*RR*TINF) - U(NP)
        END IF
        AMU = (0.00000011848D0)*DSQRT(AMOL)*(T(NG)**0.6D0)
        RT = AMU/RHO
        C1(NG)=RHO; C2(NG)=RHO*U(NG); C3(NG)=RHO*V(NG)
        E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG)); C4(NG)=RHO*E
        F1(NG)=C2(NG); F2(NG)=RHO*U(NG)*U(NG)+P(NG); F3(NG)=RHO*U(NG)*V(NG)
        F4(NG)=(C4(NG)+P(NG))*U(NG); G1(NG)=C3(NG); G2(NG)=F3(NG)
        G3(NG)=RHO*V(NG)*V(NG)+P(NG); G4(NG)=(C4(NG)+P(NG))*V(NG)
        C5(NG)=RHO*RT; F5(NG)=C5(NG)*U(NG); G5(NG)=C5(NG)*V(NG)
    END SUBROUTINE BC_BOUNDLAYER

END MODULE BOUNDARY
