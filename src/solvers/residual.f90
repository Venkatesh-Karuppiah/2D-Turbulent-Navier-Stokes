!===============================================================================
! residual.f90 - Inviscid flux residual and solution update computation
!                Supports both 2D Cartesian (IAXIS=0) and Axisymmetric (IAXIS=1)
!===============================================================================
MODULE RESIDUAL
    USE SOLVER_STATE, ONLY: NN, NELES, NC1, NC2, NC3, NC4, SC1X, SC2X, SC3X, &
        SC4X, SC1Y, SC2Y, SC3Y, SC4Y, VOL, F1, F2, F3, F4, F5, G1, G2, G3,   &
        G4, G5, R1, R2, R3, R4, R5, SOURCE3, SOURCE5, DT, D241, D242, D243, &
        D244, D245, VISU, VISV, VIST, VISR, COLD1, COLD2, COLD3, COLD4,      &
        COLD5, C1, C2, C3, C4, C5, NCC1, NCC2, NCC3, NCC4, SX1, SX2, SX3,   &
        SX4, SY1, SY2, SY3, SY4, F1I, F2I, F3I, F4I, F5I, G1I, G2I, G3I,    &
        G4I, G5I, DTBYVOL, ALPHA1, CFLTURB, I
    USE CONSTANTS, ONLY: CFLTURB_VAL
    IMPLICIT NONE

    INTEGER :: IAXIS  ! 0 = 2D Cartesian, 1 = Axisymmetric

CONTAINS

    !===========================================================================
    ! COMPUTE_RESIDUAL - Compute inviscid flux residual for all elements
    !===========================================================================
    SUBROUTINE COMPUTE_RESIDUAL
        IMPLICIT NONE
        INTEGER :: I

        DO I = 1, NELES
            NCC1=NC1(I); NCC2=NC2(I); NCC3=NC3(I); NCC4=NC4(I)
            SX1=SC1X(I); SX2=SC2X(I); SX3=SC3X(I); SX4=SC4X(I)
            SY1=SC1Y(I); SY2=SC2Y(I); SY3=SC3Y(I); SY4=SC4Y(I)
            F1I=F1(I); F2I=F2(I); F3I=F3(I); F4I=F4(I); F5I=F5(I)
            G1I=G1(I); G2I=G2(I); G3I=G3(I); G4I=G4(I); G5I=G5(I)

            R1(I) = HALF*((F1I+F1(NCC1))*SX1+(G1I+G1(NCC1))*SY1 + &
                          (F1I+F1(NCC2))*SX2+(G1I+G1(NCC2))*SY2 + &
                          (F1I+F1(NCC3))*SX3+(G1I+G1(NCC3))*SY3 + &
                          (F1I+F1(NCC4))*SX4+(G1I+G1(NCC4))*SY4)
            R2(I) = HALF*((F2I+F2(NCC1))*SX1+(G2I+G2(NCC1))*SY1 + &
                          (F2I+F2(NCC2))*SX2+(G2I+G2(NCC2))*SY2 + &
                          (F2I+F2(NCC3))*SX3+(G2I+G2(NCC3))*SY3 + &
                          (F2I+F2(NCC4))*SX4+(G2I+G2(NCC4))*SY4)
            R3(I) = HALF*((F3I+F3(NCC1))*SX1+(G3I+G3(NCC1))*SY1 + &
                          (F3I+F3(NCC2))*SX2+(G3I+G3(NCC2))*SY2 + &
                          (F3I+F3(NCC3))*SX3+(G3I+G3(NCC3))*SY3 + &
                          (F3I+F3(NCC4))*SX4+(G3I+G3(NCC4))*SY4)
            R4(I) = HALF*((F4I+F4(NCC1))*SX1+(G4I+G4(NCC1))*SY1 + &
                          (F4I+F4(NCC2))*SX2+(G4I+G4(NCC2))*SY2 + &
                          (F4I+F4(NCC3))*SX3+(G4I+G4(NCC3))*SY3 + &
                          (F4I+F4(NCC4))*SX4+(G4I+G4(NCC4))*SY4)
            R5(I) = HALF*((F5I+F5(NCC1))*SX1+(G5I+G5(NCC1))*SY1 + &
                          (F5I+F5(NCC2))*SX2+(G5I+G5(NCC2))*SY2 + &
                          (F5I+F5(NCC3))*SX3+(G5I+G5(NCC3))*SY3 + &
                          (F5I+F5(NCC4))*SX4+(G5I+G5(NCC4))*SY4) &
                          - VOL(I)*SOURCE5(I)

            ! Axisymmetric: add radial source term to R3 (momentum in Y)
            IF (IAXIS .EQ. 1) THEN
                R3(I) = R3(I) - VOL(I)*SOURCE3(I)
            END IF
        END DO

    END SUBROUTINE COMPUTE_RESIDUAL

    !===========================================================================
    ! UPDATE_SOLUTION - Advance conserved variables using RK stage
    !===========================================================================
    SUBROUTINE UPDATE_SOLUTION(ALPHA_COEFF)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: ALPHA_COEFF

        CFLTURB = CFLTURB_VAL

        DO I = 1, NELES
            DTBYVOL = DT(I)/VOL(I)
            C1(I) = COLD1(I) - ALPHA_COEFF*DTBYVOL*(R1(I) - D241(I))
            C2(I) = COLD2(I) - ALPHA_COEFF*DTBYVOL*(R2(I) - VISU(I) - D242(I))
            C3(I) = COLD3(I) - ALPHA_COEFF*DTBYVOL*(R3(I) - VISV(I) - D243(I))
            C4(I) = COLD4(I) - ALPHA_COEFF*DTBYVOL*(R4(I) - VIST(I) - D244(I))
            C5(I) = COLD5(I) - ALPHA_COEFF*DTBYVOL*CFLTURB*(R5(I) - VISR(I) - D245(I))
            C5(I) = C1(I) * DMAX1(RTOLERANCE, C5(I)/C1(I))
        END DO

    END SUBROUTINE UPDATE_SOLUTION

END MODULE RESIDUAL
