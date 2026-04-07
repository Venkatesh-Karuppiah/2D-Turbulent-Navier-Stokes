!===============================================================================
! source_terms.f90 - Turbulence source term computation
!                    Supports both 2D Cartesian (IAXIS=0) and Axisymmetric (IAXIS=1)
!===============================================================================
MODULE SOURCE_TERMS
    USE SOLVER_STATE, ONLY: NN, NELES, NC1, NC2, NC3, NC4, XCC, YCC, U, V,    &
        C1, C5, T, P, SOURCE3, SOURCE5, ITEST, NCC1, NCC2, NCC3, NCC4, XA,   &
        XB, XC, XD, YA, YB, YC, YD, UA, UB, UC, UD, VA, VB, VC, VD, RA, RB,  &
        RC, RD, QA, QB, QC, QD, DNR, ANR, DQX, DQY, DRX, DRY, DUX, DUY,     &
        DVX, DVY, RHO, RI, AMU, CHI, FMU, AMUT, PHI, DEE, DIVVEL, PK, TTT,   &
        I, AMOL, NPSEUDO, PRANDTL, PRANDTLT, DIV23, ALPHA, BETA, CONST1,     &
        CONST2, CONST3, CMU, VBYY
    USE CONSTANTS, ONLY: DIV23_VAL, ALPHA_TURB, BETA_TURB, CMU_VAL,           &
        CONST1_VAL, CONST2_VAL, CONST3_VAL, SUTH_COEF, YAXIS_THRESH, AMU_THRESH
    IMPLICIT NONE

    INTEGER :: IAXIS  ! 0 = 2D Cartesian, 1 = Axisymmetric

CONTAINS

    !===========================================================================
    ! SOURCETERMS - Compute turbulence model source terms
    !===========================================================================
    SUBROUTINE COMPUTE_SOURCE_TERMS
        IMPLICIT NONE

        DIV23 = DIV23_VAL
        ALPHA = ALPHA_TURB
        BETA  = BETA_TURB
        CONST1 = CONST1_VAL
        CONST2 = CONST2_VAL
        CONST3 = CONST3_VAL
        CMU = CMU_VAL

        SOURCE3 = 0.0D0
        SOURCE5 = 0.0D0
        ITEST = 0

        DO I = 1, NELES
            NCC1 = NC1(I); NCC2 = NC2(I); NCC3 = NC3(I); NCC4 = NC4(I)
            XA=XCC(NCC1); XB=XCC(NCC2); XC=XCC(NCC3); XD=XCC(NCC4)
            YA=YCC(NCC1); YB=YCC(NCC2); YC=YCC(NCC3); YD=YCC(NCC4)
            UA=U(NCC1); UB=U(NCC2); UC=U(NCC3); UD=U(NCC4)
            VA=V(NCC1); VB=V(NCC2); VC=V(NCC3); VD=V(NCC4)
            RA=C5(NCC1)/C1(NCC1); RB=C5(NCC2)/C1(NCC2)
            RC=C5(NCC3)/C1(NCC3); RD=C5(NCC4)/C1(NCC4)

            QA=DSQRT(UA*UA+VA*VA); QB=DSQRT(UB*UB+VB*VB)
            QC=DSQRT(UC*UC+VC*VC); QD=DSQRT(UD*UD+VD*VD)

            ! Gradient reconstruction
            DNR = (XA-XC)*(YB-YD) - (XB-XD)*(YA-YC)
            DQX = ((QA-QC)*(YB-YD)-(QB-QD)*(YA-YC))/DNR
            DRX = ((RA-RC)*(YB-YD)-(RB-RD)*(YA-YC))/DNR
            DUX = ((UA-UC)*(YB-YD)-(UB-UD)*(YA-YC))/DNR
            DVX = ((VA-VC)*(YB-YD)-(VB-VD)*(YA-YC))/DNR
            DQY = (-(QA-QC)*(XB-XD)+(QB-QD)*(XA-XC))/DNR
            DRY = (-(RA-RC)*(XB-XD)+(RB-RD)*(XA-XC))/DNR
            DUY = (-(UA-UC)*(XB-XD)+(UB-UD)*(XA-XC))/DNR
            DVY = (-(VA-VC)*(XB-XD)+(VB-VD)*(XA-XC))/DNR

            RHO = C1(I)
            RI = C5(I)/RHO
            AMU = SUTH_COEF*DSQRT(AMOL)*(T(I)**0.6D0)
            CHI = RHO*RI/(CMU*AMU)
            FMU = DTANH(ALPHA*CHI*CHI)/DTANH(BETA*CHI*CHI)
            FR2 = DTANH(DSQRT(CHI))/DTANH(CHI*DSQRT(2.0D0*DSQRT(CMU)))
            AMUT = FMU*RHO*RI
            PHI = DQX*DRX + DQY*DRY

            DEE = 0.0D0
            IF (PHI .GT. 0.0D0) DEE = DRX*DRX + DRY*DRY

            ! Divergence: differs between 2D and axisymmetric
            IF (IAXIS .EQ. 1) THEN
                IF (YCC(I) .LT. YAXIS_THRESH) THEN
                    VBYY = DVY
                ELSE
                    VBYY = V(I) / YCC(I)
                END IF
                DIVVEL = DUX + DVY + VBYY

                ! Turbulence production with axisymmetric terms
                PK = AMUT/RHO*(DUX*DUX + DVY*DVY + VBYY*VBYY + &
                     DUX*DUX + DVY*DVY + VBYY*VBYY + &
                     (DUY+DVX)*(DUY+DVX) - DIV23*DIVVEL*DIVVEL)
                PK = DMAX1(PK, 0.0D0)
                SOURCE5(I) = (RHO*(CONST1-CONST2*FR2)*DSQRT(RI*PK) - &
                             RHO*CONST3*DEE) * (AMUT/DMAX1(AMUT, AMU_THRESH*AMU))**NPSEUDO

                ! Radial momentum source term (axisymmetric only)
                TTT = (AMU+AMUT)*(VBYY+VBYY - DIV23*DIVVEL)
                SOURCE3(I) = (P(I)-TTT)/YCC(I)
            ELSE
                DIVVEL = DUX + DVY
                PK = AMUT/RHO*(DUX*DUX + DVY*DVY + DUX*DUX + DVY*DVY + &
                     (DUY+DVX)*(DUY+DVX) - DIV23*DIVVEL*DIVVEL)
                PK = DMAX1(PK, 0.0D0)
                SOURCE5(I) = (RHO*(CONST1-CONST2*FR2)*DSQRT(RI*PK) - &
                             RHO*CONST3*DEE) * (AMUT/DMAX1(AMUT, AMU_THRESH*AMU))**NPSEUDO
            END IF
        END DO

    END SUBROUTINE COMPUTE_SOURCE_TERMS

END MODULE SOURCE_TERMS
