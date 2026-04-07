!===============================================================================
! viscous.f90 - Viscous flux computation
!               Supports both 2D Cartesian (IAXIS=0) and Axisymmetric (IAXIS=1)
!===============================================================================
MODULE VISCOUS
    USE SOLVER_STATE, ONLY: NN, NODS, NELES, NGHOSTS, NOD, NC1, NC2, NC3, NC4,&
        NPARENT, NGHOST, NTYPE, NSIDE, X, Y, XCC, YCC, SC1X, SC2X, SC3X,     &
        SC4X, SC1Y, SC2Y, SC3Y, SC4Y, U, V, P, T, C1, C5, UN, VN, TN, RN,    &
        RHON, VISU, VISV, VIST, VISR, ITEST, DIV23, ALPHA, BETA, CMU, CONST1, &
        CONST2, CONST3, N1, N2, N3, N4, NG, NP, NT, NS, NCC1, NCC2, NCC3,    &
        NCC4, XA, XB, XC, XD, YA, YB, YC, YD, UA, UB, UC, UD, VA, VB, VC,   &
        VD, TA, TB, TC, TD, RHOA, RHOB, RHOC, RHOD, RA, RB, RC, RD, DNR,    &
        ANR, DUX, DUY, DVX, DVY, DTX, DTY, DRX, DRY, UAVG, VAVG, TAVG,       &
        RHOAVG, RAVG, DIVVEL, AMU, AK, AMOL, PRANDTL, PRANDTLT, CHI, FMU,    &
        AMUT, AKT, AMUEFF, AKEFF, TXX, TYY, TXY, QXX, QYY, SX, SY, UI, VI,   &
        TI, RHOI, RI, UNG, VNG, TNG, RHONG, RNG, DFLOATITEST, RTOLERANCE,    &
        RTRED, YAVG, VBYY
    USE CONSTANTS, ONLY: DIV23_VAL, ALPHA_TURB, BETA_TURB, CMU_VAL,           &
        CONST1_VAL, CONST2_VAL, CONST3_VAL, SUTH_COEF, QUARTER, YAXIS_THRESH
    IMPLICIT NONE

    INTEGER :: IAXIS  ! 0 = 2D Cartesian, 1 = Axisymmetric

CONTAINS

    !===========================================================================
    ! VISCOUS - Compute viscous fluxes for all element faces
    !===========================================================================
    SUBROUTINE COMPUTE_VISCOUS
        IMPLICIT NONE

        DIV23 = DIV23_VAL
        ALPHA = ALPHA_TURB
        BETA  = BETA_TURB
        CONST1 = CONST1_VAL
        CONST2 = CONST2_VAL
        CONST3 = CONST3_VAL
        CMU = CMU_VAL

        UN = 0.0D0; VN = 0.0D0; TN = 0.0D0; RN = 0.0D0; RHON = 0.0D0
        VISU = 0.0D0; VISV = 0.0D0; VIST = 0.0D0; VISR = 0.0D0
        ITEST = 0

        ! Nodal averaging from elements
        CALL NODAL_AVERAGE_ELEMENTS

        ! Nodal averaging from ghost cells
        CALL NODAL_AVERAGE_GHOSTS

        ! Finalize nodal averages
        DO I = 1, NODS
            DFLOATITEST = DFLOAT(ITEST(I))
            UN(I)  = UN(I)  / DFLOATITEST
            VN(I)  = VN(I)  / DFLOATITEST
            TN(I)  = TN(I)  / DFLOATITEST
            RHON(I)= RHON(I)/ DFLOATITEST
            RN(I)  = RN(I)  / DFLOATITEST
            RN(I)  = DMAX1(RTOLERANCE, RN(I))
        END DO

        ! Face loop: compute viscous fluxes
        DO I = 1, NELES
            N1 = NOD(I,1); N2 = NOD(I,2); N3 = NOD(I,3); N4 = NOD(I,4)
            NCC1 = NC1(I); NCC2 = NC2(I); NCC3 = NC3(I); NCC4 = NC4(I)

            CALL VISCOUS_FACE(1, N1, N2, NCC1)
            CALL VISCOUS_FACE(2, N2, N3, NCC2)
            CALL VISCOUS_FACE(3, N3, N4, NCC3)
            CALL VISCOUS_FACE(4, N4, N1, NCC4)
        END DO

    END SUBROUTINE COMPUTE_VISCOUS

    !---------------------------------------------------------------------------
    ! Nodal averaging: contributions from elements
    !---------------------------------------------------------------------------
    SUBROUTINE NODAL_AVERAGE_ELEMENTS
        IMPLICIT NONE
        DO I = 1, NELES
            N1 = NOD(I,1); N2 = NOD(I,2); N3 = NOD(I,3); N4 = NOD(I,4)
            UI = U(I); VI = V(I); TI = T(I); RHOI = C1(I); RI = C5(I)/RHOI
            UN(N1)=UN(N1)+UI; VN(N1)=VN(N1)+VI; TN(N1)=TN(N1)+TI
            RHON(N1)=RHON(N1)+RHOI; RN(N1)=RN(N1)+RI
            UN(N2)=UN(N2)+UI; VN(N2)=VN(N2)+VI; TN(N2)=TN(N2)+TI
            RHON(N2)=RHON(N2)+RHOI; RN(N2)=RN(N2)+RI
            UN(N3)=UN(N3)+UI; VN(N3)=VN(N3)+VI; TN(N3)=TN(N3)+TI
            RHON(N3)=RHON(N3)+RHOI; RN(N3)=RN(N3)+RI
            UN(N4)=UN(N4)+UI; VN(N4)=VN(N4)+VI; TN(N4)=TN(N4)+TI
            RHON(N4)=RHON(N4)+RHOI; RN(N4)=RN(N4)+RI
            ITEST(N1)=ITEST(N1)+1; ITEST(N2)=ITEST(N2)+1
            ITEST(N3)=ITEST(N3)+1; ITEST(N4)=ITEST(N4)+1
        END DO
    END SUBROUTINE NODAL_AVERAGE_ELEMENTS

    !---------------------------------------------------------------------------
    ! Nodal averaging: contributions from ghost cells with wall/free-slip BCs
    !---------------------------------------------------------------------------
    SUBROUTINE NODAL_AVERAGE_GHOSTS
        IMPLICIT NONE
        DO I = 1, NGHOSTS
            NP = NPARENT(I)
            NG = NGHOST(I)
            NT = NTYPE(I)
            NS = NSIDE(I)
            IF (NS .EQ. 1) THEN; N1=NOD(NP,1); N2=NOD(NP,2)
            ELSE IF (NS .EQ. 2) THEN; N1=NOD(NP,2); N2=NOD(NP,3)
            ELSE IF (NS .EQ. 3) THEN; N1=NOD(NP,3); N2=NOD(NP,4)
            ELSE IF (NS .EQ. 4) THEN; N1=NOD(NP,4); N2=NOD(NP,1)
            END IF

            UN(N1)=UN(N1)+U(NG); VN(N1)=VN(N1)+V(NG); TN(N1)=TN(N1)+T(NG)
            RHON(N1)=RHON(N1)+C1(NG); RN(N1)=RN(N1)+C5(NG)/C1(NG)
            UN(N2)=UN(N2)+U(NG); VN(N2)=VN(N2)+V(NG); TN(N2)=TN(N2)+T(NG)
            RHON(N2)=RHON(N2)+C1(NG); RN(N2)=RN(N2)+C5(NG)/C1(NG)

            IF (NT .EQ. 4) THEN
                VN(N1)=0.0D0; VN(N2)=0.0D0
            ELSE IF (NT .EQ. 5) THEN
                UN(N1)=0.0D0; UN(N2)=0.0D0; VN(N1)=0.0D0; VN(N2)=0.0D0
                RN(N1)=RTOLERANCE*RTRED; RN(N2)=RTOLERANCE*RTRED
            END IF
            ITEST(N1)=ITEST(N1)+1; ITEST(N2)=ITEST(N2)+1
        END DO
    END SUBROUTINE NODAL_AVERAGE_GHOSTS

    !---------------------------------------------------------------------------
    ! Viscous flux computation for a single face of element I
    ! FACE: face number (1-4), NN1/NN2: node indices on face, NCC: neighbor
    !---------------------------------------------------------------------------
    SUBROUTINE VISCOUS_FACE(FACE, NN1, NN2, NCC)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: FACE, NN1, NN2, NCC
        DOUBLE PRECISION :: DNR_LOCAL

        ! Set surface vector based on face
        IF (FACE .EQ. 1) THEN
            SX=SC1X(I); SY=SC1Y(I); XB=X(NN1); YB=Y(NN1); XD=X(NN2); YD=Y(NN2)
        ELSE IF (FACE .EQ. 2) THEN
            SX=SC2X(I); SY=SC2Y(I); XB=X(NN1); YB=Y(NN1); XD=X(NN2); YD=Y(NN2)
        ELSE IF (FACE .EQ. 3) THEN
            SX=SC3X(I); SY=SC3Y(I); XB=X(NN1); YB=Y(NN1); XD=X(NN2); YD=Y(NN2)
        ELSE IF (FACE .EQ. 4) THEN
            SX=SC4X(I); SY=SC4Y(I); XB=X(NN1); YB=Y(NN1); XD=X(NN2); YD=Y(NN2)
        END IF

        XA=XCC(I); YA=YCC(I); XC=XCC(NCC); YC=YCC(NCC)

        UA=U(I); VA=V(I); TA=T(I)
        UB=UN(NN1); VB=VN(NN1); TB=TN(NN1)
        UC=U(NCC); VC=V(NCC); TC=T(NCC)
        UD=UN(NN2); VD=VN(NN2); TD=TN(NN2)

        RHOA=C1(I); RHOB=RHON(NN1); RHOC=C1(NCC); RHOD=RHON(NN2)
        RA=C5(I)/C1(I); RB=RN(NN1); RC=C5(NCC)/C1(NCC); RD=RN(NN2)

        ! Gradient reconstruction
        DNR_LOCAL = (XA-XC)*(YB-YD) - (XB-XD)*(YA-YC)
        DUX = ((UA-UC)*(YB-YD) - (UB-UD)*(YA-YC)) / DNR_LOCAL
        DVX = ((VA-VC)*(YB-YD) - (VB-VD)*(YA-YC)) / DNR_LOCAL
        DTX = ((TA-TC)*(YB-YD) - (TB-TD)*(YA-YC)) / DNR_LOCAL
        DRX = ((RA-RC)*(YB-YD) - (RB-RD)*(YA-YC)) / DNR_LOCAL
        DUY = (-(UA-UC)*(XB-XD) + (UB-UD)*(XA-XC)) / DNR_LOCAL
        DVY = (-(VA-VC)*(XB-XD) + (VB-VD)*(XA-XC)) / DNR_LOCAL
        DTY = (-(TA-TC)*(XB-XD) + (TB-TD)*(XA-XC)) / DNR_LOCAL
        DRY = (-(RA-RC)*(XB-XD) + (RB-RD)*(XA-XC)) / DNR_LOCAL

        UAVG  = QUARTER*(UA+UB+UC+UD)
        VAVG  = QUARTER*(VA+VB+VC+VD)
        TAVG  = QUARTER*(TA+TB+TC+TD)
        RHOAVG= QUARTER*(RHOA+RHOB+RHOC+RHOD)
        RAVG  = QUARTER*(RA+RB+RC+RD)

        ! Divergence: differs between 2D and axisymmetric
        IF (IAXIS .EQ. 1) THEN
            ! Axisymmetric: add radial velocity divergence term
            YAVG = HALF*(YB+YD)
            IF (YAVG .LT. YAXIS_THRESH) THEN
                VBYY = DVY
            ELSE
                VBYY = VAVG / YAVG
            END IF
            DIVVEL = DUX + DVY + VBYY
        ELSE
            DIVVEL = DUX + DVY
        END IF

        ! Sutherland viscosity + Baldwin-Lomax turbulence
        AMU = SUTH_COEF*DSQRT(AMOL)*(TAVG**0.6D0)
        AK = AMU*CP/PRANDTL
        CHI = RHOAVG*RAVG/(AMU*CMU)
        FMU = DTANH(ALPHA*CHI*CHI)/DTANH(BETA*CHI*CHI)
        AMUT = RHOAVG*RAVG*FMU
        AKT = AMUT*CP/PRANDTLT
        AMUEFF = AMU + AMUT
        AKEFF = AK + AKT

        ! Stress tensor and heat flux
        TXX = AMUEFF*(DUX+DUX - DIV23*DIVVEL)
        TYY = AMUEFF*(DVY+DVY - DIV23*DIVVEL)
        TXY = AMUEFF*(DUY+DVX)
        QXX = UAVG*TXX + VAVG*TXY + AKEFF*DTX
        QYY = UAVG*TXY + VAVG*TYY + AKEFF*DTY

        ! Accumulate viscous fluxes (additive for faces 2-4)
        IF (FACE .EQ. 1) THEN
            VISU(I) = TXX*SX + TXY*SY
            VISV(I) = TXY*SX + TYY*SY
            VIST(I) = QXX*SX + QYY*SY
            VISR(I) = AMUEFF*DRX*SX + AMUEFF*DRY*SY
        ELSE
            VISU(I) = VISU(I) + TXX*SX + TXY*SY
            VISV(I) = VISV(I) + TXY*SX + TYY*SY
            VIST(I) = VIST(I) + QXX*SX + QYY*SY
            VISR(I) = VISR(I) + AMUEFF*DRX*SX + AMUEFF*DRY*SY
        END IF

    END SUBROUTINE VISCOUS_FACE

END MODULE VISCOUS
