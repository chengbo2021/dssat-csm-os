!=======================================================================
!  BIOCHAR, Subroutine
!
!  Purpose: Simulates biochar dynamics in soil.
!  Biochar is a recalcitrant form of pyrogenic carbon applied as a
!  soil amendment.  This module:
!    - Reads biochar application schedule from FILEIO
!    - Distributes applied biochar through the soil profile
!    - Simulates slow temperature- and moisture-sensitive decomposition
!    - Tracks biochar C and N pools per layer
!    - Releases mineral N proportional to C mineralised
!    - Writes daily output via OpBiochar
!
!  Key science parameters:
!    k_ref  = 0.0001 /day  (reference decomp rate at 25 deg C, opt. moist.)
!             ~ 0.3% per year; typical range 0.1 - 1% /yr (Zimmerman 2010)
!    Q10    = 1.5          (temperature sensitivity, lower than SOM ~2.0)
!    Tf_min = 0.0          (no decomposition below 0 deg C)
!    Wf     = SW/DUL       (water factor, 0 to 1; 1 at field capacity)
!
!-----------------------------------------------------------------------
!  REVISION HISTORY
!  03/23/2026 Written - initial biochar module for DSSAT-CSM
!-----------------------------------------------------------------------
!  Called  : SOIL
!  Calls   : OpBiochar, FIND, ERROR
!=======================================================================

      SUBROUTINE BIOCHAR (CONTROL, ISWITCH,
     &    SOILPROP, ST, SW,                               !Input
     &    BiochData)                                      !Output

!-----------------------------------------------------------------------
      USE ModuleDefs
      IMPLICIT NONE
      EXTERNAL OpBiochar, ERROR, FIND, BC_Apply
      SAVE

!-----------------------------------------------------------------------
!     Interface variables
!-----------------------------------------------------------------------
      TYPE (ControlType), INTENT(IN)  :: CONTROL
      TYPE (SwitchType),  INTENT(IN)  :: ISWITCH
      TYPE (SoilType),    INTENT(IN)  :: SOILPROP
      REAL, DIMENSION(NL), INTENT(IN) :: ST     !Soil temperature (deg C)
      REAL, DIMENSION(NL), INTENT(IN) :: SW     !Volumetric soil water

      TYPE (BiochType),   INTENT(OUT) :: BiochData

!-----------------------------------------------------------------------
!     Local variables
!-----------------------------------------------------------------------
      INTEGER, PARAMETER :: SRFC = 0

      CHARACTER*6  ERRKEY, SECTION
      CHARACTER*30 FILEIO
      CHARACTER*90 CHARTEST

      INTEGER DYNAMIC, YRDOY, LNUM, LUNIO, ERRNUM, LINC, FOUND
      INTEGER L, I, NLAYR
      INTEGER BCDATE_TMP

!     Biochar application schedule (read from FILEIO)
      INTEGER, PARAMETER :: MAXBCAP = 9000
      INTEGER NApSched                      !Number scheduled applications
      INTEGER BCSched_Day(MAXBCAP)          !Application dates
      REAL    BCSched_Amt(MAXBCAP)          !Amount (kg DM/ha)
      REAL    BCSched_Dep(MAXBCAP)          !Incorporation depth (cm)
      REAL    BCSched_CN (MAXBCAP)          !C:N ratio of biochar
      CHARACTER*5 BCSched_Typ(MAXBCAP)      !Feedstock type code
      DATA NApSched /0/

!     Biochar pool state variables
      REAL, DIMENSION(NL) :: BiochC_L      !Layer biochar C (kg C/ha)
      REAL, DIMENSION(NL) :: BiochN_L      !Layer biochar N (kg N/ha)
      REAL BiochC_Total                    !Total profile biochar C
      REAL BiochN_Total                    !Total profile biochar N

!     Decomposition working variables
      REAL, DIMENSION(NL) :: dBiochC       !Daily C decomposed per layer
      REAL, DIMENSION(NL) :: dBiochN       !Daily N released per layer
      REAL, DIMENSION(NL) :: Tfac          !Temperature factor
      REAL, DIMENSION(NL) :: Wfac          !Water factor
      REAL KREF                            !Reference decomp rate (/day)
      PARAMETER (KREF = 0.0001)            !~0.3% /yr at 25 deg C

!     Soil properties
      REAL, DIMENSION(NL) :: DLAYR, DUL, DS

      LOGICAL BIOC_WRITE

      BIOC_WRITE = (ISWITCH % IDETL .NE. 'N')

!-----------------------------------------------------------------------
!     Transfer soil properties
!-----------------------------------------------------------------------
      NLAYR = SOILPROP % NLAYR
      DO L = 1, NLAYR
        DLAYR(L) = SOILPROP % DLAYR(L)
        DUL(L)   = SOILPROP % DUL(L)
        DS(L)    = SOILPROP % DS(L)
      END DO

      DYNAMIC = CONTROL % DYNAMIC
      YRDOY   = CONTROL % YRDOY
      LUNIO   = CONTROL % LUNIO
      FILEIO  = CONTROL % FILEIO

      ERRKEY  = 'BIOCH '

!***********************************************************************
!***********************************************************************
!     Run initialization (once per run)
!***********************************************************************
      IF (DYNAMIC .EQ. RUNINIT) THEN
!-----------------------------------------------------------------------
        DO L = 1, NL
          BiochC_L(L) = 0.0
          BiochN_L(L) = 0.0
          dBiochC(L)  = 0.0
          dBiochN(L)  = 0.0
        END DO
        BiochData % NApBioch = 0
        BiochData % BiochDat = 0
        BiochData % CumBiochC = 0.0
        BiochData % CumBiochN = 0.0
        DO L = 1, NL
          BiochData % BiochC(L) = 0.0
          BiochData % BiochN(L) = 0.0
        END DO
        NApSched = 0

!***********************************************************************
!***********************************************************************
!     Seasonal initialization (once per season)
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. SEASINIT) THEN
!-----------------------------------------------------------------------
!       Reinitialize accumulators each season
        BiochData % NApBioch = 0
        BiochData % BiochDat = 0

!       Reset schedule counter
        NApSched = 0

!       Read biochar application schedule from FILEIO
        OPEN (LUNIO, FILE = FILEIO, STATUS = 'OLD', IOSTAT = ERRNUM)
        IF (ERRNUM .NE. 0) CALL ERROR (ERRKEY, ERRNUM, FILEIO, 0)
        LNUM = 0

        SECTION = '*BIOCH'
        CALL FIND (LUNIO, SECTION, LINC, FOUND)
        LNUM = LNUM + LINC

        IF (FOUND .GT. 0) THEN
!         Read application records
          DO I = 1, MAXBCAP
            READ (LUNIO, '(3X,I7,1X,A90)', ERR=200, END=200)
     &          BCDATE_TMP, CHARTEST
            LNUM = LNUM + 1

            READ (CHARTEST, '(F6.0,1X,F6.1,1X,F6.1,1X,A5)',
     &          IOSTAT=ERRNUM)
     &          BCSched_Amt(I), BCSched_Dep(I),
     &          BCSched_CN(I),  BCSched_Typ(I)
            IF (ERRNUM .NE. 0) GO TO 200

            BCSched_Day(I) = BCDATE_TMP
            BCSched_Amt(I) = MAX(BCSched_Amt(I), 0.0)
            BCSched_Dep(I) = MAX(BCSched_Dep(I), 0.0)
            BCSched_CN(I)  = MAX(BCSched_CN(I),  1.0)
            NApSched       = I
          END DO
        ENDIF

 200    CLOSE (LUNIO)

!       Initialize output file for this season
        IF (BIOC_WRITE) THEN
          CALL OpBiochar (CONTROL, ISWITCH,
     &        0.0, 0.0, BiochC_L, BiochN_L,
     &        dBiochC, 0, 0.0, 0.0, NLAYR)
        END IF

!***********************************************************************
!***********************************************************************
!     Daily rate calculations
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. RATE) THEN
!-----------------------------------------------------------------------
!       Initialize daily change arrays
        DO L = 1, NLAYR
          dBiochC(L) = 0.0
          dBiochN(L) = 0.0
        END DO

!       --- Check for biochar application today ---
        DO I = 1, NApSched
          IF (BCSched_Day(I) .EQ. YRDOY) THEN
            CALL BC_Apply (BiochC_L, BiochN_L, BCSched_Amt(I),
     &          BCSched_Dep(I), BCSched_CN(I), DLAYR, DS, NLAYR,
     &          BiochData)
            BiochData % NApBioch = BiochData % NApBioch + 1
            BiochData % BiochDat = YRDOY
          END IF
        END DO

!       --- Compute temperature and moisture response factors ---
        DO L = 1, NLAYR
!         Temperature factor: exponential with Q10 = 1.5
!         Tf = Q10 ^ ((T-25)/10) = exp((T-25)*ln(1.5)/10)
!         No decomposition below 0 deg C
          IF (ST(L) .GT. 0.0) THEN
            Tfac(L) = EXP((ST(L) - 25.0) * 0.04055)  ! ln(1.5)/10=0.04055
            Tfac(L) = MAX(Tfac(L), 0.0)
          ELSE
            Tfac(L) = 0.0
          END IF

!         Water factor: linear from 0 at wilting point to 1 at field cap.
!         Capped at 1 (saturated conditions similar to field capacity)
          IF (DUL(L) .GT. 1.0E-6) THEN
            Wfac(L) = MIN(SW(L) / DUL(L), 1.0)
            Wfac(L) = MAX(Wfac(L), 0.0)
          ELSE
            Wfac(L) = 0.5
          END IF

!         Daily decomposition (first-order kinetics)
          dBiochC(L) = KREF * Tfac(L) * Wfac(L) * BiochC_L(L)

!         N released proportional to C decomposed and biochar N pool
!         Guard against division by zero
          IF (BiochC_L(L) .GT. 1.0E-6) THEN
            dBiochN(L) = dBiochC(L) * (BiochN_L(L) / BiochC_L(L))
          ELSE
            dBiochN(L) = 0.0
          END IF

!         Constrain losses to available pool
          dBiochC(L) = MIN(dBiochC(L), BiochC_L(L))
          dBiochN(L) = MIN(dBiochN(L), BiochN_L(L))
        END DO

!***********************************************************************
!***********************************************************************
!     Integration of state variables
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. INTEGR) THEN
!-----------------------------------------------------------------------
        DO L = 1, NLAYR
          BiochC_L(L) = BiochC_L(L) - dBiochC(L)
          BiochN_L(L) = BiochN_L(L) - dBiochN(L)

          BiochC_L(L) = MAX(BiochC_L(L), 0.0)
          BiochN_L(L) = MAX(BiochN_L(L), 0.0)

!         Update output data type
          BiochData % BiochC(L) = BiochC_L(L)
          BiochData % BiochN(L) = BiochN_L(L)
        END DO

!***********************************************************************
!***********************************************************************
!     Daily output
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. OUTPUT) THEN
!-----------------------------------------------------------------------
        IF (BIOC_WRITE) THEN
!         Compute profile totals
          BiochC_Total = 0.0
          BiochN_Total = 0.0
          DO L = 1, NLAYR
            BiochC_Total = BiochC_Total + BiochC_L(L)
            BiochN_Total = BiochN_Total + BiochN_L(L)
          END DO

          CALL OpBiochar (CONTROL, ISWITCH,
     &        BiochC_Total, BiochN_Total, BiochC_L, BiochN_L,
     &        dBiochC, BiochData % NApBioch,
     &        BiochData % CumBiochC, BiochData % CumBiochN, NLAYR)
        END IF

!***********************************************************************
!***********************************************************************
!     End-of-season and end-of-run (nothing to do for now)
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. SEASEND .OR. DYNAMIC .EQ. ENDRUN) THEN
!-----------------------------------------------------------------------
        IF (BIOC_WRITE) THEN
          CALL OpBiochar (CONTROL, ISWITCH,
     &        0.0, 0.0, BiochC_L, BiochN_L,
     &        dBiochC, BiochData % NApBioch,
     &        BiochData % CumBiochC, BiochData % CumBiochN, NLAYR)
        END IF

      END IF  !DYNAMIC

      RETURN
      END SUBROUTINE BIOCHAR

!=======================================================================
!  BC_Apply, Subroutine (internal helper)
!
!  Distributes a biochar application through soil layers.
!  Biochar applied to the surface stays on SRFC layer (L=0 proxy
!  mapped to L=1 topsoil) when BCDEP=0, or is proportionally
!  distributed through layers up to the incorporation depth.
!=======================================================================

      SUBROUTINE BC_Apply (BiochC_L, BiochN_L,
     &    BCAMT, BCDEP, BCCN, DLAYR, DS, NLAYR, BiochData)

      USE ModuleDefs
      IMPLICIT NONE

      REAL, DIMENSION(NL), INTENT(INOUT) :: BiochC_L
      REAL, DIMENSION(NL), INTENT(INOUT) :: BiochN_L
      REAL,                INTENT(IN)    :: BCAMT  !kg DM/ha
      REAL,                INTENT(IN)    :: BCDEP  !cm
      REAL,                INTENT(IN)    :: BCCN   !C:N ratio
      REAL, DIMENSION(NL), INTENT(IN)    :: DLAYR  !Layer thicknesses
      REAL, DIMENSION(NL), INTENT(IN)    :: DS     !Layer bottom depths
      INTEGER,             INTENT(IN)    :: NLAYR
      TYPE (BiochType),    INTENT(INOUT) :: BiochData

!     Biochar carbon fraction of dry mass (~50-80% C by weight)
      REAL, PARAMETER :: BC_CF = 0.60   !Carbon fraction of DM

      REAL TotalC, TotalN, LayerFrac, DepthSoFar, TmpDep
      INTEGER L

      TotalC = BCAMT * BC_CF                    !kg C/ha
      IF (BCCN .GT. 0.0) THEN
        TotalN = TotalC / BCCN                  !kg N/ha
      ELSE
        TotalN = 0.0
      END IF

      IF (BCDEP .LE. 0.0) THEN
!       Surface application - place all in layer 1 (topsoil)
        BiochC_L(1) = BiochC_L(1) + TotalC
        BiochN_L(1) = BiochN_L(1) + TotalN
      ELSE
!       Incorporated - distribute proportionally by layer thickness
!       down to the incorporation depth
        TmpDep = 0.0
        DO L = 1, NLAYR
          IF (DS(L) .LE. BCDEP) THEN
!           Layer entirely within incorporation depth
            LayerFrac = DLAYR(L) / BCDEP
          ELSE
!           Partial layer
            IF (L .EQ. 1) THEN
              LayerFrac = MIN(BCDEP, DLAYR(L)) / BCDEP
            ELSE
              DepthSoFar = DS(L) - DLAYR(L)
              IF (DepthSoFar .GE. BCDEP) EXIT  !Below incorporation depth
              LayerFrac = (BCDEP - DepthSoFar) / BCDEP
            END IF
          END IF
          LayerFrac = MIN(MAX(LayerFrac, 0.0), 1.0)
          BiochC_L(L) = BiochC_L(L) + TotalC * LayerFrac
          BiochN_L(L) = BiochN_L(L) + TotalN * LayerFrac

          IF (DS(L) .GE. BCDEP) EXIT
        END DO
      END IF

!     Accumulate cumulative totals
      BiochData % CumBiochC = BiochData % CumBiochC + TotalC
      BiochData % CumBiochN = BiochData % CumBiochN + TotalN

      RETURN
      END SUBROUTINE BC_Apply
