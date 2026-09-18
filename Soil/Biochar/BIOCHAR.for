!=======================================================================
!  BIOCHAR, Subroutine
!
!  Purpose: Simulates biochar dynamics in soil.
!  Biochar is a recalcitrant form of pyrogenic carbon applied as a
!  soil amendment.  This module:
!    - Reads biochar application schedule from FILEIO
!    - Looks up feedstock properties from Data/BIOCHAR.CDE
!    - Distributes applied biochar through the soil profile
!    - Simulates two-pool (labile + stable) decomposition
!    - Tracks biochar C and N pools per layer
!    - Releases mineral N proportional to C mineralised
!    - Computes biochar-induced water-retention increase per layer
!    - Writes daily output via OpBiochar
!
!  Two-pool decomposition model:
!    Labile pool  : KREF_L = 0.0010/day at 25 deg C (~5-15% of C)
!                   represents physically accessible and chemically
!                   labile compounds (Zimmerman, 2010)
!    Stable pool  : KREF_S = 0.00005/day at 25 deg C (~85-95% of C)
!                   aromatic fused-ring carbon; very slow mineralisation
!    Q10 = 1.5 for both pools (pyrogenic C less T-sensitive than SOM)
!    Moisture factor: Wf = SW/DUL (0..1, capped at field capacity)
!
!  Water-retention effect:
!    DDUL_BC(L) = WR_COEF * BC_vol_frac(L)
!    BC_vol_frac = BiochMass [kg/m3] / BC_BULK_DENS [kg/m3]
!    WR_COEF = 0.04  (empirical; Omondi et al. 2016 meta-analysis)
!    BC_BULK_DENS = 250 kg/m3 (typical range 200-500)
!
!-----------------------------------------------------------------------
!  REVISION HISTORY
!  03/23/2026 Written - initial biochar module for DSSAT-CSM
!  09/18/2026 Two-pool kinetics, feedstock properties lookup,
!             water-retention effect (Omondi et al. 2016)
!-----------------------------------------------------------------------
!  Called  : SOIL
!  Calls   : OpBiochar, FIND, ERROR, BC_ReadProps, BC_Apply
!=======================================================================

      SUBROUTINE BIOCHAR (CONTROL, ISWITCH,
     &    NH4, SOILPROP, ST, SW,                          !Input
     &    BiochData)                                      !Output

!-----------------------------------------------------------------------
      USE ModuleDefs
      IMPLICIT NONE
      EXTERNAL OpBiochar, ERROR, FIND, BC_Apply, BC_ReadProps
      SAVE

!-----------------------------------------------------------------------
!     Interface variables
!-----------------------------------------------------------------------
      TYPE (ControlType), INTENT(IN)  :: CONTROL
      TYPE (SwitchType),  INTENT(IN)  :: ISWITCH
      REAL, DIMENSION(NL), INTENT(IN) :: NH4    !NH4 pool (kg N/ha/layer)
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

!     Per-application feedstock properties (from BIOCHAR.CDE lookup)
      REAL    BCSched_CF (MAXBCAP)          !C fraction of DM
      REAL    BCSched_FL (MAXBCAP)          !Labile C fraction

!     Biochar pool state variables - two pools
      REAL, DIMENSION(NL) :: BiochCL_L     !Labile C pool (kg C/ha/layer)
      REAL, DIMENSION(NL) :: BiochCS_L     !Stable C pool (kg C/ha/layer)
      REAL, DIMENSION(NL) :: BiochN_L      !Biochar N pool (kg N/ha/layer)
      REAL BiochC_Total                    !Total profile biochar C
      REAL BiochN_Total                    !Total profile biochar N

!     Decomposition working variables
      REAL, DIMENSION(NL) :: dBiochC       !Daily total C decomposed
      REAL, DIMENSION(NL) :: dBiochCL      !Daily labile C decomposed
      REAL, DIMENSION(NL) :: dBiochCS      !Daily stable C decomposed
      REAL, DIMENSION(NL) :: dBiochN       !Daily N released per layer
      REAL, DIMENSION(NL) :: Tfac          !Temperature factor
      REAL, DIMENSION(NL) :: Wfac          !Water factor

!     Decomposition rate constants (/day at 25 deg C, optimal moisture)
      REAL, PARAMETER :: KREF_L = 0.0010   !Labile pool ~0.1%/d = 30%/yr
      REAL, PARAMETER :: KREF_S = 0.00005  !Stable pool ~5e-5/d = 1.8%/yr

!     Water retention parameters
!     Biochar bulk density (kg/m3) for volume fraction calculation
      REAL, PARAMETER :: BC_BULK_DENS = 250.0
!     Empirical DUL increase per unit biochar volume fraction
      REAL, PARAMETER :: WR_COEF = 0.04

!     pH feedback parameters
!     Empirical pH increase per unit biochar mass fraction (Biederman &
!     Harpole 2013 meta-analysis; ~0.3 pH unit per 10 t/ha in top 10 cm)
      REAL, PARAMETER :: PH_COEF    = 40.0
!     Maximum biochar-induced pH increase (prevents runaway in sandy soils)
      REAL, PARAMETER :: PH_MAX_DLT = 2.0

!     NH4 sorption parameters
!     Max NH4 sorption per kg biochar DM (Chen et al. 2019 review)
      REAL, PARAMETER :: SORP_F_NH4 = 0.005  !kg N / kg biochar
!     Linear sorption coefficient (fraction of NH4 sorbed at half capacity)
      REAL, PARAMETER :: KD_NH4     = 0.10
!     Rate constant for approach to sorption equilibrium (/day)
      REAL, PARAMETER :: K_EQ_SORP  = 0.30

!     NH4 sorption state (SAVE'd - persists between calls)
      REAL, DIMENSION(NL) :: SorbNH4_L   !Sorbed NH4 per layer (kg N/ha)

!     Soil properties
      REAL, DIMENSION(NL) :: DLAYR, DUL, DS, BD

!     Water retention and pH working variables
      REAL BiochMass_L, BC_vol_frac, BC_mass_frac
      REAL SorbMax_L, SorbEq_L, dSorb_L

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
        BD(L)    = SOILPROP % BD(L)
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
          BiochCL_L(L) = 0.0
          BiochCS_L(L) = 0.0
          BiochN_L(L)  = 0.0
          dBiochC(L)   = 0.0
          dBiochN(L)   = 0.0
        END DO
        BiochData % NApBioch = 0
        BiochData % BiochDat = 0
        BiochData % CumBiochC = 0.0
        BiochData % CumBiochN = 0.0
        DO L = 1, NL
          BiochData % BiochC(L)  = 0.0
          BiochData % BiochCL(L) = 0.0
          BiochData % BiochCS(L) = 0.0
          BiochData % BiochN(L)  = 0.0
          BiochData % DDUL_BC(L)  = 0.0
          BiochData % DeltaPH(L)  = 0.0
          BiochData % SorbNH4(L)  = 0.0
          BiochData % SorbP(L)    = 0.0
          BiochData % SorbK(L)    = 0.0
          BiochData % dSorbNH4(L) = 0.0
        END DO
        DO L = 1, NL
          SorbNH4_L(L) = 0.0
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

!           Look up feedstock properties from Data/BIOCHAR.CDE
            CALL BC_ReadProps (BCSched_Typ(I),
     &          BCSched_CF(I), BCSched_FL(I), BCSched_CN(I))
          END DO
        ENDIF

 200    CLOSE (LUNIO)

!       Initialize output file for this season
        IF (BIOC_WRITE) THEN
          CALL OpBiochar (CONTROL, ISWITCH,
     &        0.0, 0.0, BiochCL_L, BiochCS_L, BiochN_L,
     &        dBiochC, 0, 0.0, 0.0, BiochData % DeltaPH,
     &        SorbNH4_L, NLAYR)
        END IF

!***********************************************************************
!***********************************************************************
!     Daily rate calculations
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. RATE) THEN
!-----------------------------------------------------------------------
!       Initialize daily change arrays
        DO L = 1, NLAYR
          dBiochC(L)  = 0.0
          dBiochCL(L) = 0.0
          dBiochCS(L) = 0.0
          dBiochN(L)  = 0.0
        END DO

!       --- Check for biochar application today ---
        DO I = 1, NApSched
          IF (BCSched_Day(I) .EQ. YRDOY) THEN
            CALL BC_Apply (BiochCL_L, BiochCS_L, BiochN_L,
     &          BCSched_Amt(I), BCSched_Dep(I), BCSched_CN(I),
     &          BCSched_CF(I),  BCSched_FL(I),
     &          DLAYR, DS, NLAYR, BiochData)
            BiochData % NApBioch = BiochData % NApBioch + 1
            BiochData % BiochDat = YRDOY
          END IF
        END DO

!       --- Compute temperature and moisture response factors ---
        DO L = 1, NLAYR
!         Temperature factor: exponential with Q10 = 1.5
!         Tf = exp((T-25) * ln(1.5)/10); no decomp below 0 deg C
          IF (ST(L) .GT. 0.0) THEN
            Tfac(L) = EXP((ST(L) - 25.0) * 0.04055)  ! ln(1.5)/10
            Tfac(L) = MAX(Tfac(L), 0.0)
          ELSE
            Tfac(L) = 0.0
          END IF

!         Water factor: linear 0 at wilting point, 1 at field capacity
          IF (DUL(L) .GT. 1.0E-6) THEN
            Wfac(L) = MIN(SW(L) / DUL(L), 1.0)
            Wfac(L) = MAX(Wfac(L), 0.0)
          ELSE
            Wfac(L) = 0.5
          END IF

!         Two-pool daily decomposition (first-order kinetics)
          dBiochCL(L) = KREF_L * Tfac(L) * Wfac(L) * BiochCL_L(L)
          dBiochCS(L) = KREF_S * Tfac(L) * Wfac(L) * BiochCS_L(L)

          dBiochCL(L) = MIN(dBiochCL(L), BiochCL_L(L))
          dBiochCS(L) = MIN(dBiochCS(L), BiochCS_L(L))

          dBiochC(L)  = dBiochCL(L) + dBiochCS(L)

!         N released proportional to C decomposed and biochar N pool
          IF ((BiochCL_L(L) + BiochCS_L(L)) .GT. 1.0E-6) THEN
            dBiochN(L) = dBiochC(L) *
     &          (BiochN_L(L) / (BiochCL_L(L) + BiochCS_L(L)))
          ELSE
            dBiochN(L) = 0.0
          END IF
          dBiochN(L) = MIN(dBiochN(L), BiochN_L(L))
        END DO

!***********************************************************************
!***********************************************************************
!     Integration of state variables
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. INTEGR) THEN
!-----------------------------------------------------------------------
        DO L = 1, NLAYR
          BiochCL_L(L) = MAX(BiochCL_L(L) - dBiochCL(L), 0.0)
          BiochCS_L(L) = MAX(BiochCS_L(L) - dBiochCS(L), 0.0)
          BiochN_L(L)  = MAX(BiochN_L(L)  - dBiochN(L),  0.0)

!         Biochar dry mass back-calculated from C pool
          BiochMass_L = 0.0
          IF (BiochCL_L(L) + BiochCS_L(L) .GT. 0.0) THEN
            BiochMass_L = (BiochCL_L(L) + BiochCS_L(L)) / 0.60
          END IF

!         Water-retention effect
!         BC volume fraction = mass [kg/ha] / (100*DLAYR [m3/ha] * rho_BC)
          IF (DLAYR(L) .GT. 0.0) THEN
            BC_vol_frac = (BiochMass_L / (100.0 * DLAYR(L)))
     &                    / BC_BULK_DENS
          ELSE
            BC_vol_frac = 0.0
          END IF

          BiochData % DDUL_BC(L) = WR_COEF * BC_vol_frac

!         pH feedback
!         BC mass fraction (kg/kg) = BC_mass [kg/ha] / soil_mass [kg/ha]
!         Soil mass = BD [g/cm3] * DLAYR [cm] * 1e5  [kg/ha per cm layer]
          IF (BD(L) .GT. 0.0 .AND. DLAYR(L) .GT. 0.0) THEN
            BC_mass_frac = BiochMass_L / (BD(L) * DLAYR(L) * 1.0E5)
            BiochData % DeltaPH(L) = MIN(PH_COEF * BC_mass_frac,
     &                                   PH_MAX_DLT)
          ELSE
            BiochData % DeltaPH(L) = 0.0
          END IF

!         NH4 sorption: linear approach to equilibrium (1-day lag)
          SorbMax_L = SORP_F_NH4 * BiochMass_L
          SorbEq_L  = MIN(SorbMax_L, KD_NH4 * NH4(L))
          dSorb_L   = K_EQ_SORP * (SorbEq_L - SorbNH4_L(L))
          IF (dSorb_L .LT. 0.0)
     &      dSorb_L = MAX(dSorb_L, -SorbNH4_L(L))
          SorbNH4_L(L)            = SorbNH4_L(L) + dSorb_L
          BiochData % SorbNH4(L)  = SorbNH4_L(L)
          BiochData % dSorbNH4(L) = dSorb_L

!         Update output data type
          BiochData % BiochCL(L) = BiochCL_L(L)
          BiochData % BiochCS(L) = BiochCS_L(L)
          BiochData % BiochC(L)  = BiochCL_L(L) + BiochCS_L(L)
          BiochData % BiochN(L)  = BiochN_L(L)
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
            BiochC_Total = BiochC_Total + BiochCL_L(L) + BiochCS_L(L)
            BiochN_Total = BiochN_Total + BiochN_L(L)
          END DO

          CALL OpBiochar (CONTROL, ISWITCH,
     &        BiochC_Total, BiochN_Total, BiochCL_L, BiochCS_L,
     &        BiochN_L, dBiochC, BiochData % NApBioch,
     &        BiochData % CumBiochC, BiochData % CumBiochN,
     &        BiochData % DeltaPH, SorbNH4_L, NLAYR)
        END IF

!***********************************************************************
!***********************************************************************
!     End-of-season and end-of-run
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. SEASEND .OR. DYNAMIC .EQ. ENDRUN) THEN
!-----------------------------------------------------------------------
        IF (BIOC_WRITE) THEN
          CALL OpBiochar (CONTROL, ISWITCH,
     &        0.0, 0.0, BiochCL_L, BiochCS_L, BiochN_L,
     &        dBiochC, BiochData % NApBioch,
     &        BiochData % CumBiochC, BiochData % CumBiochN,
     &        BiochData % DeltaPH, SorbNH4_L, NLAYR)
        END IF

      END IF  !DYNAMIC

      RETURN
      END SUBROUTINE BIOCHAR

!=======================================================================
!  BC_ReadProps, Subroutine
!
!  Looks up feedstock properties from Data/BIOCHAR.CDE.
!  If the code is not found, uses DFLT values.
!  If BCCN_in > 1 (user-supplied), keeps the user value for C:N.
!=======================================================================

      SUBROUTINE BC_ReadProps (BCTYP, CF_DM, FLABIL, BCCN)

      USE OSDefinitions
      IMPLICIT NONE
      CHARACTER*5, INTENT(IN)    :: BCTYP
      REAL,        INTENT(OUT)   :: CF_DM   !C fraction of DM
      REAL,        INTENT(OUT)   :: FLABIL  !Labile C fraction
      REAL,        INTENT(INOUT) :: BCCN    !C:N ratio (kept if user > 1)

      INTEGER         ERRNUM
      CHARACTER*5     RCODE
      REAL            R_CF, R_FL, R_CN
      CHARACTER*80    LINE
      CHARACTER*280   BCPATH
      LOGICAL         FOUND, FEXIST

!     Default fallback (generic biochar parameters)
      CF_DM  = 0.600
      FLABIL = 0.080
      FOUND  = .FALSE.

!     Locate BIOCHAR.CDE: try current dir, then DSSAT installation path
      BCPATH = 'BIOCHAR.CDE'
      INQUIRE (FILE = BCPATH, EXIST = FEXIST)
      IF (.NOT. FEXIST) THEN
        BCPATH = 'Data/BIOCHAR.CDE'
        INQUIRE (FILE = BCPATH, EXIST = FEXIST)
      END IF
      IF (.NOT. FEXIST) THEN
        BCPATH = TRIM(STDPATH) // 'Data/BIOCHAR.CDE'
        INQUIRE (FILE = BCPATH, EXIST = FEXIST)
      END IF
      IF (.NOT. FEXIST) RETURN   !File not found - use defaults silently

      OPEN (UNIT=71, FILE=BCPATH, STATUS='OLD', IOSTAT=ERRNUM)
      IF (ERRNUM .NE. 0) RETURN

      DO WHILE (.NOT. FOUND)
        READ (71, '(A)', IOSTAT=ERRNUM) LINE
        IF (ERRNUM .NE. 0) EXIT
        IF (LINE(1:1) .EQ. '!' .OR. LINE(1:1) .EQ. '*'
     &      .OR. LINE(1:1) .EQ. '@' .OR. LEN_TRIM(LINE) .EQ. 0)
     &    CYCLE
        READ (LINE, '(A5,1X,F6.3,1X,F7.3,1X,F7.1)', IOSTAT=ERRNUM)
     &      RCODE, R_CF, R_FL, R_CN
        IF (ERRNUM .NE. 0) CYCLE
        IF (RCODE .EQ. BCTYP .OR. RCODE .EQ. 'DFLT ') THEN
          CF_DM  = R_CF
          FLABIL = R_FL
!         Only use file C:N if user did not supply one (BCCN <= 1)
          IF (BCCN .LE. 1.0) BCCN = R_CN
          IF (RCODE .EQ. BCTYP) FOUND = .TRUE.
        END IF
      END DO

      CLOSE (71)
      RETURN
      END SUBROUTINE BC_ReadProps

!=======================================================================
!  BC_Apply, Subroutine (internal helper)
!
!  Distributes a biochar application through soil layers into the
!  labile and stable C pools using the feedstock labile fraction.
!=======================================================================

      SUBROUTINE BC_Apply (BiochCL_L, BiochCS_L, BiochN_L,
     &    BCAMT, BCDEP, BCCN, BC_CF, BC_FL,
     &    DLAYR, DS, NLAYR, BiochData)

      USE ModuleDefs
      IMPLICIT NONE

      REAL, DIMENSION(NL), INTENT(INOUT) :: BiochCL_L
      REAL, DIMENSION(NL), INTENT(INOUT) :: BiochCS_L
      REAL, DIMENSION(NL), INTENT(INOUT) :: BiochN_L
      REAL,                INTENT(IN)    :: BCAMT   !kg DM/ha
      REAL,                INTENT(IN)    :: BCDEP   !cm
      REAL,                INTENT(IN)    :: BCCN    !C:N ratio
      REAL,                INTENT(IN)    :: BC_CF   !C fraction of DM
      REAL,                INTENT(IN)    :: BC_FL   !Labile C fraction
      REAL, DIMENSION(NL), INTENT(IN)    :: DLAYR
      REAL, DIMENSION(NL), INTENT(IN)    :: DS
      INTEGER,             INTENT(IN)    :: NLAYR
      TYPE (BiochType),    INTENT(INOUT) :: BiochData

      REAL TotalC, TotalCL, TotalCS, TotalN
      REAL LayerFrac, DepthSoFar
      INTEGER L

      TotalC  = BCAMT * BC_CF
      TotalCL = TotalC * BC_FL          !Labile fraction
      TotalCS = TotalC * (1.0 - BC_FL)  !Stable fraction
      IF (BCCN .GT. 0.0) THEN
        TotalN = TotalC / BCCN
      ELSE
        TotalN = 0.0
      END IF

      IF (BCDEP .LE. 0.0) THEN
!       Surface application – all into layer 1
        BiochCL_L(1) = BiochCL_L(1) + TotalCL
        BiochCS_L(1) = BiochCS_L(1) + TotalCS
        BiochN_L(1)  = BiochN_L(1)  + TotalN
      ELSE
!       Incorporated – distribute by layer thickness to incorporation depth
        DO L = 1, NLAYR
          IF (DS(L) .LE. BCDEP) THEN
            LayerFrac = DLAYR(L) / BCDEP
          ELSE
            IF (L .EQ. 1) THEN
              LayerFrac = MIN(BCDEP, DLAYR(L)) / BCDEP
            ELSE
              DepthSoFar = DS(L) - DLAYR(L)
              IF (DepthSoFar .GE. BCDEP) EXIT
              LayerFrac = (BCDEP - DepthSoFar) / BCDEP
            END IF
          END IF
          LayerFrac = MIN(MAX(LayerFrac, 0.0), 1.0)
          BiochCL_L(L) = BiochCL_L(L) + TotalCL * LayerFrac
          BiochCS_L(L) = BiochCS_L(L) + TotalCS * LayerFrac
          BiochN_L(L)  = BiochN_L(L)  + TotalN  * LayerFrac

          IF (DS(L) .GE. BCDEP) EXIT
        END DO
      END IF

      BiochData % CumBiochC = BiochData % CumBiochC + TotalC
      BiochData % CumBiochN = BiochData % CumBiochN + TotalN

      RETURN
      END SUBROUTINE BC_Apply
