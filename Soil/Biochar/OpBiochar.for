!=======================================================================
!  OpBiochar, Subroutine
!
!  Purpose: Writes daily and seasonal biochar output to BIOCHAR.OUT.
!
!  Output variables written:
!    YRDOY   - Year and day of year
!    BIOCHC  - Total profile biochar C pool (kg C/ha)
!    BCLBL   - Total profile labile biochar C (kg C/ha)
!    BCSTB   - Total profile stable biochar C (kg C/ha)
!    BIOCHN  - Total profile biochar N pool (kg N/ha)
!    BCMINC  - Total C mineralised today from biochar (kg C/ha/d)
!    NAPBIO  - Number of biochar applications so far
!    CUMBCC  - Cumulative biochar C applied (kg C/ha)
!    CUMBCN  - Cumulative biochar N applied (kg N/ha)
!
!-----------------------------------------------------------------------
!  REVISION HISTORY
!  03/23/2026 Written - initial biochar output module for DSSAT-CSM
!  09/18/2026 Report labile and stable pools separately
!-----------------------------------------------------------------------
!  Called  : BIOCHAR
!=======================================================================

      SUBROUTINE OpBiochar (CONTROL, ISWITCH,
     &    BiochC_Total, BiochN_Total, BiochCL_L, BiochCS_L, BiochN_L,
     &    dBiochC, NApBioch, CumBiochC, CumBiochN, DeltaPH_L, NLAYR)

!-----------------------------------------------------------------------
      USE ModuleDefs
      IMPLICIT NONE
      EXTERNAL ERROR, GETLUN
      SAVE

!-----------------------------------------------------------------------
!     Interface variables
!-----------------------------------------------------------------------
      TYPE (ControlType), INTENT(IN) :: CONTROL
      TYPE (SwitchType),  INTENT(IN) :: ISWITCH
      REAL,               INTENT(IN) :: BiochC_Total  !kg C/ha total
      REAL,               INTENT(IN) :: BiochN_Total  !kg N/ha
      REAL, DIMENSION(NL),INTENT(IN) :: BiochCL_L     !labile C/layer
      REAL, DIMENSION(NL),INTENT(IN) :: BiochCS_L     !stable C/layer
      REAL, DIMENSION(NL),INTENT(IN) :: BiochN_L      !N pool/layer
      REAL, DIMENSION(NL),INTENT(IN) :: dBiochC       !daily C mineraliz.
      INTEGER,            INTENT(IN) :: NApBioch
      REAL,               INTENT(IN) :: CumBiochC
      REAL,               INTENT(IN) :: CumBiochN
      REAL, DIMENSION(NL),INTENT(IN) :: DeltaPH_L    !pH increment/layer
      INTEGER,            INTENT(IN) :: NLAYR

!-----------------------------------------------------------------------
!     Local variables
!-----------------------------------------------------------------------
      CHARACTER*6  ERRKEY
      PARAMETER (ERRKEY = 'OPBCH ')

      INTEGER DYNAMIC, YRDOY, FROP, DAS
      INTEGER LUNIT
      INTEGER L
      REAL    BCMinC_Total
      REAL    BCLbl_Total, BCStb_Total, DPH1

      LOGICAL FEXIST
      DATA LUNIT / 0 /

      DYNAMIC = CONTROL % DYNAMIC
      YRDOY   = CONTROL % YRDOY
      DAS     = CONTROL % DAS
      FROP    = CONTROL % FROP

!***********************************************************************
!     Seasonal initialization - open output file
!***********************************************************************
      IF (DYNAMIC .EQ. SEASINIT) THEN
        IF (LUNIT .EQ. 0) THEN
          CALL GETLUN ('BIOCH', LUNIT)
        END IF

        INQUIRE (FILE = 'BIOCHAR.OUT', EXIST = FEXIST)
        IF (FEXIST) THEN
          OPEN (LUNIT, FILE = 'BIOCHAR.OUT', STATUS = 'OLD',
     &          POSITION = 'APPEND')
        ELSE
          OPEN (LUNIT, FILE = 'BIOCHAR.OUT', STATUS = 'NEW')
          WRITE (LUNIT, '(A)')
     &      '*BIOCHAR CARBON AND NITROGEN DYNAMICS'
        END IF

        WRITE (LUNIT, '(/,A,I4)') '! Simulation run: ', CONTROL % RUN
        WRITE (LUNIT, '(A)')
     &    '@YEAR DOY   DAS   BIOCHC    BCLBL    BCSTB'//
     &    '   BIOCHN   BCMINC  DPH_L1 NAPBIO    CUMBCC    CUMBCN'

        CLOSE (LUNIT)
        RETURN
      END IF

!***********************************************************************
!     Daily output
!***********************************************************************
      IF (DYNAMIC .EQ. OUTPUT) THEN
        IF (MOD(DAS, FROP) .NE. 0 .AND. NApBioch .EQ. 0) RETURN

        BCMinC_Total = 0.0
        BCLbl_Total  = 0.0
        BCStb_Total  = 0.0
        DPH1 = DeltaPH_L(1)
        DO L = 1, NLAYR
          BCMinC_Total = BCMinC_Total + dBiochC(L)
          BCLbl_Total  = BCLbl_Total  + BiochCL_L(L)
          BCStb_Total  = BCStb_Total  + BiochCS_L(L)
        END DO

        OPEN (LUNIT, FILE = 'BIOCHAR.OUT', STATUS = 'OLD',
     &        POSITION = 'APPEND')

        WRITE (LUNIT, 100) YRDOY/1000, MOD(YRDOY,1000), DAS,
     &      BiochC_Total, BCLbl_Total, BCStb_Total,
     &      BiochN_Total, BCMinC_Total, DPH1,
     &      NApBioch, CumBiochC, CumBiochN

 100    FORMAT (I5,1X,I3,1X,I5,6(1X,F8.2),1X,I6,2(1X,F9.2))

        CLOSE (LUNIT)
        RETURN
      END IF

!***********************************************************************
!     End of season/run: write layer-level snapshot
!***********************************************************************
      IF (DYNAMIC .EQ. SEASEND .OR. DYNAMIC .EQ. ENDRUN) THEN
        OPEN (LUNIT, FILE = 'BIOCHAR.OUT', STATUS = 'OLD',
     &        POSITION = 'APPEND')

        WRITE (LUNIT, '(/,A)')
     &      '! End-of-season layer biochar pools (kg C/ha or kg N/ha):'
        WRITE (LUNIT, '(A)') '@LAYER    BCLBL    BCSTB   BIOCHC   BIOCHN'

        DO L = 1, NLAYR
          WRITE (LUNIT, 200) L,
     &        BiochCL_L(L), BiochCS_L(L),
     &        BiochCL_L(L)+BiochCS_L(L), BiochN_L(L)
        END DO
 200    FORMAT (I6, 4(1X, F8.2))

        CLOSE (LUNIT)
        RETURN
      END IF

      RETURN
      END SUBROUTINE OpBiochar
