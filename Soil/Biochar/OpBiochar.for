!=======================================================================
!  OpBiochar, Subroutine
!
!  Purpose: Writes daily and seasonal biochar output to BIOCHAR.OUT.
!
!  Output variables written:
!    YRDOY   - Year and day of year
!    BIOCHC  - Total profile biochar C pool (kg C/ha)
!    BIOCHN  - Total profile biochar N pool (kg N/ha)
!    BCMINC  - Total C mineralised today from biochar (kg C/ha/d)
!    NAPBIO  - Number of biochar applications so far
!    CUMBCC  - Cumulative biochar C applied (kg C/ha)
!    CUMBCN  - Cumulative biochar N applied (kg N/ha)
!    Layer-level C (BC1C..BC20C) - per-layer biochar C (kg C/ha)
!
!-----------------------------------------------------------------------
!  REVISION HISTORY
!  03/23/2026 Written - initial biochar output module for DSSAT-CSM
!-----------------------------------------------------------------------
!  Called  : BIOCHAR
!=======================================================================

      SUBROUTINE OpBiochar (CONTROL, ISWITCH,
     &    BiochC_Total, BiochN_Total, BiochC_L, BiochN_L,
     &    dBiochC, NApBioch, CumBiochC, CumBiochN, NLAYR)

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
      REAL,               INTENT(IN) :: BiochC_Total  !kg C/ha
      REAL,               INTENT(IN) :: BiochN_Total  !kg N/ha
      REAL, DIMENSION(NL),INTENT(IN) :: BiochC_L      !per layer kg C/ha
      REAL, DIMENSION(NL),INTENT(IN) :: BiochN_L      !per layer kg N/ha
      REAL, DIMENSION(NL),INTENT(IN) :: dBiochC       !daily C mineralized
      INTEGER,            INTENT(IN) :: NApBioch      !# applications
      REAL,               INTENT(IN) :: CumBiochC     !cumulative C applied
      REAL,               INTENT(IN) :: CumBiochN     !cumulative N applied
      INTEGER,            INTENT(IN) :: NLAYR

!-----------------------------------------------------------------------
!     Local variables
!-----------------------------------------------------------------------
      CHARACTER*6  ERRKEY
      CHARACTER*30 FILEIO
      PARAMETER (ERRKEY = 'OPBCH ')

      INTEGER DYNAMIC, YRDOY, FROP, DAS
      INTEGER LUNIT
      INTEGER L
      REAL    BCMinC_Total    !Total profile C mineralised today

      LOGICAL FEXIST
      DATA LUNIT / 0 /

      DYNAMIC = CONTROL % DYNAMIC
      YRDOY   = CONTROL % YRDOY
      DAS     = CONTROL % DAS
      FROP    = CONTROL % FROP

!***********************************************************************
!     Seasonal initialization – open output file
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
!         Write file header
          WRITE (LUNIT, '(A)')
     &      '*BIOCHAR CARBON AND NITROGEN DYNAMICS'
        END IF

!       Write run header
        WRITE (LUNIT, '(/,A,I4)') '! Simulation run: ', CONTROL % RUN
        WRITE (LUNIT, '(A)')
     &    '@YEAR DOY   DAS   BIOCHC   BIOCHN   BCMINC'//
     &    ' NAPBIO    CUMBCC    CUMBCN'

        CLOSE (LUNIT)
        RETURN
      END IF

!***********************************************************************
!     Daily output
!***********************************************************************
      IF (DYNAMIC .EQ. OUTPUT) THEN
!       Only write on reporting frequency or application days
        IF (MOD(DAS, FROP) .NE. 0 .AND. NApBioch .EQ. 0) RETURN

        BCMinC_Total = 0.0
        DO L = 1, NLAYR
          BCMinC_Total = BCMinC_Total + dBiochC(L)
        END DO

        OPEN (LUNIT, FILE = 'BIOCHAR.OUT', STATUS = 'OLD',
     &        POSITION = 'APPEND')

        WRITE (LUNIT, 100) YRDOY/1000, MOD(YRDOY,1000), DAS,
     &      BiochC_Total, BiochN_Total, BCMinC_Total,
     &      NApBioch, CumBiochC, CumBiochN

 100    FORMAT (I5,1X,I3,1X,I5,3(1X,F8.2),1X,I6,2(1X,F9.2))

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
     &      '! End-of-season layer biochar C (kg C/ha):'
        WRITE (LUNIT, '(A)') '@LAYER   BIOCHC   BIOCHN'

        DO L = 1, NLAYR
          WRITE (LUNIT, 200) L, BiochC_L(L), BiochN_L(L)
        END DO
 200    FORMAT (I6, 2(1X, F8.2))

        CLOSE (LUNIT)
        RETURN
      END IF

      RETURN
      END SUBROUTINE OpBiochar
