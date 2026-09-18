! Example *BIOCHAR section for a DSSAT FILEX file
! =====================================================================
! Column layout (same as *FERTILIZERS and *TILLAGE):
!
!   Cols 1-3  : Treatment number (right-justified, 3 chars)
!   Col  4    : Space
!   Cols 5-11 : Application date (YYDDD, 7 chars)
!   Col  12   : Space
!   Cols 13-18: Amount (kg DM ha-1, F6.0)     e.g.  5000. = 5 t/ha
!   Col  19   : Space
!   Cols 20-25: Incorporation depth (cm, F6.1) e.g.    20.0 = 20 cm
!                  0.0 = surface application, no incorporation
!   Col  26   : Space
!   Cols 27-32: C:N ratio of biochar (F6.1)
!                  Use 0.0 to let the model look up the feedstock default
!   Col  33   : Space
!   Cols 34-38: Feedstock code (A5) from Data/BIOCHAR.CDE
!                  WOOD  = Wood biochar (hardwood, high temperature)
!                  SOFWD = Softwood biochar
!                  STRAW = Crop straw biochar
!                  CROP  = Generic crop residue biochar
!                  MANUR = Manure biochar
!                  GREEN = Green waste / compost biochar
!                  SEWAG = Sewage sludge biochar
!                  DFLT  = Default (use if feedstock type unknown)
!
! NOTE: LNBIO (column 15 of the TREATMENTS section) links a treatment
!       to a *BIOCHAR section in the same way LNFER links to *FERTILIZERS.
!       Set LNBIO = 0 (or leave blank) for treatments without biochar.
!
*BIOCHAR
@TRNO BIDATE  BIAMT  BIDEP  BICNR BITYP
  1  2024092  5000.    20.0   0.0 WOOD
  1  2025092  2500.    20.0   0.0 WOOD
  2  2024092 10000.    10.0 200.0 CROP

! Explanation of the examples above:
!  Trt 1: 5 t/ha wood biochar incorporated to 20 cm, C:N from BIOCHAR.CDE
!         Applied in year 2024, day 92 and again in 2025 day 92
!  Trt 2: 10 t/ha crop-residue biochar, 10 cm depth, user C:N = 200
