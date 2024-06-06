      REAL FUNCTION half_gen()
*********************************************************
*                                                       *
* This file was generated by HUWFUN.                    *
*                                                       *
*********************************************************
*
*     Ntuple Id:      1    
*     Ntuple Title:   Half mirror. Rcol2,dhcol=   4.   0. f=200. cm
*     Creation:       05/06/2024 10.47.35
*
*********************************************************
*
      LOGICAL         CHAIN
      CHARACTER*128   CFILE
      INTEGER         IDNEVT,NCHEVT,ICHEVT
      REAL            OBS(13)
*
      COMMON /PAWIDN/ IDNEVT,OBS
      COMMON /PAWCHN/ CHAIN, NCHEVT, ICHEVT
      COMMON /PAWCHC/ CFILE
*
*--   Ntuple Variable Declarations
*
      REAL PART(6),RAY0(6),RAYCAP(6),RAYMIR(6),raytst(6),raypmt(6),DBOX
     + ,RBOX,F,RPMT,DCOL,RCOL2,DHCOL,alpha,pere,perpi
      INTEGER IREG,IP,NPART,M,intref
*
      COMMON /PAWCR4/ IREG,IP,NPART,PART,RAY0,RAYCAP,M,RAYMIR,raytst
     + ,intref,raypmt,DBOX,RBOX,F,RPMT,DCOL,RCOL2,DHCOL,alpha,pere,perpi
*

*
*--   Enter user code here
*

      half_gen = 1.
*
      END