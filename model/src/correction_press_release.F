#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"
#ifdef ALLOW_PRESSURE_RELEASE_CODE
# include "SHELFICE_OPTIONS.h"
#endif

CBOP
C     !ROUTINE: CORRECTION_STEP
C     !INTERFACE:
      SUBROUTINE CORRECTION_PRESS_RELEASE( 
     I                       bi, bj, iMin, iMax, jMin, jMax,
     I                       phiSurfX, phiSurfY,
     I                       myTime, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R CORRECTION_STEP
C     | o Corrects the horizontal flow fields with the surface
C     |   pressure (and Non-Hydrostatic pressure).
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "SURFACE.h"
#ifdef ALLOW_SHELFICE_GROUNDED_ICE
# include "SHELFICE.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine Arguments ==
C     bi, bj              :: Tile indices
C     iMin,iMax,jMin,jMax :: Loop counters range
C     k                   :: Level index
C     phiSurfX, phiSurfY  :: Surface Potential gradient
C     myTime              :: Current time in simulation
C     myThid              :: my Thread Id number
      _RL     phiSurfX(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     phiSurfY(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER bi,bj
      INTEGER iMin,iMax,jMin,jMax
      INTEGER k
      _RL     myTime
      INTEGER myThid
#ifdef ALLOW_PRESSURE_RELEASE_CODE
      _RL     curPrelVisc
#ifdef ALLOW_SHELFICE_GROUNDED_ICE
      _RL     avgGrd
#endif
#endif


#ifdef ALLOW_PRESSURE_RELEASE_CODE

C     !LOCAL VARIABLES:
C     == Local variables ==
C     i,j          :: Loop counters
C     psFac, nhFac :: Scaling parameters for supressing gradients
      INTEGER i,j
      _RL     depth, eff_depth
CEOP

C     On/off scaling parameters (including anelastic & deep-model factors)
!      psFac = pfFacMom*implicSurfPress
!     &       *recip_deepFacC(k)*recip_rhoFacC(k)

      DO j=jMin,jMax
       DO i=iMin,iMax
        pReleaseTransX (i,j,bi,bj) = 0. _d 0
        pReleaseTransY (i,j,bi,bj) = 0. _d 0
       ENDDO
      ENDDO

      DO j=jMin,jMax
       DO i=iMin,iMax
        IF (maskInC(i,j,bi,bj)*maskInC(i-1,j,bi,bj).eq.1.) THEN
         depth = depthColW(i,j,bi,bj)
         IF (depth.lt.cg2dminColumnEps) THEN

          eff_depth = cg2dminColumnEps - 
     &     2. * cg2dminColumnEps / PI  * 
     &     COS (PI * depth / (2.*cg2dminColumnEps))
#ifdef ALLOW_SHELFICE_GROUNDED_ICE
            avgGrd = min (grdFactor(i,j,bi,bj),grdFactor(i-1,j,bi,bj))
            curPrelVisc = max((0.55+avgGrd*0.45),0. _d 0)*pReleaseVisc  
#else
            curPrelVisc = pReleaseVisc
#endif

          pReleaseTransX (i,j,bi,bj) = pReleaseTransX (i,j,bi,bj)
     &       - curPrelVisc * 
     &       _recip_dxC(i,j,bi,bj)
     &       *(phiHydLowC(i,j,bi,bj)-phiHydLowC(i-1,j,bi,bj)+
     &         phi0surf(i,j,bi,bj)-phi0surf(i-1,j,bi,bj))
     &       *( eff_depth - depth )

          pReleaseTransX (i,j,bi,bj) = pReleaseTransX (i,j,bi,bj)
     &       - curPrelVisc * 
     &       (eff_depth-depth) * phiSurfX(i,j)

         ENDIF
        ENDIF
       ENDDO
      ENDDO

      DO j=jMin,jMax
       DO i=iMin,iMax
        IF (maskInC(i,j,bi,bj)*maskInC(i,j-1,bi,bj).eq.1.) THEN
         depth = depthColS(i,j,bi,bj)
         IF (depth.lt.cg2dminColumnEps) THEN

          eff_depth = cg2dminColumnEps - 
     &     2. * cg2dminColumnEps / PI  * 
     &     COS (PI * depth / (2.*cg2dminColumnEps))

#ifdef ALLOW_SHELFICE_GROUNDED_ICE
            avgGrd = min (grdFactor(i,j,bi,bj),grdFactor(i,j-1,bi,bj))
            curPrelVisc = max((0.55+avgGrd*0.45),0. _d 0)*pReleaseVisc  
#else
            curPrelVisc = pReleaseVisc
#endif

          pReleaseTransY (i,j,bi,bj) = pReleaseTransY (i,j,bi,bj)
     &       - curPrelVisc * 
     &       _recip_dyC(i,j,bi,bj)
     &       *(phiHydLowC(i,j,bi,bj)-phiHydLowC(i,j-1,bi,bj)+
     &         phi0surf(i,j,bi,bj)-phi0surf(i,j-1,bi,bj))
     &       *( eff_depth - depth )

          pReleaseTransY (i,j,bi,bj) = pReleaseTransY (i,j,bi,bj)
     &       - curPrelVisc * 
     &       (eff_depth-depth) * phiSurfY(i,j)

         ENDIF
        ENDIF
       ENDDO
      ENDDO

      _EXCH_XY_RL (pReleaseTransX,mythid)
      _EXCH_XY_RL (pReleaseTransY,mythid)

#endif

      RETURN
      END
