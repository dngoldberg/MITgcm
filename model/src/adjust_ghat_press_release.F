#include "CPP_OPTIONS.h"
#ifdef ALLOW_PRESSURE_RELEASE_CODE
# include "SHELFICE_OPTIONS.h"
#endif


CBOP
C     !ROUTINE: CALC_DIV_GHAT
C     !INTERFACE:
      SUBROUTINE ADJUST_GHAT_PRESS_RELEASE(
     U                cg2d_b,bi,bj,
     I                myThid)
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R CALC_DIV_GHAT
C     | o Form the right hand-side of the surface pressure eqn.
C     *==========================================================*
C     | Right hand side of pressure equation is divergence
C     | of veclocity tendency (GHAT) term along with a relaxation
C     | term equal to the barotropic flow field divergence.
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
#ifdef ALLOW_ADDFLUID
# include "FFIELDS.h"
#endif
#ifdef ALLOW_SHELFICE_GROUNDED_ICE
# include "SHELFICE.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     bi, bj  :: tile indices
C     k       :: Index of layer.
C     cg2d_b  :: Conjugate Gradient 2-D solver : Right-hand side vector
C     cg3d_b  :: Conjugate Gradient 3-D solver : Right-hand side vector
C     myThid  :: Instance number for this call of CALC_DIV_GHAT
      INTEGER bi,bj
      _RL cg2d_b(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER myThid

C     !LOCAL VARIABLES:
C     == Local variables ==
C     i,j    :: Loop counters
C     xA, yA :: Cell vertical face areas
C     pf     :: Intermediate array for building RHS source term.
      INTEGER i,j, k, km1
      _RL pf (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#ifdef ALLOW_PRESSURE_RELEASE_CODE
      _RL     curPrelVisc
      _RL depth, eff_depth, pinit, phiLow, phiLowM1
#ifdef ALLOW_SHELFICE_GROUNDED_ICE
      _RL     avgGrd
#endif
#endif
CEOP

#ifdef ALLOW_PRESSURE_RELEASE_CODE
#ifndef ALLOW_NONHYDROSTATIC

      IF (implicDiv2Dflow.EQ.1) THEN

C     Calculate vertical face areas
!          xA(i,j) = _dyG(i,j,bi,bj)*depthColW(i,j,bi,bj)*rhoFacC(k)
!          yA(i,j) = _dxG(i,j,bi,bj)*depthColS(i,j,bi,bj)*rhoFacC(k)

C--   Pressure equation source term
C     Term is the vertical integral of the divergence of the
C     time tendency terms along with a relaxation term that
C     pulls div(U) + dh/dt back toward zero.

      _EXCH_XY_RL (phi0surf,mythid)

      DO j=1,sNy+1
       DO i=1,sNx+1
        pf(i,j)=0.0
       ENDDO
      ENDDO
      DO j=1,sNy
       DO i=1,sNx+1
        IF (maskInC(i,j,bi,bj)*maskInC(i-1,j,bi,bj).eq.1.) THEN
!         pf(i,j) = xA(i,j)*gU(i,j,k,bi,bj) / deltaTMom
         depth = depthColW(i,j,bi,bj)
         IF (depth.lt.cg2dminColumnEps) THEN

          phiLow = phiHydLowC (i,j,bi,bj)
          phiLowM1 = phiHydLowC (i-1,j,bi,bj)
#ifdef ALLOW_SHELFICE_GROUNDED_ICE
          avgGrd = min (grdFactor(i,j,bi,bj),grdFactor(i-1,j,bi,bj))
          curPrelVisc = max((0.55+avgGrd*0.45),0. _d 0)*pReleaseVisc  
#else
          curPrelVisc = pReleaseVisc
#endif
          eff_depth = cg2dminColumnEps - 
     &     2. * cg2dminColumnEps / PI  * 
     &     COS (PI * depth / (2.*cg2dminColumnEps))

          pinit = -1. * curPrelVisc * 
     &       _recip_dxC(i,j,bi,bj) * _dyG(i,j,bi,bj) 
     &       *(phi0surf(i,j,bi,bj)-phi0surf(i-1,j,bi,bj)+
     &         phiLow-phiLowM1)
     &       /deltaTmom

          pf(i,j) = pinit * (eff_depth - depth)

         ENDIF
        ENDIF
       ENDDO
      ENDDO
      DO j=1,sNy
       DO i=1,sNx
        cg2d_b(i,j,bi,bj) = cg2d_b(i,j,bi,bj) +
     &   pf(i+1,j)-pf(i,j)
       ENDDO
      ENDDO

      DO j=1,sNy+1
       DO i=1,sNx+1
        pf(i,j)=0.0
       ENDDO
      ENDDO

      DO j=1,sNy+1
       DO i=1,sNx
        IF (maskInC(i,j,bi,bj)*maskInC(i,j-1,bi,bj).eq.1.) THEN
         depth = depthColS(i,j,bi,bj)
         IF (depth.lt.cg2dminColumnEps) THEN

          phiLow = phiHydLowC (i,j,bi,bj)
          phiLowM1 = phiHydLowC (i,j-1,bi,bj)
#ifdef ALLOW_SHELFICE_GROUNDED_ICE
          avgGrd = min (grdFactor(i,j,bi,bj),grdFactor(i,j-1,bi,bj))
          curPrelVisc = max((0.55+avgGrd*0.45),0. _d 0)*pReleaseVisc  
#else
          curPrelVisc = pReleaseVisc
#endif

          eff_depth = cg2dminColumnEps - 
     &     2. * cg2dminColumnEps / PI  * 
     &     COS (PI * depth / (2.*cg2dminColumnEps))

          pinit = -1. * curPrelVisc * 
     &       _recip_dyC(i,j,bi,bj) * _dxG(i,j,bi,bj) 
     &       *(phi0surf(i,j,bi,bj)-phi0surf(i,j-1,bi,bj)+
     &         phiLow-phiLowM1)
     &       /deltaTmom

          pf(i,j) = pinit * (eff_depth - depth)

         ENDIF
        ENDIF
       ENDDO
      ENDDO

      DO j=1,sNy
       DO i=1,sNx
        cg2d_b(i,j,bi,bj) = cg2d_b(i,j,bi,bj) +
     &   pf(i,j+1)-pf(i,j)
       ENDDO
      ENDDO

      ENDIF
    
#endif
#endif

      RETURN
      END
