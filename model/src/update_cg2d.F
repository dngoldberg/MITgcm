#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"
#ifdef ALLOW_SHELFICE
# include "SHELFICE_OPTIONS.h"
#endif

CBOP
C     !ROUTINE: UPDATE_CG2D
C     !INTERFACE:
      SUBROUTINE UPDATE_CG2D( myTime, myIter, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE UPDATE_CG2D
C     | o Update 2d conjugate gradient solver operators
C     |   account for Free-Surf effect on total column thickness
C     *==========================================================*
C     | This routine is based on INI_CG2D, and simplified. It is
C     | used when the non-linear free surface mode is activated
C     | or when bottom depth is part of the control vector.
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SURFACE.h"
#include "CG2D.h"
#ifdef ALLOW_SOLVE4_PS_AND_DRAG
# include "DYNVARS.h"
#endif /* ALLOW_SOLVE4_PS_AND_DRAG */
#ifdef ALLOW_SHELFICE_GROUNDED_ICE
# include "SHELFICE.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myTime :: Current time in simulation
C     myIter :: Current iteration number in simulation
C     myThid :: Thread number for this instance of the routine.
      _RL myTime
      INTEGER myIter
      INTEGER myThid

C     !LOCAL VARIABLES:
C-- Note : compared to "INI_CG2D", no needs to compute again
C   the solver normalisation factor or the solver tolerance
C     === Local variables ===
C     bi,bj  :: tile indices
C     i,j,k  :: Loop counters
C     faceArea :: Temporary used to hold cell face areas.
      INTEGER bi, bj
      INTEGER i, j, k, ks
      _RL     faceArea
      _RL     pW_tmp, pS_tmp
      LOGICAL updatePreCond
#ifdef ALLOW_PRESSURE_RELEASE_CODE
      _RL     curPrelVisc
      _RL depth, eff_depth, depth_fac, pinit
#ifdef ALLOW_DAMP_THINCOLUMN
      _RL drag_fac
#endif
#ifdef ALLOW_SHELFICE_GROUNDED_ICE
      _RL avgGrd
#endif
#endif
CEOP

C--   Decide when to update cg2d Preconditioner :
      IF ( cg2dPreCondFreq.EQ.0 ) THEN
        updatePreCond = .FALSE.
      ELSE
        updatePreCond = ( myIter.EQ.nIter0 )
        IF ( MOD(myIter,cg2dPreCondFreq).EQ.0 ) updatePreCond=.TRUE.
      ENDIF

C--   Initialise laplace operator
C     aW2d: integral in Z Ax/dX
C     aS2d: integral in Z Ay/dY
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          aW2d(i,j,bi,bj) = 0. _d 0
          aS2d(i,j,bi,bj) = 0. _d 0
#ifdef ALLOW_AUTODIFF
          aC2d(i,j,bi,bj) = 0. _d 0
#endif /* ALLOW_AUTODIFF */
         ENDDO
        ENDDO
#ifdef ALLOW_SOLVE4_PS_AND_DRAG
        IF ( selectImplicitDrag.EQ.2 ) THEN
         DO k=1,Nr
          DO j=1,sNy+1
           DO i=1,sNx+1
            faceArea = _dyG(i,j,bi,bj)*drF(k)*deepFacC(k)*rhoFacC(k)
     &                *_hFacW(i,j,k,bi,bj)
#if (defined ALLOW_PRESSURE_RELEASE_CODE && defined ALLOW_DAMP_THINCOLUMN)
           IF (depthColW(i,j,bi,bj).lt.cg2dminColumnEps) THEN       
             drag_fac = _recip_hFacW(i,j,k,bi,bj)**2*recip_drF(k)**2*
     &        thinColDamp * viscArNr(k) * 
     &        1./(1+exp(-10./cg2dminColumnEps*
     &         (-1.)*(depthColW(i,j,bi,bj)-cg2dminColumnEps/2)))
           ELSE
             drag_fac = 0. _d 0
           ENDIF
#endif
            aW2d(i,j,bi,bj) = aW2d(i,j,bi,bj)
     &               + faceArea*dU_psFacX(i,j,k,bi,bj)
     &                         *recip_dxC(i,j,bi,bj)
#if (defined ALLOW_PRESSURE_RELEASE_CODE && defined ALLOW_DAMP_THINCOLUMN)
     &               / (1 + drag_fac*deltaTmom)
#endif
            faceArea = _dxG(i,j,bi,bj)*drF(k)*deepFacC(k)*rhoFacC(k)
     &                *_hFacS(i,j,k,bi,bj)
#if (defined ALLOW_PRESSURE_RELEASE_CODE && defined ALLOW_DAMP_THINCOLUMN)
           IF (depthColS(i,j,bi,bj).lt.cg2dminColumnEps) THEN       
             drag_fac = _recip_hFacS(i,j,k,bi,bj)**2*recip_drF(k)**2*
     &        thinColDamp * viscArNr(k) * 
     &        1./(1+exp(-10./cg2dminColumnEps*
     &         (-1.)*(depthColS(i,j,bi,bj)-cg2dminColumnEps/2)))
           ELSE
             drag_fac = 0. _d 0
           ENDIF
#endif
            aS2d(i,j,bi,bj) = aS2d(i,j,bi,bj)
     &               + faceArea*dV_psFacY(i,j,k,bi,bj)
     &                         *recip_dyC(i,j,bi,bj)
#if (defined ALLOW_PRESSURE_RELEASE_CODE && defined ALLOW_DAMP_THINCOLUMN)
     &               / (1 + drag_fac*deltaTmom)
#endif

           ENDDO
          ENDDO
         ENDDO
        ELSE
#endif /* ALLOW_SOLVE4_PS_AND_DRAG */
         DO k=1,Nr
          DO j=1,sNy+1
           DO i=1,sNx+1
C  deep-model: *deepFacC (faceArea), /deepFacC (recip_dx,y): => no net effect
            faceArea = _dyG(i,j,bi,bj)*drF(k)
     &                *_hFacW(i,j,k,bi,bj)
#if (defined ALLOW_PRESSURE_RELEASE_CODE && defined ALLOW_DAMP_THINCOLUMN)
           IF (depthColW(i,j,bi,bj).lt.cg2dminColumnEps) THEN       
             drag_fac = _recip_hFacW(i,j,k,bi,bj)**2*recip_drF(k)**2*
     &        thinColDamp * viscArNr(k) * 
     &        1./(1+exp(-10./cg2dminColumnEps*
     &         (-1.)*(depthColW(i,j,bi,bj)-cg2dminColumnEps/2)))
           ELSE
             drag_fac = 0. _d 0
           ENDIF
#endif
            aW2d(i,j,bi,bj) = aW2d(i,j,bi,bj)
     &               + faceArea*recip_dxC(i,j,bi,bj)
#if (defined ALLOW_PRESSURE_RELEASE_CODE && defined ALLOW_DAMP_THINCOLUMN)
     &               / (1 + drag_fac*deltaTmom)
#endif
            faceArea = _dxG(i,j,bi,bj)*drF(k)
     &                *_hFacS(i,j,k,bi,bj)
#if (defined ALLOW_PRESSURE_RELEASE_CODE && defined ALLOW_DAMP_THINCOLUMN)
           IF (depthColS(i,j,bi,bj).lt.cg2dminColumnEps) THEN       
             drag_fac = _recip_hFacS(i,j,k,bi,bj)**2*recip_drF(k)**2*
     &        thinColDamp * viscArNr(k) * 
     &        1./(1+exp(-10./cg2dminColumnEps*
     &         (-1.)*(depthColS(i,j,bi,bj)-cg2dminColumnEps/2)))
           ELSE
             drag_fac = 0. _d 0
           ENDIF
#endif
            aS2d(i,j,bi,bj) = aS2d(i,j,bi,bj)
     &               + faceArea*recip_dyC(i,j,bi,bj)
#if (defined ALLOW_PRESSURE_RELEASE_CODE && defined ALLOW_DAMP_THINCOLUMN)
     &               / (1 + drag_fac*deltaTmom)
#endif
           ENDDO
          ENDDO
         ENDDO
#ifdef ALLOW_SOLVE4_PS_AND_DRAG
        ENDIF
#endif /* ALLOW_SOLVE4_PS_AND_DRAG */
        DO j=1,sNy+1
         DO i=1,sNx+1
          aW2d(i,j,bi,bj) = aW2d(i,j,bi,bj)*cg2dNorm
     &                     *implicSurfPress*implicDiv2DFlow
#ifdef ALLOW_OBCS
     &                  *maskInC(i,j,bi,bj)*maskInC(i-1,j,bi,bj)
#endif
          aS2d(i,j,bi,bj) = aS2d(i,j,bi,bj)*cg2dNorm
     &                     *implicSurfPress*implicDiv2DFlow
#ifdef ALLOW_OBCS
     &                  *maskInC(i,j,bi,bj)*maskInC(i,j-1,bi,bj)
#endif
         ENDDO
        ENDDO
#ifdef ALLOW_PRESSURE_RELEASE_CODE
        DO j=1,sNy+1
         DO i=1,sNx+1
          IF (depthColW(i,j,bi,bj).lt.cg2dminColumnEps) THEN 
           IF (maskInC(i,j,bi,bj)*maskInC(i-1,j,bi,bj).eq.1.) THEN
#ifdef ALLOW_SHELFICE_GROUNDED_ICE
            avgGrd = min (grdFactor(i,j,bi,bj),grdFactor(i-1,j,bi,bj)) 
            curPrelVisc = max((0.55+avgGrd*0.45),0. _d 0)*pReleaseVisc 
#else
            curPrelVisc = pReleaseVisc
#endif
            depth = depthColW(i,j,bi,bj)
            eff_depth = cg2dminColumnEps * 
     &       (1.-2./PI * COS (PI*depth / (2.*cg2dminColumnEps)))
            aW2d(i,j,bi,bj) = aW2d(i,j,bi,bj) + 
     &                        curPrelVisc*(eff_depth-depth)
     &                        *_dyG(i,j,bi,bj)*
     &                        recip_dxC(i,j,bi,bj) * cg2dNorm * 
     &                        implicSurfPress*implicDiv2DFlow /
     &                        deltaTmom     
           ENDIF
          ENDIF
          IF (depthColS(i,j,bi,bj).lt.cg2dminColumnEps) THEN 
           IF (maskInC(i,j,bi,bj)*maskInC(i,j-1,bi,bj).eq.1.) THEN
#ifdef ALLOW_SHELFICE_GROUNDED_ICE
            avgGrd = min (grdFactor(i,j,bi,bj),grdFactor(i,j-1,bi,bj))
            curPrelVisc = max((0.55+avgGrd*0.45),0. _d 0)*pReleaseVisc  
#else
            curPrelVisc = pReleaseVisc
#endif
            depth = depthColS(i,j,bi,bj)
            eff_depth = cg2dminColumnEps * 
     &       (1.-2./PI * COS (PI*depth / (2.*cg2dminColumnEps)))
            aS2d(i,j,bi,bj) = aS2d(i,j,bi,bj) + 
     &                        curPrelVisc*(eff_depth-depth)
     &                        *_dxG(i,j,bi,bj)*
     &                        recip_dyC(i,j,bi,bj) * cg2dNorm * 
     &                        implicSurfPress*implicDiv2DFlow /
     &                        deltaTmom
           ENDIF
          ENDIF
         ENDDO
        ENDDO
#endif
C--   compute matrix main diagonal :
        IF ( deepAtmosphere ) THEN
         DO j=1,sNy
          DO i=1,sNx
           ks = kSurfC(i,j,bi,bj)
           aC2d(i,j,bi,bj) = -(
     &       aW2d(i,j,bi,bj) + aW2d(i+1, j ,bi,bj)
     &      +aS2d(i,j,bi,bj) + aS2d( i ,j+1,bi,bj)
     &      +freeSurfFac*cg2dNorm*recip_Bo(i,j,bi,bj)*deepFac2F(ks)
     &                  *rA(i,j,bi,bj)/deltaTMom/deltaTFreeSurf
     &                        )
          ENDDO
         ENDDO
        ELSE
         DO j=1,sNy
          DO i=1,sNx
           aC2d(i,j,bi,bj) = -(
     &       aW2d(i,j,bi,bj) + aW2d(i+1, j ,bi,bj)
     &      +aS2d(i,j,bi,bj) + aS2d( i ,j+1,bi,bj)
     &      +freeSurfFac*cg2dNorm*recip_Bo(i,j,bi,bj)
     &                  *rA(i,j,bi,bj)/deltaTMom/deltaTFreeSurf
     &                        )
          ENDDO
         ENDDO
        ENDIF
C-    end bi,bj loops
       ENDDO
      ENDDO

      IF ( updatePreCond ) THEN
C--   Update overlap regions
      CALL EXCH_XY_RS(aC2d, myThid)

C--   Initialise preconditioner
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1,sNy+1
         DO i=1,sNx+1
          IF ( aC2d(i,j,bi,bj) .EQ. 0. ) THEN
            pC(i,j,bi,bj) = 1. _d 0
          ELSE
           pC(i,j,bi,bj) =  1. _d 0 / aC2d(i,j,bi,bj)
          ENDIF
          pW_tmp = aC2d(i,j,bi,bj)+aC2d(i-1,j,bi,bj)
          IF ( pW_tmp .EQ. 0. ) THEN
           pW(i,j,bi,bj) = 0.
          ELSE
           pW(i,j,bi,bj) =
     &     -aW2d(i,j,bi,bj)/((cg2dpcOffDFac *pW_tmp)**2 )
          ENDIF
          pS_tmp = aC2d(i,j,bi,bj)+aC2d(i,j-1,bi,bj)
          IF ( pS_tmp .EQ. 0. ) THEN
           pS(i,j,bi,bj) = 0.
          ELSE
           pS(i,j,bi,bj) =
     &     -aS2d(i,j,bi,bj)/((cg2dpcOffDFac *pS_tmp)**2 )
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO
C-    if update Preconditioner : end
      ENDIF


      RETURN
      END
