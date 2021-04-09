#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
      SUBROUTINE PRESSURE_RELEASE_SALT( 
     U                       gS_arr,
     I                       iMin,iMax,jMin,jMax, k, bi,bj,
     I                       myTime, myIter, myThid )
C     *============================================================*
C     | SUBROUTINE PRESSURE_RELEASE_SALT
C     | o Transport salt with darcy flux
C     *============================================================*
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "SURFACE.h"
#include "FFIELDS.h"

C     === Routine arguments ===
C     myThid -  Number of this instance 

      _RL     gS_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER k,bi,bj
      INTEGER iMin,iMax,jMin,jMax
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

CEndOfInterface

#ifdef ALLOW_PRESSURE_RELEASE_CODE

C     === Local Variables ===
      INTEGER i,j,k_e,k_ce,k_s,k_cs,k_w,k_cw,k_n,k_cn
      _RL     S_trans_west,S_trans_east,S_trans_south,S_trans_north
      _RL     vol_trans_west,vol_trans_east,
     &         vol_trans_south,vol_trans_north

      DO j=jMin+1,jMax-1
        DO i=iMin+1,iMax-1

         S_trans_west = 0.0 
         S_trans_north = 0.0
         S_trans_east = 0.0 
         S_trans_south = 0.0
         vol_trans_west = 0.0 
         vol_trans_north = 0.0
         vol_trans_east = 0.0 
         vol_trans_south = 0.0

C calculate the k cells the tracers are transferred between in north,
C south east and west cells.
C Need to find if adjacent cells are deeper or shallower
          IF (kLowC(i,j,bi,bj) .GT. ksurfC(i+1,j,bi,bj)) THEN
             k_e = kLowC(i+1,j,bi,bj)
             k_ce = kSurfC(i,j,bi,bj)
          ELSE
             k_e = kSurfC(i+1,j,bi,bj)
             k_ce = kLowC(i,j,bi,bj)
          ENDIF

          IF (kLowC(i,j,bi,bj) .GT. ksurfC(i-1,j,bi,bj)) THEN
             k_w = kLowC(i-1,j,bi,bj)
             k_cw = kSurfC(i,j,bi,bj)
          ELSE
             k_w = kSurfC(i-1,j,bi,bj)
             k_cw = kLowC(i,j,bi,bj)
          ENDIF

          IF (kLowC(i,j,bi,bj) .GT. ksurfC(i,j+1,bi,bj)) THEN
             k_n = kLowC(i,j+1,bi,bj)
             k_cn = kSurfC(i,j,bi,bj)
          ELSE
             k_n = kSurfC(i,j+1,bi,bj)
             k_cn = kLowC(i,j,bi,bj)
          ENDIF

          IF (kLowC(i,j,bi,bj) .GT. ksurfC(i,j-1,bi,bj)) THEN
             k_s = kLowC(i,j-1,bi,bj)
             k_cs = kSurfC(i,j,bi,bj)
          ELSE
             k_s = kSurfC(i,j-1,bi,bj)
             k_cs = kLowC(i,j,bi,bj)
          ENDIF

C calculate the net tracer flux through north, south east and west
C faces.


          IF (k .EQ. k_cw) THEN
           if (k_cw.gt.0 .and. k_w.gt.0) then
            S_trans_west =0.5 _d 0  * 
     &                  ( pReleaseTransX(i,j,bi,bj) *
     &                  (salt(i-1,j,k_w,bi,bj)+salt(i,j,k_cw,bi,bj))
     &                   +abs(pReleaseTransX(i,j,bi,bj)) * 
     &                  (salt(i-1,j,k_w,bi,bj)-salt(i,j,k_cw,bi,bj)))
     &                  *recip_rA(i,j,bi,bj)*_dyG(i,j,bi,bj)
     &                  *recip_deepFac2C(k_cw)*recip_rhoFacC(k_cw)
     &                  *recip_drF(k_cw)*_recip_hFacC(i,j,k_cw,bi,bj)
            vol_trans_west = pReleaseTransX(i,j,bi,bj)
     &                  *recip_rA(i,j,bi,bj)*_dyG(i,j,bi,bj)
     &                  *recip_deepFac2C(k_cw)*recip_rhoFacC(k_cw)
     &                  *recip_drF(k_cw)*_recip_hFacC(i,j,k_cw,bi,bj)
           endif
          ENDIF
          IF (k .EQ. k_ce) THEN
           if (k_ce.gt.0 .and. k_e.gt.0) then
            S_trans_east =0.5 _d 0  *
     &                  ( pReleaseTransX(i+1,j,bi,bj) *
     &                  (salt(i,j,k_ce,bi,bj)+salt(i+1,j,k_e,bi,bj))
     &                   +abs(pReleaseTransX(i+1,j,bi,bj)) *
     &                  (salt(i,j,k_ce,bi,bj)-salt(i+1,j,k_e,bi,bj)))
     &                  *recip_rA(i,j,bi,bj)*_dyG(i+1,j,bi,bj)
     &                  *recip_deepFac2C(k_ce)*recip_rhoFacC(k_ce)
     &                  *recip_drF(k_ce)*_recip_hFacC(i,j,k_ce,bi,bj)
            vol_trans_east = pReleaseTransX(i+1,j,bi,bj)
     &                  *recip_rA(i,j,bi,bj)*_dyG(i+1,j,bi,bj)
     &                  *recip_deepFac2C(k_ce)*recip_rhoFacC(k_ce)
     &                  *recip_drF(k_ce)*_recip_hFacC(i,j,k_ce,bi,bj)
!            S_trans_east =pReleaseTransX(i+1,j,bi,bj)*
!     &                  (salt(i,j,k_ce,bi,bj) -salt(i+1,j,k_e,bi,bj))
!     &                  *recip_dxG(i,j,bi,bj)
!     &                  *recip_drF(k_ce)*_recip_hFacC(i,j,k_ce,bi,bj)
           endif
          ENDIF
          IF (k .EQ. k_cs) THEN
           if (k_cs.gt.0 .and. k_s.gt.0) then
            S_trans_south =0.5 _d 0  *
     &                  ( pReleaseTransY(i,j,bi,bj) *
     &                  (salt(i,j-1,k_s,bi,bj)+salt(i,j,k_cs,bi,bj))
     &                   +abs(pReleaseTransY(i,j,bi,bj)) *
     &                  (salt(i,j-1,k_s,bi,bj)-salt(i,j,k_cs,bi,bj)))
     &                  *recip_rA(i,j,bi,bj)*_dxG(i,j,bi,bj)
     &                  *recip_deepFac2C(k_cs)*recip_rhoFacC(k_cs)
     &                  *recip_drF(k_cs)*_recip_hFacC(i,j,k_cs,bi,bj)
            vol_trans_south = pReleaseTransY(i,j,bi,bj)
     &                  *recip_rA(i,j,bi,bj)*_dxG(i,j,bi,bj)
     &                  *recip_deepFac2C(k_cs)*recip_rhoFacC(k_cs)
     &                  *recip_drF(k_cs)*_recip_hFacC(i,j,k_cs,bi,bj)
!            S_trans_south =pReleaseTransY(i,j,bi,bj)*
!     &                  (salt(i,j-1,k_s,bi,bj) -salt(i,j,k_cs,bi,bj))
!     &                  *recip_dyG(i,j,bi,bj)
!     &                  *recip_drF(k_cs)*_recip_hFacC(i,j,k_cs,bi,bj)
           endif
          ENDIF
          IF (k .EQ. k_cn) THEN
           if (k_cn.gt.0 .and. k_n.gt.0) then
            S_trans_north =0.5 _d 0  *
     &                  ( pReleaseTransY(i,j+1,bi,bj) *
     &                  (salt(i,j,k_cn,bi,bj)+salt(i,j+1,k_n,bi,bj))
     &                   +abs(pReleaseTransY(i,j+1,bi,bj)) *
     &                  (salt(i,j,k_cn,bi,bj)-salt(i,j+1,k_n,bi,bj)))
     &                  *recip_rA(i,j,bi,bj)*_dxG(i,j+1,bi,bj)
     &                  *recip_deepFac2C(k_cn)*recip_rhoFacC(k_cn)
     &                  *recip_drF(k_cn)*_recip_hFacC(i,j,k_cn,bi,bj)
            vol_trans_north = pReleaseTransY(i,j+1,bi,bj)
     &                  *recip_rA(i,j,bi,bj)*_dxG(i,j+1,bi,bj)
     &                  *recip_deepFac2C(k_cn)*recip_rhoFacC(k_cn)
     &                  *recip_drF(k_cn)*_recip_hFacC(i,j,k_cn,bi,bj)
!            S_trans_north =pReleaseTransY(i,j+1,bi,bj)*
!     &                  (salt(i,j,k_cn,bi,bj) -salt(i,j+1,k_n,bi,bj))
!     &                  *recip_dyG(i,j,bi,bj)
!     &                  *recip_drF(k_cn)*_recip_hFacC(i,j,k_cn,bi,bj)
           endif
          ENDIF


C Add to get total tracer tendency.
          gS_arr(i,j) = gS_arr(i,j) + S_trans_west - S_trans_east
     &                 + S_trans_south - S_trans_north  
          if (k.eq.k_cn) then
           gS_arr(i,j) = gS_arr(i,j)  
     &         - salt(i,j,k_cn,bi,bj) * (-vol_trans_north)
          endif
          if (k.eq.k_cs) then
           gS_arr(i,j) = gS_arr(i,j)    
     &         - salt(i,j,k_cs,bi,bj) * (vol_trans_south)
          endif
          if (k.eq.k_ce) then
           gS_arr(i,j) = gS_arr(i,j)  
     &         - salt(i,j,k_ce,bi,bj) * (-vol_trans_east)
          endif
          if (k.eq.k_cw) then
           gS_arr(i,j) = gS_arr(i,j)  
     &         - salt(i,j,k_cw,bi,bj) * (vol_trans_west)
          endif
         ENDDO
       ENDDO


#endif /* ALLOW_PRESSURE_RELEASE_CODE */

      RETURN
      END
