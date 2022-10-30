#include "STREAMICE_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

CBOP
      SUBROUTINE STREAMICE_ADVECT_SETUP ( 
     I           myThid,
     I           myIter,
     I           flux_bdry_x,
     I           flux_bdry_y,
     o           utrans,
     o           vtrans,
     o           BCMASKX,
     o           BCMASKY,
     o           BCVALX,
     o           BCVALY )


C     *============================================================*
C     | SUBROUTINE                                                 |
C     | o                                                          |
C     *============================================================*
C     |                                                            |
C     *============================================================*
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "STREAMICE.h"
#include "STREAMICE_ADV.h"
#ifdef ALLOW_AUTODIFF_TAMC
# include "tamc.h"
#endif
#ifdef ALLOW_SHELFICE
# include "SHELFICE.h"
#endif

      INTEGER myThid, myIter
      _RL time_step
      _RL flux_bdry_x(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL flux_bdry_y(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL BCVALY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL BCVALX(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL BCVALY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS BCMASKX(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS BCMASKY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL utrans(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vtrans(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

#ifdef ALLOW_STREAMICE

      INTEGER i, j, bi, bj, Gi, Gj, k

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1,sNy+1
         DO i=1,sNx+1

          IF (STREAMICE_ufacemask(i,j,bi,bj).eq.3.0) THEN
           BCMASKX(i,j,bi,bj) = 3.0
           BCVALX(i,j,bi,bj) = h_ubdry_values_SI(i,j,bi,bj)
           IF (STREAMICE_vel_ext_cgrid) THEN
            utrans(i,j,bi,bj) = u_streamice_ext(i,j,bi,bj)
           ELSEIF (STREAMICE_vel_ext) THEN
            utrans(i,j,bi,bj) = .5 * (
     &       u_streamice_ext(i,j,bi,bj)+u_streamice_ext(i,j+1,bi,bj))
           ELSE
            utrans(i,j,bi,bj) = .5 * (
     &       u_streamice(i,j,bi,bj)+u_streamice(i,j+1,bi,bj))
           ENDIF
          ELSEIF (STREAMICE_ufacemask(i,j,bi,bj).eq.4.0) THEN
           IF (STREAMICE_hmask(i,j,bi,bj).eq.1.0) THEN
            uflux = flux_bdry_x(i,j,bi,bj)
            BCMASKX(i,j,bi,bj) = 3.0
            BCVALX(i,j,bi,bj) = uflux
            utrans(i,j,bi,bj) = 1.0
           ELSEIF (STREAMICE_hmask(i-1,j,bi,bj).eq.1.0) THEN
            uflux =flux_bdry_x(i,j,bi,bj)
            BCMASKX(i,j,bi,bj) = 3.0
            BCVALX(i,j,bi,bj) = uflux
            utrans(i,j,bi,bj) = -1.0
           ENDIF
          ELSEIF (.not.(
     &     STREAMICE_hmask(i,j,bi,bj).eq.1.0.OR.
     &     STREAMICE_hmask(i-1,j,bi,bj).eq.1.0)) THEN
           BCMASKX(i,j,bi,bj) = 0.0
           BCVALX(i,j,bi,bj) = 0. _d 0
           utrans(i,j,bi,bj) = 0. _d 0
          ELSE
           BCMASKX(i,j,bi,bj) = 0.0
           BCVALX(i,j,bi,bj) = 0. _d 0
           IF (STREAMICE_vel_ext_cgrid) THEN
            utrans(i,j,bi,bj) = u_streamice_ext(i,j,bi,bj)
           ELSEIF (STREAMICE_vel_ext) THEN
            utrans(i,j,bi,bj) = .5 * (
     &       u_streamice_ext(i,j,bi,bj)+u_streamice_ext(i,j+1,bi,bj))
           ELSE
            utrans(i,j,bi,bj) = .5 * (
     &       u_streamice(i,j,bi,bj)+u_streamice(i,j+1,bi,bj))
           ENDIF
          ENDIF

          IF (STREAMICE_vfacemask(i,j,bi,bj).eq.3.0) THEN
           BCMASKy(i,j,bi,bj) = 3.0
           BCVALy(i,j,bi,bj) = h_vbdry_values_SI(i,j,bi,bj)
           IF (STREAMICE_vel_ext_cgrid) THEN
            vtrans(i,j,bi,bj) = v_streamice_ext(i,j,bi,bj)
           ELSEIF (STREAMICE_vel_ext) THEN
            vtrans(i,j,bi,bj) = .5 * (
     &       v_streamice_ext(i,j,bi,bj)+v_streamice_ext(i+1,j,bi,bj))
           ELSE
            vtrans(i,j,bi,bj) = .5 * (
     &       v_streamice(i,j,bi,bj)+v_streamice(i+1,j,bi,bj))
           ENDIF
          ELSEIF (STREAMICE_vfacemask(i,j,bi,bj).eq.4.0) THEN
           IF (STREAMICE_hmask(i,j,bi,bj).eq.1.0) THEN
            vflux = flux_bdry_y(i,j,bi,bj)
            BCMASKY(i,j,bi,bj) = 3.0
            BCVALY(i,j,bi,bj) = vflux
            vtrans(i,j,bi,bj) = 1.0
           ELSEIF (STREAMICE_hmask(i,j-1,bi,bj).eq.1.0) THEN
            vflux = flux_bdry_y(i,j,bi,bj)
            BCMASKY(i,j,bi,bj) = 3.0
            BCVALY(i,j,bi,bj) = vflux
            vtrans(i,j,bi,bj) = -1.0
           ENDIF
          ELSEIF (.not.(
     &     STREAMICE_hmask(i,j,bi,bj).eq.1.0.OR.
     &     STREAMICE_hmask(i,j-1,bi,bj).eq.1.0)) THEN
           BCMASKY(i,j,bi,bj) = 0.0
           BCVALY(i,j,bi,bj) = 0. _d 0
           vtrans(i,j,bi,bj) = 0. _d 0
          ELSE
           BCMASKy(i,j,bi,bj) = 0.0
           BCVALy(i,j,bi,bj) = 0. _d 0
           IF (STREAMICE_vel_ext_cgrid) THEN
            vtrans(i,j,bi,bj) = v_streamice_ext(i,j,bi,bj)
           ELSEIF (STREAMICE_vel_ext) THEN
            vtrans(i,j,bi,bj) = .5 * (
     &       v_streamice_ext(i,j,bi,bj)+v_streamice_ext(i+1,j,bi,bj))
           ELSE
            vtrans(i,j,bi,bj) = .5 * (
     &       v_streamice(i,j,bi,bj)+v_streamice(i+1,j,bi,bj))
           ENDIF
          ENDIF

         ENDDO
        ENDDO
       ENDDO
      ENDDO 

      _EXCH_XY_RL(utrans,myThid)
      _EXCH_XY_RL(vtrans,myThid)
      _EXCH_XY_RS(BCMASKx,myThid)
      _EXCH_XY_RS(BCMASKy,myThid)
      _EXCH_XY_RL(BCVALX,myThid)
      _EXCH_XY_RL(BCVALY,myThid)     

#endif
      RETURN
      END 
