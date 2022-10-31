#include "STREAMICE_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

CBOP
      SUBROUTINE STREAMICE_ADVECT_CREVASSE ( myThid,myIter,time_step )

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

#if (defined(ALLOW_STREAMICE) && defined(STREAMICE_ALLOW_DAMAGE))
      INTEGER i, j, bi, bj, Gi, Gj, k
      _RL thick_bd, uflux, vflux, irho, max_icfl, loc_icfl
      _RL time_step_full, time_step_rem
      _RL sec_per_year, time_step_loc, MR, SMB
      _RL BCVALX(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL BCVALY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS BCMASKX(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS BCMASKY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL utrans(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vtrans(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL d_prev(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL d_after_uflux_SI(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL d_after_vflux_SI(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL dflux_x_SI(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL dflux_y_SI(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL d_bdry_flux_x(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL d_bdry_flux_y(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      CHARACTER*(MAX_LEN_MBUF) msgBuf

      CALL TIMER_START ('STREAMICE_ADVECT_CREVASSE',myThid)

      sec_per_year = 365.*86400.

      time_step_loc = time_step / sec_per_year
      time_step_full = time_step_loc
      time_step_rem = time_step_loc

      if (useStreamiceDamage) then

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE streamice_hmask  = comlev1, key=ikey_dynamics
#endif

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1,sNy+1
         DO i=1,sNx+1

          d_prev(i,j,bi,bj) =
     &     crevasse_depth_adv(i,j,bi,bj)

          d_bdry_flux_x(i,j,bi,bj) = 0. _d 0
          d_bdry_flux_y(i,j,bi,bj) = 0. _d 0

          dflux_x_SI (i,j,bi,bj) = 0. _d 0
          dflux_y_SI (i,j,bi,bj) = 0. _d 0

          IF (STREAMICE_hmask(i,j,bi,bj).eq.1.0) THEN
           d_after_uflux_SI(i,j,bi,bj) = 0. _d 0
           d_after_vflux_SI(i,j,bi,bj) = 0. _d 0
          ELSE
           d_after_uflux_SI(i,j,bi,bj) = d_prev(i,j,bi,bj)
           d_after_vflux_SI(i,j,bi,bj) = d_prev(i,j,bi,bj)
          ENDIF

         ENDDO
        ENDDO
       ENDDO
      ENDDO

      CALL STREAMICE_ADVECT_SETUP (
     &           myThid,
     &           myIter,
     &           d_bdry_flux_x,
     &           d_bdry_flux_y,
     &           utrans,
     &           vtrans,
     &           BCMASKX,
     &           BCMASKY,
     &           BCVALX,
     &           BCVALY )

      _EXCH_XY_RL(d_after_uflux_SI,myThid)
      _EXCH_XY_RL(d_after_vflux_SI,myThid)
      _EXCH_XY_RL(dflux_x_SI,myThid)
      _EXCH_XY_RL(dflux_y_SI,myThid)

#ifndef ALLOW_AUTODIFF

      max_icfl = 1.e-20

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1,sNy
         DO i=1,sNx
          IF (streamice_hmask(i,j,bi,bj).eq.1.0) THEN
           loc_icfl=max(abs(utrans(i,j,bi,bj)),
     &                  abs(utrans(i+1,j,bi,bj))) / dxF(i,j,bi,bj)
           loc_icfl=max(loc_icfl,max(abs(vtrans(i,j,bi,bj)),
     &                  abs(vtrans(i,j+1,bi,bj))) / dyF(i,j,bi,bj))
           if (loc_icfl.gt.max_icfl) then
            max_icfl = loc_icfl
           ENDIF
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      CALL GLOBAL_MAX_R8 (max_icfl, myThid)

#endif /* ALLOW_AUTODIFF */

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE streamice_hmask  = comlev1, key=ikey_dynamics
CADJ STORE H_streamice  = comlev1, key=ikey_dynamics
#endif

#ifndef ALLOW_AUTODIFF
      do while (time_step_rem .gt. 1.e-15)
       time_step_loc = min (
     &  streamice_cfl_factor / max_icfl,
     &  time_step_rem )
       if (time_step_loc .lt. time_step_full) then
        PRINT *, "TAKING PARTIAL TIME STEP", time_step_loc
       endif
#endif /* ALLOW_AUTODIFF */

      CALL STREAMICE_ADV_FLUX_FL_X ( myThid ,
     I   utrans,
     I   crevasse_depth_adv,
     I   BCMASKX,
     I   BCVALX,
     O   dflux_x_SI,
     I   time_step_loc )

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-3,sNy+3
         DO i=1,sNx
          Gi = (myXGlobalLo-1)+(bi-1)*sNx+i
          Gj = (myYGlobalLo-1)+(bj-1)*sNy+j
          IF (((Gj .ge. 1) .and. (Gj .le. Ny))
     &       .or.STREAMICE_NS_PERIODIC) THEN

          IF (STREAMICE_hmask(i,j,bi,bj).eq.1.0) THEN
           d_after_uflux_SI (i,j,bi,bj) = 
     &      crevasse_depth_adv(i,j,bi,bj) -
     &      (dflux_x_SI(i+1,j,bi,bj)*dyG(i+1,j,bi,bj) -
     &        dflux_x_SI(i,j,bi,bj)*dyG(i,j,bi,bj))
     &        * recip_rA (i,j,bi,bj) * time_step_loc
           IF ( d_after_uflux_SI (i,j,bi,bj).le.0.0) THEN
            PRINT *, "h neg after x", i,j,dflux_x_SI(i+1,j,bi,bj),
     &        dflux_x_SI(i,j,bi,bj)
           ENDIF
          ENDIF

          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE streamice_hmask  = comlev1, key=ikey_dynamics
#endif


      CALL STREAMICE_ADV_FLUX_FL_Y ( myThid ,
     I   vtrans ,
     I   d_after_uflux_si ,
     I   BCMASKY,
     I   BCVALY,
     O   dflux_y_SI,
     I   time_step_loc )

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1,sNy
         DO i=1,sNx
          Gi = (myXGlobalLo-1)+(bi-1)*sNx+i
          Gj = (myYGlobalLo-1)+(bj-1)*sNy+j
          IF (((Gj .ge. 1) .and. (Gj .le. Ny))
     &       .or.STREAMICE_EW_PERIODIC) THEN

          IF (STREAMICE_hmask(i,j,bi,bj).eq.1.0) THEN
           d_after_vflux_SI (i,j,bi,bj) = d_after_uflux_SI(i,j,bi,bj) -
     &      (dflux_y_SI(i,j+1,bi,bj)*dxG(i,j+1,bi,bj) -
     &        dflux_y_SI(i,j,bi,bj)*dxG(i,j,bi,bj)) *
     &        recip_rA (i,j,bi,bj) * time_step_loc
           IF ( d_after_vflux_SI (i,j,bi,bj).le.0.0) THEN
            PRINT *, "h neg after y", i,j,dflux_y_SI(i,j+1,bi,bj),
     &        dflux_y_SI(i,j,bi,bj)
           ENDIF

          ENDIF
         ENDIF

         ENDDO
        ENDDO
       ENDDO
      ENDDO

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1,sNy
         DO i=1,sNx
          IF (STREAMICE_hmask(i,j,bi,bj).eq.1.0) THEN
            crevasse_depth_adv(i,j,bi,bj) =
     &      d_after_vflux_SI (i,j,bi,bj)
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO


#ifndef ALLOW_AUTODIFF
      time_step_rem = time_step_rem - time_step_loc
      enddo
#endif /* ALLOW_AUTODIFF */

! NOW WE APPLY MELT RATES !!
! THIS MAY BE MOVED TO A SEPARATE SUBROUTINE
      irho = 1./streamice_density

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1,sNy
         DO i=1,sNx
          Gi = (myXGlobalLo-1)+(bi-1)*sNx+i
          Gj = (myYGlobalLo-1)+(bj-1)*sNy+j
          IF (STREAMICE_hmask(i,j,bi,bj).eq.1.0 .or.
     &       STREAMICE_hmask(i,j,bi,bj).eq.2.0) THEN

           IF (STREAMICE_allow_cpl) THEN
#ifdef ALLOW_SHELFICE
            MR = -1. * (1.-float_frac_streamice(i,j,bi,bj)) *
     &        shelfIceFreshWaterFlux(i,j,bi,bj) * irho *
     &        sec_per_year

#else
            STOP 'SHELFICE IS NOT ENABLED'
#endif
           ELSE
              MR = (1.-float_frac_streamice(i,j,bi,bj)) *
     &             (BDOT_STREAMICE(i,j,bi,bj) +
     &              BDOT_pert(i,j,bi,bj))
           ENDIF

           SMB = ADOT_STREAMICE(i,j,bi,bj)

           crevasse_depth_adv(i,j,bi,bj) = 
     &      crevasse_depth_adv(i,j,bi,bj) - 
     &      (max(SMB,0.0) + max(MR,0.0)) * 
     &      d_prev(i,j,bi,bj)/H_streamice(i,j,bi,bj)      

          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      _EXCH_XY_RL(crevasse_depth_adv,myThid)

      WRITE(msgBuf,'(A)') 'END STREAMICE_ADVECT_CREVASSE'
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT , 1)
      CALL TIMER_STOP ('STREAMICE_ADVECT_CREVASSE',myThid)

      endif

#endif
      RETURN
      END
