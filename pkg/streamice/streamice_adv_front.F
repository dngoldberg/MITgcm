#include "STREAMICE_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

CBOP
      SUBROUTINE STREAMICE_ADV_FRONT (
     & myThid,
     & time_step,
     & hflux_x_si_v,
     & hflux_y_si_v,
     & hflux_x_si_c,
     & hflux_y_si_c)

C     /============================================================\
C     | SUBROUTINE                                                 |
C     | o                                                          |
C     |============================================================|
C     |                                                            |
C     \============================================================/
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

      INTEGER myThid
      _RL time_step
      _RL hflux_x_SI_v (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hflux_y_SI_v (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hflux_x_SI_c (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hflux_y_SI_c (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

#ifdef ALLOW_STREAMICE

      INTEGER i, j, bi, bj, k, iter_count, iter_rpt
      INTEGER Gi, Gj
      INTEGER new_partial(4), face_neg_flux(4), face_pos_flux(4)

      _RL iter_flag
      _RL n_flux_1, n_flux_2, n_neg_flux, n_pos_flux
      _RL href, rho, partial_vol, tot_flux, hpot
      _RL face_flux, cell_vol, vol_budget
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      _RL hflux_x_SI (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hflux_y_SI (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hflux_x_SI2 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hflux_y_SI2 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hmask_new (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      
      CALL DIAGNOSTICS_FILL(hflux_x_SI_v,'VFLUX_X ',
     &      0,1,0,1,1,myThid)
      CALL DIAGNOSTICS_FILL(hflux_y_SI_v,'VFLUX_Y ',
     &      0,1,0,1,1,myThid)
      CALL DIAGNOSTICS_FILL(hflux_x_SI_c,'CFLUX_X ',
     &      0,1,0,1,1,myThid)
      CALL DIAGNOSTICS_FILL(hflux_y_SI_c,'CFLUX_Y ',
     &      0,1,0,1,1,myThid)

      rho = streamice_density

      iter_flag = 1. _d 0
      iter_rpt = 0

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
            hflux_x_SI2(i,j,bi,bj) = 0. _d 0
            hflux_y_SI2(i,j,bi,bj) = 0. _d 0
            hflux_x_SI(i,j,bi,bj) = hflux_x_SI_v(i,j,bi,bj) + 
     &        hflux_x_SI_c(i,j,bi,bj)
            hflux_y_SI(i,j,bi,bj) = hflux_y_SI_v(i,j,bi,bj) + 
     &        hflux_y_SI_c(i,j,bi,bj)
            hmask_new(i,j,bi,bj) = STREAMICE_hmask(i,j,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      
      DO iter_count = 0, 1


      IF ( iter_flag .GT. 0. ) THEN

       iter_flag = 0. _d 0

       IF (iter_count .gt. 0) then
        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            hflux_x_SI(i,j,bi,bj)=hflux_x_SI2(i,j,bi,bj)
            hflux_y_SI(i,j,bi,bj)=hflux_y_SI2(i,j,bi,bj)
            hflux_x_SI2(i,j,bi,bj) = 0. _d 0
            hflux_y_SI2(i,j,bi,bj) = 0. _d 0
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDIF

!       iter_count = iter_count + 1
       iter_rpt = iter_rpt + 1

       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)

         DO j=1-1,sNy+1
          Gj = (myYGlobalLo-1)+(bj-1)*sNy+j
          IF ((Gj .ge. 1) .and. (Gj .le. Ny)) THEN
           DO i=1-1,sNx+1
            Gi = (myXGlobalLo-1)+(bi-1)*sNx+i

            IF (.not. STREAMICE_calve_to_mask .OR.
     &       STREAMICE_calve_mask (i,j,bi,bj) .eq. 1.0) THEN

             IF ((Gi .ge. 1) .and. (Gi .le. Nx) .and.
     &          (STREAMICE_Hmask(i,j,bi,bj).eq.0.0 .or.
     &           STREAMICE_Hmask(i,j,bi,bj).eq.2.0)) THEN
     
              n_flux_1 = 0. _d 0
              n_pos_flux = 0. _d 0
              n_neg_flux = 0. _d 0
              href = 0. _d 0
              tot_flux = 0. _d 0
              DO k=1,4
                face_pos_flux(k) = 0
                face_neg_flux(k) = 0
              ENDDO

              DO k=0,1
               face_flux = hflux_x_SI(i+k,j,bi,bj) * (1.0-2.0*k)
               IF (face_flux.ne.0. _d 0) THEN
                n_flux_1 = n_flux_1 + 1. _d 0
                tot_flux = tot_flux + face_flux * 
     &           dyG(i+k,j,bi,bj) * time_step
                IF (face_flux.gt.0. _d 0) THEN
                 href = href + H_streamice(i+2*k-1,j,bi,bj)
                 n_pos_flux = n_pos_flux + 1. _d 0
                 face_pos_flux(k+1) = 1
                ELSE
                 n_neg_flux = n_neg_flux + 1. _d 0
                 face_neg_flux(k+1) = 1                 
                ENDIF
                hflux_x_SI(i+k,j,bi,bj) = 0. _d 0
               ENDIF
               
               face_flux = hflux_y_SI(i,j+k,bi,bj) * (1.0-2.0*k)
               IF (face_flux.ne.0. _d 0) THEN
                n_flux_1 = n_flux_1 + 1. _d 0
                tot_flux = tot_flux + face_flux * 
     &           dxG(i,j+k,bi,bj) * time_step
                IF (face_flux.gt.0. _d 0) THEN
                 href = href + H_streamice(i,j+2*k-1,bi,bj)
                 n_pos_flux = n_pos_flux + 1. _d 0
                 face_pos_flux(k+3) = 1
                ELSE
                 n_neg_flux = n_neg_flux + 1. _d 0
                 face_neg_flux(k+3) = 1
                ENDIF
                hflux_y_SI(i,j+k,bi,bj) = 0. _d 0
               ENDIF
              ENDDO 

              IF (n_flux_1 .gt. 0.) THEN

               IF (n_pos_flux .gt. 0.0) THEN
                 href = href / n_pos_flux
               ELSE
                 href = H_streamice(i,j,bi,bj)
               ENDIF
               
               partial_vol = H_streamice (i,j,bi,bj) *
     &          area_shelf_streamice (i,j,bi,bj) + tot_flux
               hpot = partial_vol * recip_rA(i,j,bi,bj)
               
               IF (hpot .eq. 0. _d 0) THEN
               
                hmask_new (i,j,bi,bj) = 0.0
                H_streamice (i,j,bi,bj) = 0.0
                area_shelf_streamice(i,j,bi,bj) = 0.0

               ELSEIF (hpot .lt. 0.0) THEN ! cell has emptied -- contributing
                                      ! cells must become partial
                hmask_new (i,j,bi,bj) = 0.0
                H_streamice (i,j,bi,bj) = 0.0
                area_shelf_streamice(i,j,bi,bj) = 0.0
                
                PRINT *, "GOT HERE UNDERFLOW ", i,j, tot_flux
                 
                n_flux_2 = 0. _d 0 ;
                DO k=1,4
                 new_partial (k) = 0
                ENDDO

                DO k=0,1
                 IF ( STREAMICE_hmask(i+2*k-1,j,bi,bj).eq.1.0
     &              ) THEN  ! at a permanent calving boundary - no advance allowed
                    n_flux_2 = n_flux_2 + 1. _d 0
                    new_partial (k+1) = 1
                 ENDIF
                 IF ( STREAMICE_hmask(i,j+2*k-1,bi,bj).eq.1.0
     &              ) THEN  ! at a permanent calving boundary - no advance allowed
                    n_flux_2 = n_flux_2 + 1. _d 0
                    new_partial (k+3) = 1
                 ENDIF 
                ENDDO

                IF (n_flux_2 .gt. 0.) THEN ! there is nowhere to put the extra ice!
                 H_streamice(i,j,bi,bj) = href

                 DO k=0,1
                 
                  IF (new_partial(k+1) .eq. 1) THEN
                   
                   cell_vol = area_shelf_streamice(i+2*k-1,j,bi,bj) * 
     &              H_streamice(i+2*k-1,j,bi,bj)
                   vol_budget = 1.0 / n_flux_2 * partial_vol
                  
                   IF (cell_vol + vol_budget .lt. 0.0) THEN

                     WRITE(msgBuf,'(A,E15.7)') 'RETREAT PAST CELL',
     &                       cell_vol + vol_budget
                     CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , 1)

                     hmask_new(i+2*k-1,j,bi,bj) = 0.0
                     area_shelf_streamice(i+2*k-1,j,bi,bj) = 0.0
                     H_streamice(i+2*k-1,j,bi,bj) = 0.0
  
                   ELSE

                     cell_vol = cell_vol + vol_budget
                     hmask_new(i+2*k-1,j,bi,bj) = 2.0
                     area_shelf_streamice(i+2*k-1,j,bi,bj) = 
     &                   cell_vol / H_streamice(i+2*k-1,j,bi,bj)

                   ENDIF

                  ENDIF
                  
                  IF (new_partial(k+3) .eq. 1) THEN
                   
                   cell_vol = area_shelf_streamice(i,j+2*k-1,bi,bj) * 
     &              H_streamice(i,j+2*k-1,bi,bj)
                   vol_budget = 1.0 / n_flux_2 * partial_vol
                  
                   IF (cell_vol + vol_budget .lt. 0.0) THEN

                     WRITE(msgBuf,'(A,E15.7)') 'RETREAT PAST CELL',
     &                       cell_vol + vol_budget
                     CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , 1)

                     hmask_new(i,j+2*k-1,bi,bj) = 0.0
                     area_shelf_streamice(i,j+2*k-1,bi,bj) = 0.0
                     H_streamice(i,j+2*k-1,bi,bj) = 0.0
  
                   ELSE

                     cell_vol = cell_vol + vol_budget
                     hmask_new(i,j+2*k-1,bi,bj) = 2.0
                     area_shelf_streamice(i,j+2*k-1,bi,bj) = 
     &                   cell_vol / H_streamice(i,j+2*k-1,bi,bj)

                   ENDIF

                  ENDIF
                  
                 ENDDO

                ENDIF
                 
               ELSEIF (hpot .eq. href) THEN ! cell is exactly covered, no overflow
               
                hmask_new (i,j,bi,bj) = 1.0
                H_streamice (i,j,bi,bj) = href
                area_shelf_streamice(i,j,bi,bj) =
     &           rA(i,j,bi,bj)
     
               ELSEIF (hpot .lt. href) THEN ! cell still unfilled


                hmask_new (i,j,bi,bj) = 2.0
                area_shelf_streamice (i,j,bi,bj) = partial_vol / href
                IF (tot_flux .gt. 0. _d 0) THEN
                 H_streamice (i,j,bi,bj) = href
                ENDIF
                
               ELSE ! cell is filled - do overflow

                PRINT *, "GOT HERE OVERFLOW ", i,j, hpot, href,
     &           area_shelf_streamice(i,j,bi,bj), h_streamice(i,j,bi,bj)
                iter_flag = 1. _d 0
                partial_vol = partial_vol - href * rA(i,j,bi,bj)

                hmask_new (i,j,bi,bj) = 1.0
                area_shelf_streamice(i,j,bi,bj) =
     &           rA(i,j,bi,bj)


                iter_flag  = 1. _d 0

                n_flux_2 = 0. _d 0 ;
                DO k=1,4
                 new_partial (k) = 0
                ENDDO

                DO k=0,1
                 IF ( (STREAMICE_ufacemask(i+k,j,bi,bj).eq.2.0) .or.
     &             (STREAMICE_calve_to_mask .and.
     &              STREAMICE_calve_mask(i+2*k-1,j,bi,bj).ne.1.0)
     &              ) THEN  ! at a permanent calving boundary - no advance allowed
                    n_flux_2 = n_flux_2 + 1. _d 0
                 ELSEIF (STREAMICE_hmask(i+2*k-1,j,bi,bj).eq.0 _d 0) THEN ! adjacent cell is completely ice free
                    n_flux_2 = n_flux_2 + 1. _d 0
                    new_partial (k+1) = 1
                 ENDIF

                 IF ( (STREAMICE_vfacemask (i,j+k,bi,bj).eq.2.0) .or.
     &             (STREAMICE_calve_to_mask .and.
     &              STREAMICE_calve_mask(i,j+2*k-1,bi,bj).ne.1.0)
     &              ) THEN  ! at a permanent calving boundary - no advance allowed
                     n_flux_2 = n_flux_2 + 1. _d 0
                 ELSEIF (STREAMICE_hmask(i,j+2*k-1,bi,bj).eq.0 _d 0) THEN
                     n_flux_2 = n_flux_2 + 1. _d 0
                     new_partial (k+3) = 1
                 ENDIF
                ENDDO

                IF (n_flux_2 .eq. 0.) THEN ! there is nowhere to put the extra ice!
                 H_streamice(i,j,bi,bj) = href + partial_vol *
     &             recip_rA(i,j,bi,bj)
                ELSE
                 H_streamice(i,j,bi,bj) = href

                 DO k=0,1
                  IF (new_partial(k+1) .eq. 1) THEN
                   hflux_x_SI2(i+k,j,bi,bj) =
     &              partial_vol/time_step/n_flux_2/
     &                dyG(i+k,j,bi,bj)
                  ENDIF

                  IF (new_partial(k+3) .eq. 1) THEN
                   hflux_y_SI2(i,j+k,bi,bj) =
     &              partial_vol/time_step/n_flux_2/
     &                dxG(i,j+k,bi,bj)
                  ENDIF
                 ENDDO

                ENDIF
               ENDIF
              ENDIF

             ENDIF
            ENDIF
           ENDDO
          ENDIF
         ENDDO
c
        ENDDO
       ENDDO
c
      ENDIF
      ENDDO
      
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-1,sNy+1
         DO i=1-1,sNx+1
          STREAMICE_hmask(i,j,bi,bj) = hmask_new(i,j,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      IF (iter_rpt.gt.1) THEN
       WRITE(msgBuf,'(A,I5,A)') 'FRONT ADVANCE: ',iter_rpt,
     &  ' ITERATIONS'
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT , 1)
      ENDIF

#endif
      RETURN
      END


