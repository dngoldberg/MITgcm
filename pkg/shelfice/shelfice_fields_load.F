#include "SHELFICE_OPTIONS.h"

CStartOfInterface
      SUBROUTINE SHELFICE_FIELDS_LOAD(
     I                            myTime, myIter, myThid )
C     *==========================================================*
C     | SUBROUTINE SHELFICE_FIELDS_LOAD
C     | o Control reading of AIM fields from external source.
C     *==========================================================*
C     | Loads surface condition datasets for SHELFICE.
C     | The routine is called every timetep and periodically
C     | loads a set of external fields.
C     | Files read are
C     | a) shelfice_mass_timedep
C     |  Varies linearly through interval similarly to other 
C     |  forcing fields (e.g. those for AIM)
C     | b) shelfice_tendency_timedep
C     |  The non-thermodynamic mass tendency term. Constant
C     |  over a forcing interval.
C     | If shelfice_mass_timedep is specified, 
C     | shelfice_tendency_timedep is not used
C     *==========================================================*
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SHELFICE.h"



C     === Routine arguments ===
C     myTime :: Simulation time
C     myIter :: Simulation timestep number
C     myThid :: Thread no. that called this routine.
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEndOfInterface

C     === Functions ===

#ifdef ALLOW_SHELFICE
#ifdef ALLOW_SHELFICE_TIMEDEP_FORCING
C     === Local variables ===
C     bi,bj, i,j  :: Loop counters
C     tYear       :: Fraction within year of myTime
C     mnthIndex   :: Current time in whole months
C     prevMnthIndex
C     fNam        :: Strings used in constructing file names
C     mnthNam
C     pfact       :: used to convert Pot.Temp. to in-situ Temp.
C     loadNewData :: true when need to load new data from file
      INTEGER bi,bj,i,j
c     _RL pfact
      LOGICAL first, changed
C--   for use with useMMsurfFc:
      CHARACTER*(MAX_LEN_FNAM) fNam

      INTEGER nm0, nm1, nmP
      _RL myRelTime, fac, tmpFac
      
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C--   find which month to use for surface BC
C     aim_surfForc_TimePeriod :: Length of forcing time period (e.g. 1 month)
C     aim_surfForc_NppCycle   :: Number of time period per Cycle (e.g. 12)

      myRelTime = myTime - startTime
      first = (myRelTime .lt. 0.5*deltaTClock)
      if ( shelfice_forcing_period .eq. 0.D0 ) THEN
!     &     .or. externForcingCycle .eq. 0.D0 ) then
C     control parameter is constant in time and only needs to be updated
C     once in the beginning
       changed = .false.
       nm0  = 1
       nm1  = 1
       fac     = 1.D0
      else

C--   Now calculate whether it is time to update the forcing arrays
       if (externForcingCycle .eq. 0.0 ) THEN
        CALL GET_PERIODIC_INTERVAL(
     O                   nmP, nm0, nm1, tmpFac, fac,
     I                   externForcingCycle, shelfice_forcing_period,
     I                   deltaTclock, 
     I                   myTime+0.5*shelfice_forcing_period, 
     I                   myThid )
        fac = 1.D0 - fac
       else
        CALL GET_PERIODIC_INTERVAL(
     O                   nmP, nm0, nm1, tmpFac, fac,
     I                   externForcingCycle, shelfice_forcing_period,
     I                   deltaTclock, myTime, 
     I                   myThid )
       endif

       IF ( nm0.NE.nmP ) THEN
        changed = .true.
       ELSE
        changed = .false.
       ENDIF
       IF ( first ) changed = .false.
      endif

!      print *, "GOT HERE FIELDS LOAD", fac, first, changed,nm0
 
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C-     Load new data:

         
C-    Only one thread updates parameter in common block
C-    Wait for everyone to set loadNewData before Master updates prevMnthIndex
        _BARRIER

        IF ( first ) THEN

         IF ( SHELFICEMassTimeDepFile .NE. ' '  ) THEN
          CALL READ_REC_XY_RL(SHELFICEMassTimeDepFile,
     &                        mass_shelfice1, 
     &                        nm0,myIter,myThid)
         ENDIF

C-      endif 1rst iter.
        ENDIF

        IF ( first .OR. changed) THEN

         DO bj=myByLo(myThid),myByHi(myThid)
          DO bi=myBxLo(myThid),myBxHi(myThid)
           DO j=1,sNy
            DO i=1,sNx
             mass_shelfice0(i,j,bi,bj) = 
     &        mass_shelfice1(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO

         IF ( SHELFICEMassTimeDepFile .NE. ' '  ) THEN
          CALL READ_REC_XY_RL(SHELFICEMassTimeDepFile,
     &                        mass_shelfice1, 
     &                        nm1,myIter,myThid)
         ELSEIF ( SHELFICEMassTendTimeDepFile .NE. ' '  ) THEN
          CALL READ_REC_XY_RL(SHELFICEMassTendTimeDepFile,
     &                        shelfIceMassDynTendency, 
     &                        nm0,myIter,myThid)
         ENDIF

C-      endif 1rst iter.
        ENDIF


        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
          DO j=1,sNy
           DO i=1,sNx
            IF ( SHELFICEMassTimeDepFile .NE. ' '  ) THEN
             shelficeMass (i,j,bi,bj) = 
     &        fac    * mass_shelfice0(i,j,bi,bj)
     &     + (1-fac) * mass_shelfice1(i,j,bi,bj)
              shelfIceMassDynTendency(i,j,bi,bj) = 0. _d 0
            ENDIF
           ENDDO
          ENDDO
         ENDDO
        ENDDO


      CALL EXCH_XY_RL (shelficemass, myThid)

#endif
#endif

      RETURN
      END
