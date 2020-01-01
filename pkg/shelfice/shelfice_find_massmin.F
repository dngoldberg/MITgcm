#include "SHELFICE_OPTIONS.h"

CBOP
      SUBROUTINE SHELFICE_FIND_MASSMIN( 
     U                             massMin, 
     I                             myThid )
C     *============================================================*
C     | SUBROUTINE SHELFICE_CALC_GRD_FAC
C     *============================================================*
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SHELFICE.h"
#ifdef ALLOW_STREAMICE
#include "STREAMICE.h"
#endif
#include "DYNVARS.h"
#ifdef ALLOW_COST
# include "SHELFICE_COST.h"
#endif /* ALLOW_COST */

C     === Routine arguments ===
C     myThid -  Number of this instance 
C     
C     This subroutine calculates the minimum mass to allow MWCT space
C     under the grounded ice, either by assuming a uniform density or
C     by calling another S/R which calculates an average density pro-
C     file.

      _RL massMin(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER myThid 

CEndOfInterface

#ifdef ALLOW_SHELFICE
C     === Local variables ===
C     i,j,bi,bj - Loop counters
      INTEGER i, j, k, bi, bj
      _RL SEALEVEL, ETA
      _RL oce_density, mass
      _RL R_min(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CEOP

#ifdef ALLOW_SHELFICE_GROUNDED_ICE

#ifdef ALLOW_STREAMICE
      IF (.not.usestreamice) THEN
#endif
         oce_density = 1028.
#ifdef ALLOW_STREAMICE
      ELSE
         oce_density = streamice_density_ocean_avg
      ENDIF
#endif

      IF (shelfice_massmin_trueDens) then
 
        CALL SHELFICE_MASSMIN(massMin,myThid)
      
      ELSE
  
       SEALEVEL = 0. _d 0
       ETA = 0. _d 0
       CALL SHELFICE_SEA_LEVEL_AVG( SEALEVEL,ETA, myThid )
   
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j = 1-OLy, sNy+OLy
          DO i = 1-OLx, sNx+OLx

           massMin(i,j,bi,bj) = oce_density *
     &        (SEALEVEL-(R_low(i,j,bi,bj)+R_MWCT(i,j,bi,bj)))

          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDIF
      
#endif
#endif /* ALLOW_SHELFICE */

      RETURN
      END
