!
!                    **********************
                     SUBROUTINE INTERPMETEO
!                    **********************
!
     &(WW,WINDX,WINDY,TAIR,PATM,HREL,NEBU,RAINFALL,EVAPORATION,AT,NFO)
!
!***********************************************************************
! TELEMAC2D   V7P2
!***********************************************************************
!
!brief    READS AND INTERPOLATES VARIABLES IN AN ASCII FILE
!+
!
!history  R. SAMIE, E. RAZAFINDRAKOTO, C.-T. PHAM (EDF-LNHE)
!+        09/07/2014
!+        V7P0
!+        FROM LAPLUIE AND LECENT
!+
!
!history  A. LEROY (EDF-LNHE)
!+        25/11/2015
!+        V7P1
!+        The interpolation of the norm of the wind velocity
!+        is now correct and the Secchi length is added in the
!+        variables to be read
!+
!history  R.ATA (EDF-LNHE)
!+        25/03/2016
!+        V7P2
!+        harmonization with meteo file, last column is evaporation
!+         not secchi length
!+        all out are changed to inout
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| AT             |-->| TIME OF TIME STEP
!| HREL           |<--| RELATIVE HUMIDITY
!| NEBU           |<--| NEBULOSITY
!| NFO            |-->| LOGICAL UNIT OF THE FORMATTED DATA FILE
!| PATM           |<--| ATMOSPHERIC PRESSURE
!| RAINFALL       |<--| RAINFALL
!| EVAPORATION    |<--| EVAPORATION RATE
!| TAIR           |<--| AIR TEMPERATURE
!| WINDX          |<--| WIND ALONG X
!| WINDY          |<--| WIND ALONG Y
!| WW             |<--| MAGNITUDE OF WIND VELOCITY
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE DECLARATIONS_TELEMAC2D, ONLY: DEJA_IPM,DUMMY_IPM,TABENT_IPM,
     &                                  POSTAB_IPM,NBENR_IPM
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE
!
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN)             :: NFO
      DOUBLE PRECISION, INTENT(IN)    :: AT
      DOUBLE PRECISION, INTENT(INOUT) :: WW,WINDX,WINDY,TAIR,PATM,HREL
      DOUBLE PRECISION, INTENT(INOUT) :: NEBU,RAINFALL,EVAPORATION
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER I
!     NUMBER OF VARIABLES IN THE ASCII FILE NFO
      INTEGER, PARAMETER :: NINPUTVAR = 9
      INTEGER, PARAMETER :: NLINESTOSKIP = 3
!
      DOUBLE PRECISION DELTAT,ALPHA
      DOUBLE PRECISION DTR
!
!-----------------------------------------------------------------------
!
      DTR = ATAN(1.D0)/45.D0
!
!-----------------------------------------------------------------------
!
!  READS INPUT DATA FILE
!  AT THE FIRST TIME STEP AND FILLS IN TABENT_IPM ARRAY
!
      IF (.NOT.DEJA_IPM) THEN
!
        IF(LNG.EQ.1) THEN
          WRITE(LU,*)'====================================='
          WRITE(LU,*)'DEBUT DE LECTURE DU FICHIER D''ENTREE'
          WRITE(LU,*)'====================================='
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*)'=================================='
          WRITE(LU,*)'BEGINNING OF READING OF INPUT FILE'
          WRITE(LU,*)'=================================='
        ENDIF
!
        REWIND NFO
!
!  READS THE HEADLINE OF THE DATA FILE
!
        DO I=1,NLINESTOSKIP
          READ(NFO,*)
          WRITE(LU,*)
        ENDDO
!
        NBENR_IPM = 0
!
        ALLOCATE(DUMMY_IPM(NINPUTVAR))
!
!  READS VARIABLES AND FILLS IN TABENT_IPM ARRAY
!
 100    READ(NFO,*,END=20) DUMMY_IPM(1:NINPUTVAR)
        NBENR_IPM = NBENR_IPM + 1
        GO TO 100
!
 20     CONTINUE
!
!  +2 TO STORE 2 EXTRA DATA (THE X AND Y COMPONENTS OF WIND VELOCITY)
        ALLOCATE(TABENT_IPM(NBENR_IPM,NINPUTVAR+2))
        DEJA_IPM = .TRUE.
!
!-----------------------------------------------------------------------
!
        REWIND NFO
!
!  READS THE HEADLINE OF THE DATA FILE
!
        DO I=1,NLINESTOSKIP
          READ(NFO,*)
        ENDDO
!
!  READS VARIABLES AND FILLS IN TABENT_IPM ARRAY
!
        DO I=1,NBENR_IPM
          READ(NFO,*)TABENT_IPM(I,1:NINPUTVAR)
        ENDDO
!
        IF(LNG.EQ.1) THEN
          WRITE(LU,*)'======================================='
          WRITE(LU,*)'  FIN DE LECTURE DU FICHIER D''ENTREE  '
          WRITE(LU,*)'     IL Y A ',NBENR_IPM,' ENREGISTREMENTS  '
          WRITE(LU,*)'  DE T = ',TABENT_IPM(1,1), ' A  = ',
     &               TABENT_IPM(NBENR_IPM,1),' SECONDES '
          WRITE(LU,*)'======================================='
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*)'======================================='
          WRITE(LU,*)'  END OF READING OF INPUT DATA         '
          WRITE(LU,*)'  THERE ARE ',NBENR_IPM,' RECORDS          '
          WRITE(LU,*)'  FROM T = ',TABENT_IPM(1,1), ' TO = ',
     &               TABENT_IPM(NBENR_IPM,1),' SECONDS '
          WRITE(LU,*)'======================================='
        ENDIF
!
!  FILLS IN WITH X AND Y COMPONENTS OF WIND VELOCITY
!
        DO I=1,NBENR_IPM
          TABENT_IPM(I,NINPUTVAR+1) =
     &             -TABENT_IPM(I,2)*SIN(TABENT_IPM(I,3)*DTR)
          TABENT_IPM(I,NINPUTVAR+2) =
     &             -TABENT_IPM(I,2)*COS(TABENT_IPM(I,3)*DTR)
        ENDDO
!
!  END TO FILL IN TABENT_IPM ARRAY
!  INITIALISATION OF THE POINTER POSTAB_IPM
!
        POSTAB_IPM = 1
!
      ENDIF
!
!-----------------------------------------------------------------------
!
!  INTERPOLATES DATA AT EACH TIME STEP
!  POINTER POSTAB_IPM POSITIONNED AT EACH TIME STEP
!
120   IF(AT.LT.TABENT_IPM(POSTAB_IPM,1) .OR.
     &   AT.GE.TABENT_IPM(POSTAB_IPM+1,1))
     &  THEN
        IF(AT.LT.TABENT_IPM(POSTAB_IPM,1))   POSTAB_IPM = POSTAB_IPM - 1
        IF(AT.GE.TABENT_IPM(POSTAB_IPM+1,1)) POSTAB_IPM = POSTAB_IPM + 1
        IF(POSTAB_IPM.GT.NBENR_IPM) THEN
          IF(LNG.EQ.1) THEN
            WRITE(LU,*)'==============================================='
            WRITE(LU,*)'ATTENTION : LE TEMPS DU CALCUL AT = ', AT
            WRITE(LU,*)'EST SUPERIEUR AU TEMPS MAXIMUM DE VOTRE FICHIER'
            WRITE(LU,*)'DE DONNEES D''ENTREE T = ',
     &                  TABENT_IPM(NBENR_IPM,1)
            WRITE(LU,*)'==============================================='
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*)'==============================================='
            WRITE(LU,*)'WARNING: TIME OF CALCULATION AT = ', AT
            WRITE(LU,*)'IS BIGGER THAN MAXIMUM TIME IN YOUR INPUT DATA '
            WRITE(LU,*)'FILE T = ', TABENT_IPM(NBENR_IPM,1)
            WRITE(LU,*)'==============================================='
          ENDIF
          CALL PLANTE(1)
        ELSEIF(POSTAB_IPM.LT.1) THEN
          IF(LNG.EQ.1) THEN
            WRITE(LU,*)'==============================================='
            WRITE(LU,*)'ATTENTION : LE TEMPS DU CALCUL AT = ', AT
            WRITE(LU,*)'EST INFERIEUR AU TEMPS MINIMUM DE VOTRE FICHIER'
            WRITE(LU,*)'DE DONNEES D''ENTREE T = ', TABENT_IPM(1,1)
            WRITE(LU,*)'==============================================='
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*)'==============================================='
            WRITE(LU,*)'WARNING: TIME OF CALCULATION AT = ', AT
            WRITE(LU,*)'IS LOWER THAN MINIMUM TIME IN YOUR INPUT DATA '
            WRITE(LU,*)'FILE T = ', TABENT_IPM(1,1)
            WRITE(LU,*)'==============================================='
          ENDIF
          CALL PLANTE(1)
        ENDIF
        GO TO 120
      ENDIF
!
      DELTAT = TABENT_IPM(POSTAB_IPM+1,1)-TABENT_IPM(POSTAB_IPM,1)
      ALPHA  = (AT-TABENT_IPM(POSTAB_IPM,1))/DELTAT
!
!-----------------------------------------------------------------------
!
      WINDX =  TABENT_IPM(POSTAB_IPM,NINPUTVAR+1)
     &      + (TABENT_IPM(POSTAB_IPM+1,NINPUTVAR+1)
     &      -  TABENT_IPM(POSTAB_IPM,NINPUTVAR+1))*ALPHA
      WINDY =  TABENT_IPM(POSTAB_IPM,NINPUTVAR+2)
     &      + (TABENT_IPM(POSTAB_IPM+1,NINPUTVAR+2)
     &      -  TABENT_IPM(POSTAB_IPM,NINPUTVAR+2))*ALPHA
      WW    =  SQRT(WINDX**2+WINDY**2)
!
      TAIR  =  TABENT_IPM(POSTAB_IPM,4)
     &      + (TABENT_IPM(POSTAB_IPM+1,4)
     &         -TABENT_IPM(POSTAB_IPM,4))*ALPHA
      PATM  =  TABENT_IPM(POSTAB_IPM,5)
     &      + (TABENT_IPM(POSTAB_IPM+1,5)
     &         -TABENT_IPM(POSTAB_IPM,5))*ALPHA
      HREL  =  TABENT_IPM(POSTAB_IPM,6)
     &      + (TABENT_IPM(POSTAB_IPM+1,6)
     &         -TABENT_IPM(POSTAB_IPM,6))*ALPHA
      NEBU  =  TABENT_IPM(POSTAB_IPM,7)
     &      + (TABENT_IPM(POSTAB_IPM+1,7)
     &         -TABENT_IPM(POSTAB_IPM,7))*ALPHA
!
      RAINFALL = TABENT_IPM(POSTAB_IPM+1,8)/DELTAT
!
      EVAPORATION   = TABENT_IPM(POSTAB_IPM,9)
     &     + (TABENT_IPM(POSTAB_IPM+1,9)-TABENT_IPM(POSTAB_IPM,9))*ALPHA
!
!-----------------------------------------------------------------------
!
      RETURN
      END

!                    **********************
                     SUBROUTINE SOURCE_TRAC
!                    **********************
     & (LT)
!
!
!***********************************************************************
! TELEMAC3D   V7P1
!***********************************************************************
!
!brief    PREPARES SOURCE TERMS FOR DIFFUSION OF TRACERS.
!
!history  CDG/SOGREAH
!+        **/06/2001
!+
!+   TRACER SOURCES
!
!history  J-M HERVOUET (LNHE)
!+        21/10/2004
!+        V5P5
!+
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        13/07/2010
!+        V6P0
!+   Translation of French comments within the FORTRAN sources into
!+   English comments
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        21/08/2010
!+        V6P0
!+   Creation of DOXYGEN tags for automated documentation and
!+   cross-referencing of the FORTRAN sources
!
!history  A. GINEAU, N. DURAND, N. LORRAIN, C.-T. PHAM (LNHE)
!+        09/07/2014
!+        V7P0
!+   Adding an example of the penetration of the solar radiation
!+   for exchange with atmosphere
!
!history  A. LEROY (LNHE)
!+        25/11/2015
!+        V7P1
!+   Remove the call to INTERPMETEO: all the meteo variables are now
!+   read in meteo.f and stored in variables of waqtel
!+
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| LT             |-->| ITERATION NUMBER
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE DECLARATIONS_TELEMAC,ONLY : COUPLING
      USE DECLARATIONS_TELEMAC3D
!     HEAT EXCHANGE WITH ATMOSPHERE
      USE EXCHANGE_WITH_ATMOSPHERE
      USE DECLARATIONS_WAQTEL,ONLY : NEBU,ZSD,WAQPROCESS,RAYEFF
      USE INTERFACE_WAQTEL
      USE INTERFACE_TELEMAC3D,EX_SOURCE_TRAC =>SOURCE_TRAC
!
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN)           :: LT
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!----------------------------------------------------------------------
!
      INTEGER ITRAC
      DOUBLE PRECISION DUMM(2)
!
      INTEGER IPLAN,I,J,IHOUR
      DOUBLE PRECISION LATITUDE,LONGITUDE,RADSOL,RAY_SOL
      DOUBLE PRECISION KD
!
!----------------------------------------------------------------------
!
!     SETS SOURCE TERMS TO ZERO
!
      IF(NTRAC.GE.1) THEN
!
!       CALL OS ( 'X=C     ' , X=S0TA , C=0.D0 )
!       CALL OS ( 'X=C     ' , X=S1TA , C=0.D0 )
!
!       SOURCE TERMS SIMPLY MARKED
!
!       BEWARE, PUT Q INSTEAD OF 0 IN TYPR IF NOT NIL
!
        LATITUDE  = 50.33
        LONGITUDE = -3.55
!        WRITE(LU,*) 'LAT: ',LATITUDE
        CALL SOLRAD(RAY_SOL,NEBU,MARDAT,MARTIM,AT,LATITUDE,LONGITUDE)
        IHOUR=MARTIM(1)
        ZSD = 0.9D0 ! MAY BE UNCOMMENTED TOO
        KD  = 1.7D0/ZSD ! 83% OF THE INCIDENT ENERGY IS ABSORBED
        DO ITRAC=1,NTRAC
          S0TA%ADR(ITRAC)%P%TYPR='0'
          S1TA%ADR(ITRAC)%P%TYPR='0'
          IF(ITRAC.GE.3) THEN
             S1TA%ADR(3)%P%TYPR='Q'
             DO I = 1,NPOIN2
               DO IPLAN = 1,NPLAN
                 J = I + (IPLAN-1)*NPOIN2
!                 write (*,*) "salinity: ", TA%ADR(IND_S)%P%R(IPOIN3)         
!                 S1TA%ADR(3)%P%R(IPOIN3) =(0.036D0+0.02D0*TA%ADR(IND_S)%P%R(IPOIN3))*
!      &           1.07D0*(TA%ADR(IND_T)%P%R(IPOIN3)-20.D0)
                 RADSOL=RAY_SOL*EXP(KD*(Z(J)-Z(I+(NPLAN-1)*NPOIN2)))
                 IF(IHOUR.GE.6.AND.IHOUR.LE.18) THEN
                    S1TA%ADR(3)%P%R(J) = 1/3600.D0*0.27D0
                 ELSE
                    S1TA%ADR(3)%P%R(J) =1/3600.D0*0.075D0
                 ENDIF
               ENDDO
             ENDDO
          ENDIF
        ENDDO
!
!       EXAMPLE OF RADIOACTIVE DECAY E**(-KT) ON FIRST TRACER, HERE C=K
!
!       S1TA%ADR(1)%P%TYPR='Q'
!       CALL OS('X=C     ',S1TA%ADR(1)%P,C=1.D0)
!
!!       EXAMPLE OF PENETRATION OF THE SOLAR RADIATION
!!       IF EXCHANGE WITH ATMPOSPHERE IS USED, DO NOT FORGET TO
!!       UNCOMMENT THE DECLARATIONS OF THE MODULE + THE VARIABLES ABOVE
!!
!!       SOURCE IN TEMPERATURE NOT EQUAL TO ZERO
!        S0TA%ADR(IND_T)%P%TYPR='Q'
!!       INCIDENT SOLAR RADIATION
!!       LATITUDE AND LONGITUDE TO BE CHANGED DEPENDING ON THE MEAN LOCATION


!        CALL SOLRAD(RAY_SOL,NEBU,MARDAT,MARTIM,AT,LATITUDE,LONGITUDE)
!       FORMULA FOR TURBID WATER WITH SECCHI LENGTH!!       ZSD: SECCHI LENGTH DECLARED IN THE WAQTEL MODULE
!       ZSD: SECCHI LENGTH DECLARED IN THE WAQTEL MODULE
!       IT CAN BE READ IN A FILE LIKE OTHER METEO DATA
!       IF NOT PROVIDED, THE USER CAN TRY TO USE A CONSTANT HERE, E.G.
!        ZSD = 0.9D0 ! MAY BE UNCOMMENTED TOO
!        KD  = 1.7D0/ZSD ! 83% OF THE INCIDENT ENERGY IS ABSORBED
!        DO I=1,NPOIN2
!          DO IPLAN=1,NPLAN
!            J = I + (IPLAN-1)*NPOIN2
!            TREEL=TA%ADR(IND_T)%P%R(NPOIN3-NPOIN2+I)
!            IF (IND_S.NE.0) THEN
!              SAL = TA%ADR(IND_S)%P%R(NPOIN3-NPOIN2+I)
!            ENDIF
!            S0TA%ADR(IND_T)%P%R(J) =
!     &       KD*EXP(KD*(Z(J)-Z(I+(NPLAN-1)*NPOIN2)))
!     &      *RAY_SOL/LAMB
!            RADSOL=RAY_SOL*EXP(KD*(Z(J)-Z(I+(NPLAN-1)*NPOIN2)))
            !WRITE(LU,*) 'RADSOL: ',RADSOL
            !WRITE(LU,*) 'RAY_SOL: ',RAY_SOL
            !WRITE(LU,*) 'DEPTH: ',Z(J)-Z(I+(NPLAN-1)*NPOIN2)
!
!           EXAMPLE OF FORMULA FOR TURBID WATER
!           ALL CONSTANTS MAY BE TUNED
!           0.22D0 = 1.D0-0.78D0
!           S0TA%ADR(IND_T)%P%R(J) =
!           ( 0.78D0*0.66D0 *EXP(0.66D0* (Z(J)-Z(I+(NPLAN-1)*NPOIN2)))
!     &      +0.22D0*0.125D0*EXP(0.125D0*(Z(J)-Z(I+(NPLAN-1)*NPOIN2))))
!     &     *RAY_SOL/LAMB
!!
!!           EXAMPLE OF FORMULA FOR CLEAR WATER
!!           ALL CONSTANTS MAY BE TUNED
!!           0.42D0 = 1.D0-0.58D0
!!           S0TA%ADR(IND_T)%P%R(J) =
!!           ( 0.58D0/0.35D0*EXP((Z(J)-Z(I+(NPLAN-1)*NPOIN2))/0.35D0)
!!     &      +0.42D0/23.D0 *EXP((Z(J)-Z(I+(NPLAN-1)*NPOIN2))/23.D0 ))
!!     &     *RAY_SOL/LAMB
!          ENDDO
!        ENDDO
!
      ENDIF
!***********************************************************************
!     WATER QUALITY COUPLING
!*********************************************************************** 
      IF(INCLUS(COUPLING,'WAQTEL'))THEN
!
!       ACTIVATE IMPLICIT SOURCE TERMS
!        IF(LT.EQ.1) CALL YASMI_WAQ(NTRAC,YASMI) 
!
!       MAIN ROUTINE FOR WATER QUALITY
        IF(DEBUG.GT.0) WRITE(LU,*) 'CALL OF SOURCE_WAQ'
!                     TEXP,TIMP,TN
!        CALL SOURCE_WAQ
!     &  (NPOIN3,NPOIN2,S0TA,S1TA,TA,NTRAC,WAQPROCESS,RAYEFF,IND_T,IND_S,
!     &   HN,HPROP,U,V,CF,T3_01,T3_02,T3_03,T3_04,T3_05,T3_07,T3_08,
!     &   T3_09,T3_10,T3_11,T3_12,T3_13,
!     &   T2_01,T2_02,T2_03,
!     &   PATMOS,LISTIN,GRAV,ZF,DEBUG,DUMM,DT,3,VOLU2D,NPLAN,LATIT,
!     &   LONGIT,AT,MARDAT,MARTIM,ZPROP)
!        IF(DEBUG.GT.0) WRITE(LU,*) 'BACK FROM SOURCE_WAQ'
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
!                    *****************************
                      SUBROUTINE CALCS3D_THERMICS
!                    ****************************
     & (NPOIN2,NPOIN3,IND_T,IND_S,TA,ATABOS,BTABOS,PATMOS,ATMOSEXCH,
     &  WIND,LISTIN)
!
!
!***********************************************************************
! TELEMAC2D   V7P0                                        21/09/2014
!***********************************************************************
!
!brief   COMPUTES BOUNDARY CONDITIONS FOR  WAQ THERMIC PROCESS
!         COUPLED WITH T3D
!
!history  R. ATA
!+        21/02/2016
!+        V7P2
!+       CREATION
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| AT             |-->| TIME IN SECONDS
!| DT             |-->| TIME STEP
!| DIMM           |-->| 2D OR 3D
!| IND_T          |-->| INDEX OF THE TEMPERATURE IN THE TRACER TABLE
!| LISTIN         |-->| LOGICAL FOR LISTING
!| MASSOU         |<--| MASS OF TRACER ADDED BY SOURCE TERM
!| MAXSCE         |-->| MAXIMUM NUMBER OF SOURCES
!| MAXTRA         |-->| MAXIMUM NUMBER OF TRACERS
!| NPOIN          |-->| NUMBER OF NODES IN THE MESH
!| NTRAC          |-->| NUMBER OF TRACERS
!| PATMOS         |-->| ATMOSPHERIC PRESSURE
!| TETAT          |-->| COEFFICIENT OF IMPLICITATION FOR TRACERS.
!| TEXP           |-->| EXPLICIT SOURCE TERM.
!| TIMP           |-->| IMPLICIT SOURCE TERM.
!| TN             |-->| TRACERS AT TIME N
!| TSCE           |-->| PRESCRIBED VALUES OF TRACERS AT POINT SOURCES
!| TSCEXP         |<--| EXPLICIT SOURCE TERM OF POINT SOURCES
!|                |   | IN TRACER EQUATION, EQUAL TO:
!|                |   | TSCE - ( 1 - TETAT ) TN
!| VOLU2D         |-->| BASES AREA (NON ASSEMBLED)
!| YASMI          |<--| IF YES, THERE ARE IMPLICIT SOURCE TERMS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!-----------------------------------------------------------------------
!***********************************************************************
!
      USE BIEF
      USE DECLARATIONS_SPECIAL
      USE DECLARATIONS_WAQTEL,ONLY:C_ATMOS,HREL,TAIR,NEBU,CP_EAU,RO0
      USE EXCHANGE_WITH_ATMOSPHERE
      USE INTERFACE_WAQTEL, EX_CALCS3D_THERMICS => CALCS3D_THERMICS
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN)             :: NPOIN2,NPOIN3
      INTEGER, INTENT(IN)             :: IND_T,IND_S,ATMOSEXCH
      TYPE(BIEF_OBJ), INTENT(IN)      :: TA,WIND
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: ATABOS,BTABOS
      TYPE(BIEF_OBJ), INTENT(IN)      :: PATMOS
      LOGICAL,        INTENT(IN)      :: LISTIN
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!  LOCAL VARIABLES
!
      INTEGER                     :: IPOIN2
      DOUBLE PRECISION            :: TREEL,SAL,RO,LAMB
      DOUBLE PRECISION            :: FACT,WW,WW2,A
      DOUBLE PRECISION            :: RAY_ATM,RAY_EAU,FLUX_EVAP
      DOUBLE PRECISION            :: FLUX_SENS,DEBEVAP
!
! ----------------------------------------------------------------
!
!     INITIALISATION
      CALL OS( 'X=0     ' ,X=ATABOS%ADR(IND_T)%P)
      CALL OS( 'X=0     ' ,X=BTABOS%ADR(IND_T)%P)
      IF(ATMOSEXCH.EQ.0) RETURN
      IF(ATMOSEXCH.EQ.1.OR.ATMOSEXCH.EQ.2) THEN
!
        FACT=LOG(1.D4)/LOG(5.D4)
        DO IPOIN2=1,NPOIN2
          TREEL=TA%ADR(IND_T)%P%R(NPOIN3-NPOIN2+IPOIN2)
          IF (IND_S.EQ.0) THEN
            SAL = 0.D0
          ELSE
            SAL = TA%ADR(IND_S)%P%R(NPOIN3-NPOIN2+IPOIN2)
          ENDIF
          RO = RO0*(1.D0-(7.D0*(TREEL-4.D0)**2-750.D0*SAL)*1.D-6)
          LAMB=RO*CP_EAU
!
          WW = SQRT(WIND%ADR(1)%P%R(IPOIN2)*WIND%ADR(1)%P%R(IPOIN2)
     &       + WIND%ADR(2)%P%R(IPOIN2)*WIND%ADR(2)%P%R(IPOIN2))
!         LOG LAW FOR WIND AT 2 METERS
!          WW2 = WW * LOG(2.D0/0.0002D0)/LOG(10.D0/0.0002D0)
!         WRITTEN BELOW AS:
          WW2 = WW * FACT
!         ALTERNATIVE LAW FOR WIND AT 2 METERS
!          WW2 = 0.6D0*WW
          IF(ATMOSEXCH.EQ.1) THEN
            A=(4.48D0+0.049D0*TREEL)+2021.5D0*C_ATMOS*(1.D0+WW)*
     &        (1.12D0+0.018D0*TREEL+0.00158D0*TREEL**2)/LAMB
            ATABOS%ADR(IND_T)%P%R(IPOIN2)=-A
            BTABOS%ADR(IND_T)%P%R(IPOIN2)= A*TAIR%R(IPOIN2)
          
          ELSEIF(ATMOSEXCH.EQ.2) THEN
!
!     SENSIBLE HEAT FLUXES
!
            CALL EVAPO(TREEL,TAIR%R(IPOIN2),WW2,PATMOS%R(IPOIN2),HREL,
     &                 RO,FLUX_EVAP,FLUX_SENS,DEBEVAP,C_ATMOS)
!
!     LONGWAVE HEAT FLUXES
!
            CALL SHORTRAD(TREEL,TAIR%R(IPOIN2),NEBU,HREL,RAY_ATM,
     &                    RAY_EAU)
!
!     BOUNDARY CONDITION FOR TEMPERATURE AT SURFACE
!
            ATABOS%ADR(IND_T)%P%R(IPOIN2) = 0.D0
            BTABOS%ADR(IND_T)%P%R(IPOIN2) = (RAY_ATM-RAY_EAU-FLUX_EVAP
     &                                      -FLUX_SENS)/LAMB
          ENDIF
        ENDDO
      ELSE
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) "CALCS3D_THERMICS: MODELE D'ECHANGE AVEC "
          WRITE(LU,*) "        L'ATMOSPHERE NON ENCORE PROGRAMME"
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) "CALCS3D_THERMICS: MODELE EXCHANGE WITH  "
          WRITE(LU,*) "        THE ATMOSPHERE NOT IMPLEMENTED YET"
        ENDIF
        CALL PLANTE(1)
        STOP
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN

      END SUBROUTINE
      
