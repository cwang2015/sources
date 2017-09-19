    SUBROUTINE KISUPERCAM(NTOTAL,IMAT,PARAM,S,DEDT,DOMGDT,DSDT,   &
    ZR,ZRS,ZZETA,ZHARBETA,TNONE,VP0,EVI)

    IMPLICIT REAL*8 (A-H,O-Z)

    COMMON /SUPER1/ZR(IMZ),ZRS(IMZ),ZZETA(IMZ),ZHARBETA(6,IMZ)
    COMMON /CRITICAL/S1S3(IMZ)
    COMMON /TNSAND/DEPSON(6),E(3),ZKT0
    COMMON /TN1/TNONE(IMZ),VP(IMZ),VP0(IMZ)
    DIMENSION DEP(6,6),SIG(6),DE(6)
    DIMENSION HARBETA(6),ETA(6),ETASTAR(6),DHARBETA(6),ETABB(6)
    REAL MM

    ZRAMDA=param(1,IMAT)
    ZKAPA=param(2,IMAT)
    EPS0=param(3,IMAT)
    POI=param(4,IMAT)
    ZM=param(5,IMAT)
    ZMS=param(6,IMAT)
    ZBR=param(7,IMAT)
    ZMB=param(8,IMAT)
    ZMF=param(9,IMAT)
    OCR=param(10,IMAT)
    ZETA0=param(11,IMAT)
    RS0=param(12,IMAT)
    SIG00=param(13,IMAT)

    DO ID=1,NTOTAL
        SIGMM=-(S(1,1,ID)+S(2,2,ID)+S(3,3,ID))/3.0D0
        MM=3.0D0*(ZMF-1.0D0)/(ZMF+2.0D0)
        RT3=1.7320508100D0
        RT2=1.4142135624D0
        ZMB=0.95
        HEE=(ZRAMDA-ZKAPA)/(1.0D0+EPS0)
        E0=3.0D0*(1.0D0-2.0D0*POI)*(1.0D0+EPS0)*SIGMM/ZKAPA
        DLAMDA=E0*POI/(1.0D0+POI)/(1.0D0-2.0D0*POI)
        DMU=E0/2.0D0/(1.0D0+POI)
        DMU2=DMU*2.0D0

        DO I=1,6
            DO J=1,6
                DEP(I,J)=0.00
            ENDDO
        ENDDO
        DEP(1,1)=DLAMDA+DMU2
        DEP(2,2)=DLAMDA+DMU2
        DEP(3,3)=DLAMDA+DMU2
        DEP(1,2)=DLAMDA
        DEP(1,3)=DLAMDA
        DEP(2,1)=DLAMDA
        DEP(2,3)=DLAMDA
        DEP(3,1)=DLAMDA
        DEP(3,2)=DLAMDA
        DEP(4,4)=DMU
        DEP(5,5)=DMU
        DEP(6,6)=DMU

        SX=S(1,1,ID)
        SY=S(2,2,ID)
        SZ=S(3,3,ID)
        DSX=SX-SIGMM
        DSY=SY-SIGMM
        DSZ=SZ-SIGMM
        DSXY=SIG(4)
        DSYZ=SIG(5)
        DSZX=SIG(6)

        SJJ =0.500*(DSX**2+DSY**2+DSZ**2)+DSXY**2
        &    +DSYZ**2+ DSZX**2
        SJ = DSQRT(SJJ)

        ETA(1)=DSX/SIGMM
        ETA(2)=DSY/SIGMM
        ETA(3)=DSZ/SIGMM
        ETA(4)=DSXY/SIGMM
        ETA(5)=DSYZ/SIGMM
        ETA(6)=DSZX/SIGMM

        ETAVAL=0.D0
        DO I=1,3
            ETAVAL=ETAVAL+ETA(I)*ETA(I)+2.D0*ETA(I+3)*ETA(I+3)
        END DO
        ETAVAL=DSQRT(3.0/2.0*ETAVAL)

        DO I=1,6
            HARBETA(I)=ZHARBETA(I,ID)  !ƒÀij
            ETASTAR(I)=ETA(I)-HARBETA(I)  !ƒÅ*ij=ƒÅij-ƒÀij
        END DO

        ETASTARVAL=0.D0
        DO I=1,3
            ETASTARVAL=ETASTARVAL+ETASTAR(I)**2+2.D0*ETASTAR(I+3)**2
        END DO
        ETASTARVAL=DSQRT(3.D0/2.D0*ETASTARVAL)

        EBS=0.D0
        DO I=1,3
            EBS=EBS+ETASTAR(I)*HARBETA(I)+2.D0*ETASTAR(I+3)*HARBETA(I+3)
        END DO    ! EBS=ƒÅ*ijƒÀij
        ZBO=ZMB*RT2/RT3*MM

        DO I=1,6
            IF(ETASTARVAL.LE.1.0E-6) THEN
                ETABB(I)=0.D0
            ELSE
                ETABB(I)=ZBO*RT3/RT2*ETASTAR(I)/ETASTARVAL-HARBETA(I)
            END IF
        END DO
        ETABBVALUE=0.D0

        DO I=1,3
            ETABBVALUE=ETABBVALUE+ETABB(I)*ETABB(I)
            &    +2.D0*ETABB(I+3)*ETABB(I+3)
        END DO
        ETABBVALUE=DSQRT(ETABBVALUE)

        IF(ETABBVALUE.NE.0.0D0)THEN
            DO I=1,6
                ETABB(I)=ETABB(I)/ETABBVALUE
            END DO
        END IF

        SETAB=0.D0
        DO I=1,3
            SETAB=SETAB+ETASTAR(I)*ETABB(I)+2.D0*ETASTAR(I+3)*ETABB(I+3)
        END DO

        SBETAB=0.D0
        DO I=1,3
            SBETAB=SBETAB+HARBETA(I)*ETABB(I)+2.D0*HARBETA(I+3)*ETABB(I+3)
        END DO

        PM=SIGMM*(MM**2-ZZETA(ID)**2+ETASTARVAL**2)
        &/(MM**2-ZZETA(ID)**2)
        VPT=VP(ID)+VP0(ID)
        PMSTAR=TNONE(ID)*EXP(VPT/HEE)
        PMBAR=PMSTAR/ZRS(ID)
        ZR(ID)=PM/PMBAR

        TT1=1.0D0/(MM*MM-ZZETA(ID)*ZZETA(ID)
        &+ETASTARVAL*ETASTARVAL)/SIGMM
        TT2=MM*MM-ETAVAL*ETAVAL
        F11=TT1*(3.D0*ETASTAR(1)+TT2/3.D0)
        F22=TT1*(3.D0*ETASTAR(2)+TT2/3.D0)
        F33=TT1*(3.D0*ETASTAR(3)+TT2/3.D0)
        F12=6.D0*TT1*ETASTAR(4)
        F23=6.D0*TT1*ETASTAR(5)
        F31=6.D0*TT1*ETASTAR(6)

        HDD=TT1*TT2    !Ýf/ÝƒÐm

        SS1=MM*MM*SETAB-ETASTARVAL**2*SBETAB-ZZETA(ID)**2*SETAB
        STRSTD=(SIGMM/98.D0)**2+0.01

        TMS=TT2+6.D0*MM*ZBR*(1.D0-(ETAVAL/MM))*(ZMB*MM-ZZETA(ID))
        &*SS1*ETASTARVAL/(MM*MM-ZZETA(ID)**2)
        &/(MM*MM-ZZETA(ID)**2+ETASTARVAL**2)
        &-2.D0*ZMS*MM*(1.D0-ZRS(ID))*ETASTARVAL
        &-ZM*MM*DLOG(ZR(ID))/ZR(ID)
        &*DSQRT(6.D0*ETASTARVAL**2+TT2*TT2/3.D0)*(STRSTD/(STRSTD+1.D0))

        HHH=TT1*TMS/HEE
        FIJ=2.0D0*(F11**2+F22**2+F33**2)+F12**2+F23**2+F31**2
        HHH=DLAMDA*(HDD**2)+DMU*FIJ+HHH
        D=1.0D0/HHH
        DF=DLAMDA*HDD
        DG=DLAMDA*HDD
        
        RAMDA=(DMU2*F11+DLAMDA*HDD)*DEDT(1,1,ID)
        &+(DMU2*F22+DLAMDA*HDD)*DEDT(2,2,ID)
        &+(DMU2*F33+DLAMDA*HDD)*DEDT(3,3,ID)
        &+ DMU*(F12*DEDT(1,2,ID)+F23*DEDT(2,3,ID)+F31*DEDT(3,1,ID))
        GAMA=RAMDA*D

        IF(GAMA >= 0.D0) THEN
            
            DRS=2.0D0*MM*ZMS*ZRS(ID)*(1.D0-ZRS(ID))*ETASTARVAL
            &*GAMA*TT1/HEE
            ZRS(ID)=ZRS(ID)+DRS
            
            TT3=2.D0*MM*ZBR*TT1*GAMA/HEE*ETASTARVAL*(ZMB*MM-ZZETA(ID))
            DO I=1,6
                DHARBETA(I)=TT3*ETABB(I)
                ZHARBETA(I,ID)=ZHARBETA(I,ID)+DHARBETA(I)
            END DO
            
            ZZETA(ID)=0.D0
            DO I=1,3
                ZZETA(ID)=ZZETA(ID)+ZHARBETA(I,ID)**2
                &    +2.D0*ZHARBETA(I+3,ID)**2
            END DO
            ZZETA(ID)=DSQRT(3.D0/2.D0*ZZETA(ID))
            
            DEPSON(1)=GAMA*F11
            DEPSON(2)=GAMA*F22
            DEPSON(3)=GAMA*F33
            DEPSON(4)=GAMA*F12
            DEPSON(5)=GAMA*F23
            DEPSON(6)=GAMA*F31
            DVP=DEPSON(1)+DEPSON(2)+DEPSON(3)
            VP(ID)=VP(ID)+DVP
            
            QQ1=DG+DMU2*F11
            QQ2=DG+DMU2*F22
            QQ3=DG+DMU2*F33
            QQ4=DMU*F12
            QQ5=DMU*F23
            QQ6=DMU*F31
            RR1=DF+DMU2*F11
            RR2=DF+DMU2*F22
            RR3=DF+DMU2*F33
            RR4=DMU*F12
            RR5=DMU*F23
            RR6=DMU*F31
            
            DEP(1,1)=DEP(1,1)-D*QQ1*RR1
            DEP(1,2)=DEP(1,2)-D*QQ1*RR2
            DEP(1,3)=DEP(1,3)-D*QQ1*RR3
            DEP(1,4)=DEP(1,4)-D*QQ1*RR4
            DEP(1,5)=DEP(1,5)-D*QQ1*RR5
            DEP(1,6)=DEP(1,6)-D*QQ1*RR6
            DEP(2,1)=DEP(2,1)-D*QQ2*RR1
            DEP(2,2)=DEP(2,2)-D*QQ2*RR2
            DEP(2,3)=DEP(2,3)-D*QQ2*RR3
            DEP(2,4)=DEP(2,4)-D*QQ2*RR4
            DEP(2,5)=DEP(2,5)-D*QQ2*RR5
            DEP(2,6)=DEP(2,6)-D*QQ2*RR6
            DEP(3,1)=DEP(3,1)-D*QQ3*RR1
            DEP(3,2)=DEP(3,2)-D*QQ3*RR2
            DEP(3,3)=DEP(3,3)-D*QQ3*RR3
            DEP(3,4)=DEP(3,4)-D*QQ3*RR4
            DEP(3,5)=DEP(3,5)-D*QQ3*RR5
            DEP(3,6)=DEP(3,6)-D*QQ3*RR6
            DEP(4,1)=DEP(4,1)-D*QQ4*RR1
            DEP(4,2)=DEP(4,2)-D*QQ4*RR2
            DEP(4,3)=DEP(4,3)-D*QQ4*RR3
            DEP(4,4)=DEP(4,4)-D*QQ4*RR4
            DEP(4,5)=DEP(4,5)-D*QQ4*RR5
            DEP(4,6)=DEP(4,6)-D*QQ4*RR6
            DEP(5,1)=DEP(5,1)-D*QQ5*RR1
            DEP(5,2)=DEP(5,2)-D*QQ5*RR2
            DEP(5,3)=DEP(5,3)-D*QQ5*RR3
            DEP(5,4)=DEP(5,4)-D*QQ5*RR4
            DEP(5,5)=DEP(5,5)-D*QQ5*RR5
            DEP(5,6)=DEP(5,6)-D*QQ5*RR6
            DEP(6,1)=DEP(6,1)-D*QQ6*RR1
            DEP(6,2)=DEP(6,2)-D*QQ6*RR2
            DEP(6,3)=DEP(6,3)-D*QQ6*RR3
            DEP(6,4)=DEP(6,4)-D*QQ6*RR4
            DEP(6,5)=DEP(6,5)-D*QQ6*RR5
            DEP(6,6)=DEP(6,6)-D*QQ6*RR6
            
            DSDT(1,1,ID)=DEP(1,1)*(-DEDT(1,1,ID))+DEP(1,2)*(-DEDT(2,2,ID))+DEP(1,3)*(-DEDT(3,3,ID))+DEP(1,4)*2.D0*DEDT(1,2,ID)+DEP(1,5)*2.D0*DEDT(2,3,ID)+DEP(1,6)*2.D0*DEDT(3,1,ID)
            DSDT(1,2,ID)=DEP(4,1)*(-DEDT(1,1,ID))+DEP(4,2)*(-DEDT(4,2,ID))+DEP(4,3)*(-DEDT(3,3,ID))+DEP(4,4)*2.D0*DEDT(1,2,ID)+DEP(4,5)*2.D0*DEDT(2,3,ID)+DEP(4,6)*2.D0*DEDT(3,1,ID)
            DSDT(2,1,ID)=DSDT(1,2,ID)
            DSDT(2,2,ID)=DEP(2,1)*(-DEDT(1,1,ID))+DEP(2,2)*(-DEDT(2,2,ID))+DEP(2,3)*(-DEDT(3,3,ID))+DEP(2,4)*2.D0*DEDT(1,2,ID)+DEP(2,5)*2.D0*DEDT(2,3,ID)+DEP(2,6)*2.D0*DEDT(3,1,ID)
            DSDT(3,3,ID)=DEP(3,1)*(-DEDT(1,1,ID))+DEP(3,2)*(-DEDT(2,2,ID))+DEP(3,3)*(-DEDT(3,3,ID))+DEP(3,4)*2.D0*DEDT(1,2,ID)+DEP(3,5)*2.D0*DEDT(2,3,ID)+DEP(3,6)*2.D0*DEDT(3,1,ID)
        ENDIF
    ENDDO
    RETURN
    END
