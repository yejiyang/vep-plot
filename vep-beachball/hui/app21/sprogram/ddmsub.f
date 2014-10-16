C------ TO GET THE ELEMENT SIZE
        SUBROUTINE SIZEELEM(NNE,NCOOR,KNODE,NODI,COOR,ELEMSIZE)
        IMPLICIT real*8 (A-H,O-Z)
        DIMENSION NODI(4),COOR(NCOOR,KNODE)
6       FORMAT(2X,10I6)
7       FORMAT(2X,6E12.5)
 
        ELEMSIZE=1.E+20
        DO INOD=1,NNE
          DO JNOD=INOD+1,NNE
            XSKP=0.0
            DO NR=1,NCOOR
              XSKP=XSKP+(COOR(NR,NODI(INOD))-COOR(NR,NODI(JNOD)))**2.
            ENDDO
            XSKP=SQRT(XSKP)
            ELEMSIZE=MIN(ELEMSIZE,XSKP)
          ENDDO
        ENDDO
CCC        WRITE(*,*) 'NE,ELEMSIZE ===',NE,ELEMSIZE
        RETURN
        END
 
C------ TO GET THE MINIUM AND MAXMIUM OF THE ELEMENT.
        SUBROUTINE MINMAXELEM(MULNO,KELEM,NNODE,NNE,NE,NCOOR,
     *  KNODE,NODI,NODE,COOR,XE,XMIN,XMAX)
        IMPLICIT real*8 (A-H,O-Z)
        DIMENSION NODI(4),NODE(KELEM),COOR(NCOOR,KNODE),XE(3,4),
     *  XMIN(3),XMAX(3)
6       FORMAT(2X,8I8)
7       FORMAT(2X,6E12.5)
 
C------ MULNO STANDS FOR THE NUMBER OF MULTIPLIERS.
        MULNO=NODE(NE*NNODE)
        DO INOD=1,NNE
          NODI(INOD)=NODE((NE-1)*NNODE+INOD)
          DO NR=1,NCOOR
          XE(NR,INOD)=COOR(NR,NODI(INOD))
          ENDDO
        ENDDO
C        WRITE(*,*) 'NE,NODI =======',NE
C        WRITE(*,6) (NODI(INOD),INOD=1,NNE)
C-------- XMIN(NR) STAND FOR THE MINIMIUM VALUE OF THE ELEMENT.
C-------  XMAX(NR) STAND FOR THE MAXMIUM  VALUE OF THE ELEMENT.
C--------------------------------------------------------------
C------- TO GET XMIN(NR) AND XMAX(NR)
        DO NR=1,NCOOR
          XMIN(NR)=XE(NR,1)
          XMAX(NR)=XE(NR,1)
        ENDDO
        DO INOD=2,NNE
          DO NR=1,NCOOR
          XMIN(NR)=MIN(XMIN(NR),XE(NR,INOD))
          XMAX(NR)=MAX(XMAX(NR),XE(NR,INOD))
          ENDDO
        ENDDO
        DO NR=1,NCOOR
          ABSMIN=XMIN(NR)/10.
          IF (ABSMIN.LT.0) ABSMIN=-ABSMIN
          ABSMAX=XMAX(NR)/10.
          IF (ABSMAX.LT.0) ABSMAX=-ABSMAX
          XSKIP=MAX(ABSMIN,dble(0.1))
          XSKIP=MAX(ABSMAX,XSKIP)
          XMIN(NR)=XMIN(NR)-XSKIP
          XMAX(NR)=XMAX(NR)+XSKIP
        ENDDO
        RETURN
        END
 
 
C...... TO GET THE DIRECTION OF THE FORCES OR THE LOCAL COORDINATE
        SUBROUTINE GETD(NNE,NCOOR,XE,DE)
        IMPLICIT real*8 (A-H,O-Z)
        DIMENSION XE(3,4),DE(3,3),XD(3),YD(3),ZD(3)
6       FORMAT(1X,8I6)
7       FORMAT(1X,6E12.5)
        KGOTO=NNE-1
        GOTO (1,2,3) KGOTO
1       CONTINUE
        DO NR=1,NCOOR
          XD(NR)=XE(NR,2)-XE(NR,1)
        ENDDO
        YD(1)=-XD(2)
        YD(2)= XD(1)
        XSUM=0.0
        DO NR=1,NCOOR
          XSUM=XSUM+YD(NR)*YD(NR)
        ENDDO
        XSUM=SQRT(XSUM)
        IF (XSUM.LE.1.E-9) THEN
          WRITE(*,*) 'THE LENGTH OF LAGRANGE ELEM IS ZERO =='
          WRITE(*,*) 'XSUM ===',XSUM
          STOP 2222
        ENDIF
        DO NR=1,NCOOR
          XD(NR)=-XD(NR)/XSUM
          YD(NR)= YD(NR)/XSUM
        ENDDO
        DO NR=1,NCOOR
          DE(NR,1)=YD(NR)
          DE(NR,2)=XD(NR)
        ENDDO
        RETURN
2       CONTINUE
CC        WRITE(*,*) 'NNE EQ 3 ======',NNE
        DO NR=1,NCOOR
          XD(NR)=XE(NR,2)-XE(NR,1)
          YD(NR)=XE(NR,3)-XE(NR,1)
        ENDDO
        CALL GETE3(XD,YD,ZD)
        CALL GETE3(ZD,XD,YD)
        CALL GETEC(XD,YD,ZD)
        DO NR=1,NCOOR
          DE(NR,1)=ZD(NR)
          DE(NR,2)=XD(NR)
          DE(NR,3)=YD(NR)
        ENDDO
        RETURN
3       CONTINUE
C        WRITE(*,*) 'NNE EQ 4 =========',NNE
        DO NR=1,NCOOR
          XD(NR)=XE(NR,2)-XE(NR,1)
          YD(NR)=XE(NR,3)-XE(NR,1)
        ENDDO
C        WRITE(*,*) 'XD,YD ======='
C        WRITE(*,7) (XD(II),II=1,NCOOR)
C        WRITE(*,7) (YD(II),II=1,NCOOR)
        CALL GETE3(XD,YD,ZD)
        CALL GETE3(ZD,XD,YD)
        CALL GETEC(XD,YD,ZD)
        DO NR=1,NCOOR
          DE(NR,1)=ZD(NR)
          DE(NR,2)=XD(NR)
          DE(NR,3)=YD(NR)
        ENDDO
        RETURN
        END
 
C------ TO GET THE PROJECT POINT OF XMAST TO THE ELEM ------
C------ XMAST(NR)=XPROJ(NR)+DE(NR,1)*DIST
C------ DIST=(XMAST(NR)-XE(NR,1))*DE(NR,1)
        SUBROUTINE PROJECT(NNE,NCOOR,XE,DE,XMAST,XPROJ,DIST)
        IMPLICIT real*8(A-H,O-Z)
        DIMENSION XE(3,4),DE(3,3),XMAST(3),XPROJ(3)
6       FORMAT(1X,10I5)
7       FORMAT(1X,6E12.5)
        DIST=0.0
        DIST4=0.0
        DO NR=1,NCOOR
          DIST =DIST +(XMAST(NR)-XE(NR,1))*DE(NR,1)
          IF (NNE.EQ.4) DIST4=DIST4+(XE(NR,4)-XE(NR,1))*DE(NR,1)
        ENDDO
        DO NR=1,NCOOR
          XPROJ(NR)=XMAST(NR)-DE(NR,1)*DIST
          IF (NNE.EQ.4) XE(NR,4)=XE(NR,4)-DE(NR,1)*DIST4
        ENDDO
        RETURN
        END
 
        SUBROUTINE POINTINLINE(KEY1,NCOOR,XELINE,XPROJ,
     *  KID,DISTT,KSHAP)
        IMPLICIT real*8 (A-H,O-Z)
        DIMENSION XELINE(3,2),XPROJ(3),XPROJ2(3),TDIRECT(3)
C------ TDIRECT=XELINE(2)-XELINE(1)
C------ TDIRECT=TDIRECT/|TDIRECT|
C------ XPROJ2=XLINE(1)+DISTT*TDIRECT
C------ (XPROJ-XPROJ2)*TDIRECT=0
C------ DISTT=(XPROJ-XELINE(1))*TDIRECT
C------ DISTN=SQRT(XPROJ-XPROJ2)
 
6       FORMAT(2X,8I8)
7       FORMAT(2X,6E12.5)
 
        DO NR=1,NCOOR
          TDIRECT(NR)=XELINE(NR,2)-XELINE(NR,1)
        ENDDO
        TDIRECTLEN=0.0
        DO NR=1,NCOOR
          TDIRECTLEN=TDIRECTLEN+TDIRECT(NR)**2.
        ENDDO
        TDIRECTLEN=SQRT(TDIRECTLEN)
        DO NR=1,NCOOR
          TDIRECT(NR)=TDIRECT(NR)/TDIRECTLEN
        ENDDO
        DISTT=0.0
        DO NR=1,NCOOR
          DISTT=DISTT+(XPROJ(NR)-XELINE(NR,1))*TDIRECT(NR)
        ENDDO
        DO NR=1,NCOOR
          XPROJ2(NR)=XELINE(NR,1)+DISTT*TDIRECT(NR)
        ENDDO
        DISTN=0.0
        DO NR=1,NCOOR
          DISTN=DISTN+(XPROJ(NR)-XPROJ2(NR))**2.
        ENDDO
        DISTN=SQRT(DISTN)
        DISTT=DISTT/TDIRECTLEN
        DISTN=DISTN/TDIRECTLEN
CC        WRITE(*,*) 'DISTN,DISTT ==',DISTN,DISTT
CC        WRITE(*,*) 'XPROJ,XPROJ2 ==='
CC        WRITE(*,7) (XPROJ(NR),NR=1,NCOOR),(XPROJ2(NR),NR=1,NCOOR)
        IF (DISTN.LE.2.E-5.AND.DISTT.GE.-2.E-5.AND.
     *  DISTT.LE.(1.0+2.E-5)) THEN
          KEY1=0
          IF (ABS(DISTT).LE.2.E-5) THEN
            DO NR=1,NCOOR
              XPROJ(NR)=XELINE(NR,1)
            ENDDO
            DISTT=0.0
          ENDIF
          IF (ABS(DISTT-1.0).LE.2.E-5) THEN
            DO NR=1,NCOOR
              XPROJ(NR)=XELINE(NR,2)
            ENDDO
            DISTT=1.
          ENDIF
          IF (DISTT.GT.2.E-5.AND.DISTT.LE.(1.0-2.E-5)) THEN
            DO NR=1,NCOOR
              XPROJ(NR)=XPROJ2(NR)
            ENDDO
          ENDIF
        ELSEIF (KID.EQ.-1.AND.DISTN.LE.0.01.AND.
     *  DISTT.GE.-0.01.AND.DISTT.LE.(1.0+0.01)) THEN
          KEY1=0
          IF (KSHAP.EQ.2) THEN
            IF (ABS(DISTT).LE.0.01) THEN
              DO NR=1,NCOOR
                XPROJ(NR)=XELINE(NR,1)
              ENDDO
              DISTT=0.0
            ENDIF
            IF (ABS(DISTT-1.).LE.0.01) THEN
              DO NR=1,NCOOR
                XPROJ(NR)=XELINE(NR,2)
              ENDDO
              DISTT=1.
            ENDIF
            IF (DISTT.GT.0.01.AND.DISTT.LE.(1.0-0.01)) THEN
              DO NR=1,NCOOR
                XPROJ(NR)=XPROJ2(NR)
              ENDDO
            ENDIF
          ENDIF
        ENDIF
        RETURN
        END
 
        SUBROUTINE POINTINTRI(KEY1,NNE,NCOOR,XE,XPROJ,UE)
        IMPLICIT real*8 (A-H,O-Z)
        DIMENSION XE(3,4),XPROJ(3),XP(3,3),UE(4)
6       FORMAT(2X,10I5)
7       FORMAT(2X,6E12.5)
 
        DO INOD=1,NNE
          INOD1=INOD+1
          IF (INOD1.GT.NNE) INOD1=INOD1-NNE
          DO NR=1,NCOOR
          XP(NR,1)=XPROJ(NR)
          XP(NR,2)=XE(NR,INOD)
          XP(NR,3)=XE(NR,INOD1)
          ENDDO
          CALL AREATRI(NCOOR,XP,SS)
          INOD2=INOD1+1
          IF (INOD2.GT.NNE) INOD2=INOD2-NNE
          UE(INOD2)=SS
        ENDDO
        DO NR=1,NCOOR
          XP(NR,1)=XE(NR,1)
          XP(NR,2)=XE(NR,2)
          XP(NR,3)=XE(NR,3)
        ENDDO
        CALL AREATRI(NCOOR,XP,SS)
        SSUM=0.0
        DO INOD=1,NNE
          SSUM=SSUM+UE(INOD)
        ENDDO
        XKEY=ABS(SSUM-SS)
        IF (SS.GT.1) XKEY=XKEY/SS
CC        WRITE(*,*) 'UE,SSUM,SS,XKEY =='
CC        WRITE(*,7) (UE(INOD),INOD=1,NNE),SSUM,SS,XKEY
        IF (XKEY.LE.1.E-6) THEN
          KEY1=0
          DO INOD=1,NNE
            UE(INOD)=UE(INOD)/SS
          ENDDO
        ENDIF
        RETURN
        END
 
        SUBROUTINE POINTIN(KEY1,NNE,NCOOR,KID,XE,XPROJ,UE,KSHAP)
        IMPLICIT real*8 (A-H,O-Z)
        DIMENSION XE(3,4),XPROJ(3),XE3(3,4),UE(4),XO(3),
     *  XELINE(3,2),UE3(4)
6       FORMAT(1X,8I6)
7       FORMAT(1X,6E12.5)
        KGOTO=NNE-1
 
        GOTO(1,2,3) KGOTO
 
1       CONTINUE
        DO INOD=1,NNE
          DO NR=1,NCOOR
            XELINE(NR,INOD)=XE(NR,INOD)
          ENDDO
        ENDDO
        CALL POINTINLINE(KEY1,NCOOR,XELINE,XPROJ,KID,DISTT,KSHAP)
        IF (KEY1.EQ.0.AND.KSHAP.EQ.2) THEN
          UE(1)=1.-DISTT
          UE(2)=DISTT
CC          WRITE(*,*) 'UE ======'
CC          WRITE(*,7) (UE(I),I=1,NNE)
        ENDIF
        RETURN
 
2       CONTINUE
C        WRITE(*,*) 'NNE EQ 3 ===',NNE
        DO INOD=1,NNE
          INOD1=INOD+1
          IF (INOD1.GT.NNE) INOD1=INOD1-NNE
          DO NR=1,NCOOR
            XELINE(NR,1)=XE(NR,INOD)
            XELINE(NR,2)=XE(NR,INOD1)
          ENDDO
          CALL POINTINLINE(KEY1,NCOOR,XELINE,XPROJ,KID,DISTT,KSHAP)
          IF (KEY1.EQ.0.) THEN
            IF (KSHAP.EQ.2) THEN
              DO JNOD=1,NNE
                UE(JNOD)=0.0
              ENDDO
              UE(INOD)=1.-DISTT
              UE(INOD1)=DISTT
            ENDIF
            GOTO 15
          ENDIF
        ENDDO
15      CONTINUE
        IF (KEY1.EQ.0) GOTO 20
        CALL POINTINTRI(KEY1,NNE,NCOOR,XE,XPROJ,UE)
20      CONTINUE
        RETURN
 
3       CONTINUE
CC        WRITE(*,*) 'NNE EQ 4 ===',NNE
        DO INOD=1,NNE
          INOD1=INOD+1
          IF (INOD1.GT.NNE) INOD1=INOD1-NNE
          DO NR=1,NCOOR
            XELINE(NR,1)=XE(NR,INOD)
            XELINE(NR,2)=XE(NR,INOD1)
          ENDDO
          CALL POINTINLINE(KEY1,NCOOR,XELINE,XPROJ,KID,DISTT,KSHAP)
          IF (KEY1.EQ.0) THEN
            IF (KSHAP.EQ.2) THEN
              DO JNOD=1,NNE
                UE(JNOD)=0.0
              ENDDO
CCC              WRITE(*,*) 'INOD,DISTT ===',INOD,DISTT
              UE(INOD)=1.-DISTT
              UE(INOD1)=DISTT
            ENDIF
            GOTO 25
          ENDIF
        ENDDO
25      CONTINUE
        IF (KEY1.EQ.0) GOTO 30
        DO NR=1,NCOOR
          XO(NR)=0.0
        ENDDO
        DO INOD=1,NNE
          DO NR=1,NCOOR
          XO(NR)=XO(NR)+XE(NR,INOD)
          ENDDO
        ENDDO
        DO NR=1,NCOOR
          XO(NR)=XO(NR)/NNE
        ENDDO
        DO 100 INOD=1,NNE
          INOD1=INOD+1
          IF (INOD1.GT.NNE) INOD1=INOD1-NNE
          DO NR=1,NCOOR
            XE3(NR,1)=XE(NR,INOD)
            XE3(NR,2)=XE(NR,INOD1)
            XE3(NR,3)=XO(NR)
          ENDDO
          NNE1=NNE-1
CC          WRITE(*,*) 'CALL POINTINTRI,INOD =',INOD
          CALL POINTINTRI(KEY1,NNE1,NCOOR,XE3,XPROJ,UE3)
CC          WRITE(*,*) 'INOD,KEY1 ===',INOD,KEY1
          IF (KEY1.EQ.0) THEN
            IF (KSHAP.EQ.2) THEN
              DO JNOD=1,NNE
                UE(JNOD)=UE3(3)/NNE
              ENDDO
              UE(INOD)=UE(INOD)+UE3(1)
              UE(INOD1)=UE(INOD1)+UE3(2)
            ENDIF
            GOTO 99
          ENDIF
100     CONTINUE
99      CONTINUE
30      CONTINUE
CC        WRITE(*,*) 'UE =========='
CC        WRITE(*,7) (UE(INOD),INOD=1,NNE)
        RETURN
        END
 
C...... TO CALCULATE THE AREA OF THE TRIANGLE
C...... |a*b|=|a|.|b|.sin(a^b)
C...... S=|a*b|/2
        SUBROUTINE AREATRI(NCOOR,XP,SS)
        IMPLICIT real*8 (A-H,O-Z)
        DIMENSION XP(3,3),XA(3),XB(3),XC(3)
6       FORMAT(1X,10I5)
7       FORMAT(1X,6E12.5)
        DO NR=1,NCOOR
        XA(NR)=XP(NR,2)-XP(NR,1)
        XB(NR)=XP(NR,3)-XP(NR,1)
        ENDDO
        CALL GETE3(XA,XB,XC)
        SS=0.0
        DO NR=1,NCOOR
          SS=SS+XC(NR)*XC(NR)
        ENDDO
        SS=SQRT(SS)/2.
        RETURN
        END
 
        SUBROUTINE SHORTDIST(KEY3,NCOOR,LMSUMC,NMDOF,DISTNEW,
     *  IC1,LGNE,GAP,DIRECT,DE)
        IMPLICIT real*8 (A-H,O-Z)
        DIMENSION GAP(NMDOF,LMSUMC),DIRECT(9,LMSUMC),DE(3,3)
6       FORMAT(2X,10I5)
7       FORMAT(2X,6E12.5)
 
        IF (LGNE.eq.1) THEN
          KEY3=2
        ELSE
C------   TO CALCULATE THE ANGLE BETWEEN TWO NORMAL DIRECTION.
          PROD=0.0
          DO NR=1,NCOOR
            PROD=PROD+DE(NR,1)*DIRECT(NR,IC1)
          ENDDO
          IF (PROD.GT.1.0) PROD= 1.
          ANGLESUB=ACOS(PROD)
          ANGLEKEY=5.*3.1415926/180.
CC          WRITE(*,*) 'ANGLESUB,ANGLEKEY ==',ANGLESUB,ANGLEKEY
C------   THE STANDAND OF DETERMINION IS 5 DEGREE.
          IF (ANGLESUB.LE.ANGLEKEY) THEN
            KEY3=2
            IF (DISTNEW.LT.GAP(1,IC1)) KEY3=1
          ENDIF
        ENDIF
CC        WRITE(*,*) '*****  KEYDEAL, KEY3 *****',KEYDEAL,KEY3
        RETURN
        END
 
        SUBROUTINE GETIDL(LMSUMC,KDGOF,KNODE,NNE,NCOOR,
     *  NMDOF,IC,KDGOFL,KBSYS0,
     *  NODVAR,LSYS0,IDLM,IDLS,NODI,DE,UE)
        IMPLICIT real*8 (A-H,O-Z)
        DIMENSION NODVAR(KDGOF,KNODE),LSYS0(KBSYS0,LMSUMC),IDLM(3),
     *  IDLS(3),NODI(4),IDSC(3),UE(4),DE(3,3)
C------ TO DEAL WITH THE MULTI. WHICH ARE RELATED TO
C------ CONSTRAINED POINTS.FIRSTLY, WE DEAL WITH THE
C------ MASTER POINT, THEN WE DEAL WITH THE SLAVE POINTS.
6       FORMAT(2X,8I6)
7       FORMAT(2X,6E12.5)
        DO NDF=1,NMDOF
          IDLM(NDF)=1
          IDLS(NDF)=1
        ENDDO
        IF (KDGOFL.EQ.1) THEN
          IF (LSYS0(4,IC).LT.0) THEN
            DO NDF=1,NMDOF
              IDLM(NDF)=-1
            ENDDO
          ENDIF
          SHAPVSUM=0.0
          DO I=1,NNE
            IF (NODVAR(1,NODI(I)).LE.0) SHAPVSUM=SHAPVSUM+UE(I)
          ENDDO
          RNNE=1.-1./(NNE*1.)
CC          WRITE(*,*) 'UE,SHAPVSUM,RNNE =='
CC          WRITE(*,7) (UE(INOD),INOD=1,NNE),SHAPVSUM,RNNE
          IF (SHAPVSUM.GE.RNNE) THEN
            DO NDF=1,NMDOF
              IDLS(NDF)=-1
            ENDDO
          ENDIF
        ENDIF
        IF (KDGOFL.GT.1) THEN
C--------------  THE MASTER POINT -----------------------
        KEY4=0
        DO NR=1,NCOOR
          IF (LSYS0(3+NR,IC).LE.0) KEY4=KEY4+1
        ENDDO
CC        WRITE(*,*) 'KEY4 =======',KEY4
        IF (KEY4.EQ.NCOOR) THEN
           DO NDF=1,NMDOF
             IDLM(NDF)=-1
           ENDDO
        ENDIF
        IF (KEY4.GT.1.AND.KEY4.LT.NCOOR) THEN
C-------  TO DEAL WITH ONLY SINGLE DIRECTION
          DO NR=1,NCOOR
            IF (LSYS0(3+NR,IC).LE.0) THEN
              DO NDF=1,NMDOF
                XSKP=ABS(ABS(DE(NR,NDF))-1)
                IF (XSKP.LE.0.01) IDLM(NDF)=-1
              ENDDO
            ELSE
              JNOD=NR
            ENDIF
          ENDDO
C-------  TO DEAL WITH THE TWO DIRECTIONS
          DO NDF=1,NMDOF
            XSKP=ABS(DE(JNOD,NDF))
            IF (XSKP.LE.0.01) IDLM(NDF)=-1
          ENDDO
        ENDIF
        IF (KEY4.EQ.1) THEN
          DO NR=1,NCOOR
            IF (LSYS0(3+NR,IC).LE.0) THEN
              DO NDF=1,NMDOF
                XSKP=ABS(ABS(DE(NR,NDF))-1)
                IF (XSKP.LE.0.01) IDLM(NDF)=-1
              ENDDO
            ENDIF
          ENDDO
        ENDIF
        IF (KEY4.EQ.0) THEN
          DO NDF=1,NMDOF
            IDLM(NDF)=1
          ENDDO
        ENDIF
C
C--------------  THE SLAVE POINTS -----------------------
C...... IDSC(1): THE NUMBER OF THE X-COOR CONSTRAINED
C...... IDSC(2): THE NUMBER OF THE Y-COOR CONSTRAINED
C...... IDSC(3): THE NUMBER OF THE Z-COOR CONSTRAINED
        DO NR=1,NCOOR
          IDSC(NR)=0.0
        ENDDO
        DO I=1,NNE
          INOD=NODI(I)
          DO NR=1,NCOOR
            IF (NODVAR(NR,INOD).LE.0) IDSC(NR)=IDSC(NR)+1
          ENDDO
        ENDDO
        IDSSUM=0
        DO NR=1,NCOOR
          IDSSUM=IDSSUM+IDSC(NR)
        ENDDO
C....... IF IDDSUM.EQ.0, THEN THERE IS NO POINT CONSTRAINED
        IF (IDSSUM.GT.0) THEN
C........ TO DEAL WITH THE SINGLE DIRECTION
          DO 100 NR=1,NCOOR
            IF (IDSC(NR).EQ.0) GOTO 100
            SHAPVSUM=0.0
            DO I=1,NNE
              IF (NODVAR(NR,NODI(I)).LE.0) SHAPVSUM=SHAPVSUM+UE(I)
            ENDDO
            RNNE=1.-1./(NNE*1.)
            IF (SHAPVSUM.GE.RNNE) THEN
              DO NDF=1,NMDOF
                XSKP=ABS(ABS(DE(NR,NDF))-1)
                IF (XSKP.LE.0.01) IDLS(NDF)=-1
              ENDDO
            ENDIF
100       CONTINUE
C........ TO DEAL WITH TWO DIRECTIONS
          DO 200 NR=1,NCOOR
            SHAPVSUM=0.0
            DO I=1,NNE
              KEY5=0
              DO JNR=1,NCOOR
                IF ((JNR.EQ.NR.AND.NODVAR(JNR,NODI(I)).LE.0)
     *              .OR.(JNR.NE.NR.AND.NODVAR(JNR,NODI(I)).GT.0)) THEN
                  KEY5=1
                ENDIF
              ENDDO
              IF (KEY5.EQ.0) SHAPVSUM=SHAPVSUM+UE(I)
            ENDDO
            RNNE=1.-1./(NNE*1.)
            IF (SHAPVSUM.GE.RNNE) THEN
              DO NDF=1,NMDOF
                XSKP=ABS(DE(NR,NDF))
                IF (XSKP.LE.0.01) IDLS(NDF)=-1
              ENDDO
            ENDIF
200       CONTINUE
C........ TO DEAL WITH THE CASE THAT ALL DIRECTIONS ARE CONSTRAINED.
          SHAPVSUM=0.0
          DO I=1,NNE
            KEY6=0
            DO JNR=1,NCOOR
              IF(NODVAR(JNR,NODI(I)).GT.0) KEY6=1
            ENDDO
            IF (KEY6.EQ.0) SHAPVSUM=SHAPVSUM+UE(I)
          ENDDO
          RNNE=1.-1./(NNE*1.)
          IF (SHAPVSUM.GE.RNNE) THEN
            DO NDF=1,NMDOF
              IDLS(NDF)=-1
            ENDDO
          ENDIF
        ENDIF
        ENDIF
C----------------------------------------------------------
          RETURN
          END
 
        SUBROUTINE VALUECONTACT(NBLK,NNE,NCOOR,LMSUMC,KDGOF,
     *  KNODE,NMDOF,IC,KDGOFL,KBSYS0,KBSYS1,NMAT,JMATL,
     *  KDFRC,NCOOR2,MULNO,ILC1,ELEMSIZE,
     *  LSYS1,LSYS0,NODI,XMAST,XPROJ,DIRECT,DE,SHAPV,UE,GAP,
     *  FRCCOF,XMAT,COORL1,NODVAR,IDL,STAND,AREAP0)
        IMPLICIT real*8 (A-H,O-Z)
        DIMENSION LSYS1(KBSYS1,LMSUMC),LSYS0(KBSYS0,LMSUMC),NODI(4),
     *  XMAST(3),XPROJ(3),DIRECT(9,LMSUMC),DE(3,3),SHAPV(4,LMSUMC),
     *  UE(4),GAP(NMDOF,LMSUMC),FRCCOF(KDFRC,LMSUMC),XMAT(NMAT),
     *  NODVAR(KDGOF,KNODE),IDLM(3),IDLS(3),IDL(NMDOF,LMSUMC),
     *  COORL1(NCOOR2,LMSUMC),STAND(LMSUMC),AREAP0(LMSUMC)
6       FORMAT(2X,8I6)
7       FORMAT(2X,6E12.5)
 
        DO II=1,2
          LSYS1(II,ILC1)=LSYS0(II,IC)
        ENDDO
        LSYS1(3,ILC1)=NBLK
        LSYS1(4,ILC1)=NNE
        DO JNOD=1,NNE
          LSYS1(JNOD+4,ILC1)=NODI(JNOD)
        ENDDO
CCCC        IF (LSYS0(3,IC).GT.0) LSYS1(9,ILC1)= 1
CCCC        IF (LSYS0(3,IC).LT.0) LSYS1(9,ILC1)=-1
        LSYS1(9,ILC1)=LSYS0(3,IC)
        DO NR=1,NCOOR
          COORL1(NR,ILC1)=XMAST(NR)
          COORL1(NR+NCOOR,ILC1)=XPROJ(NR)
        ENDDO
C------ THE ARRAY DIRECT IS TO SAVE THE LOCAL COORDINATES.
        IJ=0
        DO JJ=1,NCOOR
          DO II=1,NCOOR
            IJ=IJ+1
            DIRECT(IJ,ILC1)=DE(II,JJ)
          ENDDO
        ENDDO
C------ THE ARRAY SHAPV IS TO SAVE THE SHAP FUNCTION VALUE.
        DO JNOD=1,NNE
          SHAPV(JNOD,ILC1)=UE(JNOD)
        ENDDO
C------ THE ARRAY GAP IS TO SAVE THE GAP .
        DO NDF=1,NMDOF
          GAP(NDF,ILC1)=0.0
          DO NR=1,NCOOR
          IF (NDF.EQ.1)
     *    GAP(NDF,ILC1)=GAP(NDF,ILC1)+(XMAST(NR)-XPROJ(NR))*DE(NR,NDF)
          ENDDO
          IF (NDF.GT.1) GAP(NDF,ILC1)=-GAP(NDF,ILC1)
        ENDDO
C------ THE ARRAY FRCCOF IS TO SAVE THE FRICTION COFFI.
        MULNO1=ABS(MULNO)
        FRCCOF(1,ILC1)=XMAT((MULNO1-1)*JMATL+1)
        FRCCOF(2,ILC1)=XMAT((MULNO1-1)*JMATL+2)
        FRCCOF(3,ILC1)=AREAP0(IC)
        STAND(ILC1)=ELEMSIZE
C---------------------------------------------------------
        KEYIDL=1
        IF (LSYS1(1,ILC1).EQ.LSYS1(3,ILC1)) THEN
          INODM=LSYS1(2,ILC1)
          DO JNOD=1,NNE
            INODS=LSYS1(JNOD+4,ILC1)
            KEYIDLINOD=1
            DO KDF=1,KDGOFL
              NVMASTER=NODVAR(KDF,INODM)
              NVSLAVE =NODVAR(KDF,INODS)
              IF (NVMASTER.GT.0.AND.NVSLAVE.GT.0.AND.NVMASTER.EQ.
     *            NVSLAVE) KEYIDLINOD=0
CC              WRITE(*,*) 'INODM,INODS,KDF,NVMASTER,NVSLAVE,KEYIDLINOD='
CC              WRITE(*,6)  INODM,INODS,KDF,NVMASTER,NVSLAVE,KEYIDLINOD
            ENDDO
            IF (KEYIDLINOD.EQ.0) KEYIDL=0
          ENDDO
        ENDIF
CC        WRITE(*,*) 'KEYIDL *******',KEYIDL
        IF (KEYIDL.EQ.0) THEN
          DO NDF=1,NMDOF
            IDL(NDF,ILC1)=-1
          ENDDO
          GOTO 100
        ENDIF
        CALL GETIDL(LMSUMC,KDGOF,KNODE,NNE,NCOOR,
     *  NMDOF,IC,KDGOFL,KBSYS0,
     *  NODVAR,LSYS0,IDLM,IDLS,NODI,DE,UE)
        DO NDF=1,NMDOF
          IDL(NDF,ILC1)=MAX(IDLM(NDF),IDLS(NDF))
        ENDDO
        IF (IDL(1,ILC1).EQ.-1) THEN
          DO NDF=2,NMDOF
            IDL(NDF,ILC1)=-1
          ENDDO
        ENDIF
100     CONTINUE
CC        WRITE(*,*) 'IC,ILC1 ==',IC,ILC1
CC        WRITE(*,*) 'IDLM,IDLS,IDL ==='
CC        WRITE(*,6) (IDLM(NDF),NDF=1,NMDOF),
CC     *  (IDLS(NDF),NDF=1,NMDOF),IDL(1,ILC1)
         RETURN
         END
 
        SUBROUTINE POSTVALUECONTACT(JLC1,ILC1,KBSYS1,LMSUMC,NCOOR,NNE,
     *  KDFRC,NCOOR2,NMDOF,LGNE,
     *  LSAME,LSYS1,COORL1,DIRECT,SHAPV,GAP,FRCCOF,STAND,IDL)
        IMPLICIT real*8 (A-H,O-Z)
        DIMENSION LSYS1(KBSYS1,LMSUMC),COORL1(NCOOR2,LMSUMC),
     *  DIRECT(9,LMSUMC),SHAPV(4,LMSUMC),GAP(NMDOF,LMSUMC),
     *  FRCCOF(KDFRC,LMSUMC),STAND(LMSUMC),IDL(NMDOF,LMSUMC),
     *  LSAME(LMSUMC)
        DIMENSION LSAMI(100),DE(3,3),XNORM1(3),XNORM2(3),POINT(3)
6       FORMAT(2X,10I6)
7       FORMAT(2X,6E12.5)
 
        IF (LGNE.eq.0) THEN
          DO IC1=1,ILC1
            IF (LSYS1(9,IC1).GT. 1) LSYS1(9,IC1)= 1
            IF (LSYS1(9,IC1).LT.-1) LSYS1(9,IC1)=-1
          ENDDO
          RETURN
        ENDIF
c----------------------------------------------------------------
C------ TO FIND OUT THE MULTIPLIERS WHICH CONNECT THE SAME POINTS.
C------ THE NUMBERS OF THE SAME MULTIPLIERS ARE SAVED AS lsam.
C------ AS FOR THE ARRAY LSAME, ITS ELEMENTS ARE ZERO IF THEY HAVE
C------ NO SAME MULTIPLIERS. IF THE ELEMENTS OF THE ARRAY LSAME ARE
C------ NEGATIVE, THEY DENOTE THAT THESE MULTIPLIERS BELONG TO THE
C------ SAME MULTIPLIER FACE AND THE SAME POINT.
c----------------------------------------------------------------
 
        LSAM=0
        DO 100 IC1=1,ILC1
          KEY=0
          DO 200 JC1=1,IC1-1
            SIZEMIN=MIN(STAND(IC1),STAND(JC1))
            SIZEMIN=SIZEMIN*0.005
            KEYDIST=1
            ERRDISTM=0.0
            ERRDISTS=0.0
            DO NR=1,NCOOR
              ERRDISTM=ERRDISTM+(COORL1(NR,IC1)-COORL1(NR,JC1))**2.
              ERRDISTS=ERRDISTS+
     *                (COORL1(NR+NCOOR,IC1)-COORL1(NR+NCOOR,JC1))**2.
            ENDDO
            ERRDISTM=SQRT(ERRDISTM)
            ERRDISTS=SQRT(ERRDISTS)
            IF (ERRDISTM.LE.SIZEMIN.AND.ERRDISTS.LE.SIZEMIN) KEYDIST=0
            IF (KEYDIST.EQ.1) THEN
              ERRDISTM=0.0
              ERRDISTS=0.0
              DO NR=1,NCOOR
                ERRDISTM=ERRDISTM+
     *                   (COORL1(NR,IC1)-COORL1(NR+NCOOR,JC1))**2.
                ERRDISTS=ERRDISTS+
     *                   (COORL1(NR+NCOOR,IC1)-COORL1(NR,JC1))**2.
              ENDDO
              ERRDISTM=SQRT(ERRDISTM)
              ERRDISTS=SQRT(ERRDISTS)
            IF (ERRDISTM.LE.SIZEMIN.AND.ERRDISTS.LE.SIZEMIN) KEYDIST=0
            ENDIF
            IF (KEYDIST.EQ.0) THEN
              KEY=1
              GOTO 210
            ENDIF
200       CONTINUE
210       CONTINUE
          IF (KEY.EQ.1) THEN
            IF (LSAME(JC1).NE.0) THEN
              LSAME(IC1)=LSAME(JC1)
            ELSE
              KSAMEP1=LSYS1(2,IC1)
              KSAMEP2=LSYS1(2,JC1)
              MULN0=LSYS1(9,IC1)
              MULN1=LSYS1(9,JC1)
              IF (KSAMEP1.EQ.KSAMEP2.AND.MULN0.EQ.MULN1) THEN
                LSAM=LSAM+1
                LSAME(JC1)=-LSAM
                LSAME(IC1)=-LSAM
              ELSE
                LSAM=LSAM+1
                LSAME(JC1)=LSAM
                LSAME(IC1)=LSAM
              ENDIF
            ENDIF
          ENDIF
100     CONTINUE
CC        WRITE(*,*) 'LSAM,ILC1,LSAME ==',LSAM
CC        WRITE(*,6) (LSAME(IC),IC=1,ILC1)
C-------------------------------------------------------------------
C------ DEAL WITH THE SAME MULTIPLIERS BY DIFFIRENT METHODS.
C------ FOR EQ. CASE, OMIT THE FIRST MULTIPLIERS OF THE SAME
C------ MULTIPLIERS; FOR INEQ. CASE, HERE ONLY THE MULTIPLIERS OF
C------ THE SAME POINTS AND THAT WHOSE SUM IS TWO ARE DEALED WITH.
C------ THE DEALING METHOD IS TO AVERAGE THE NORMAL DIRECTION.
C-------------------------------------------------------------------
        IF (LGNE.eq.1) THEN
          DO I=1,LSAM
            KEY=0
            DO IC1=1,ILC1
              IF (LSAME(IC1).EQ.I.AND.KEY.EQ.1) LSAME(IC1)=0
              IF (LSAME(IC1).EQ.I.AND.KEY.EQ.0) KEY=1
            ENDDO
          ENDDO
        ELSE
CC          DO 300 I=1,LSAM
CC            LSUM=0
CC            DO IC1=1,ILC1
CC              ISAM=LSAME(IC1)
CC              IF (ISAM.LT.0) ISAM=-ISAM
CC              IF (ISAM.EQ.I) THEN
CC                LSUM=LSUM+1
CC                LSAMI(LSUM)=IC1
CC              ENDIF
CC            ENDDO
CC            IF (LSUM.GT.0) THEN
CC              IC1=LSAMI(1)
CC              ISAM=LSAME(IC1)
CC              IF (ISAM.LT.0.OR.LSUM.EQ.2) THEN
CC                DO 310 J=2,LSUM
CC                  JC1=LSAMI(J)
CC                  DO NR=1,NCOOR
CC                    XNORM1(NR)=DIRECT(NR,IC1)
CC                    XNORM2(NR)=DIRECT(NR,JC1)
CC                    POINT(NR)=COORL1(NR+NCOOR,IC1)
CC                  ENDDO
CC                  PROD=0.0
CC                  DO NR=1,NCOOR
CC                    PROD=PROD+XNORM1(NR)*XNORM2(NR)
CC                  ENDDO
CC                  IF (PROD.LT.0.0) THEN
CC                    DO NR=1,NCOOR
CC                      XNORM2(NR)=-XNORM2(NR)
CC                    ENDDO
CC                  ENDIF
CC                  CALL SECTIONANGLE(NCOOR,XNORM1,XNORM2,DE,POINT)
CC                  IJ=0
CC                  DO JJ=1,NCOOR
CC                    DO II=1,NCOOR
CC                      IJ=IJ+1
CC                      DIRECT(IJ,IC1)=DE(II,JJ)
CC                    ENDDO
CC                  ENDDO
CC310             CONTINUE
CC              ENDIF
CC              LSAME(IC1)=0
CC              FRCCOF(1,IC1)=1.E+10
CC            ENDIF
CC300       CONTINUE
          DO IC=1,ILC1
            LSAME(IC)=0
          ENDDO
        ENDIF
CC        WRITE(*,*) 'LSAM,ILC1,LSAME AFTER DEALING ==',LSAM
CC        WRITE(*,6) (LSAME(IC),IC=1,ILC1)
        DO IC1=1,ILC1
          IF (LSYS1(9,IC1).GT. 1) LSYS1(9,IC1)= 1
          IF (LSYS1(9,IC1).LT.-1) LSYS1(9,IC1)=-1
        ENDDO
        JLC1=0
        DO 400 IC1=1,ILC1
          IF (LSAME(IC1).NE.0) GOTO 400
          JLC1=JLC1+1
          DO J=1,9
            LSYS1(J,JLC1)=LSYS1(J,IC1)
          ENDDO
          DO NR=1,NCOOR
            COORL1(NR,JLC1)=COORL1(NR,IC1)
            COORL1(NR+NCOOR,JLC1)=COORL1(NR+NCOOR,IC1)
          ENDDO
C-------- THE ARRAY DIRECT IS TO SAVE THE LOCAL COORDINATES.
          IJ=0
          DO JJ=1,NCOOR
            DO II=1,NCOOR
              IJ=IJ+1
              DIRECT(IJ,JLC1)=DIRECT(IJ,IC1)
            ENDDO
          ENDDO
C-------- THE ARRAY SHAPV IS TO SAVE THE SHAP FUNCTION VALUE.
          DO JNOD=1,NNE
            SHAPV(JNOD,JLC1)=SHAPV(JNOD,IC1)
          ENDDO
C-------- THE ARRAY GAP IS TO SAVE THE GAP .
          DO NDF=1,NMDOF
            GAP(NDF,JLC1)=GAP(NDF,IC1)
          ENDDO
C-------- THE ARRAY FRCCOF IS TO SAVE THE FRICTION COFFI.
          DO NDF=1,KDFRC
            FRCCOF(NDF,JLC1)=FRCCOF(NDF,IC1)
          ENDDO
          STAND(JLC1)=STAND(IC1)
          DO NDF=1,NMDOF
            IDL(NDF,JLC1)=IDL(NDF,IC1)
          ENDDO
400     CONTINUE
        RETURN
        END
 
        SUBROUTINE SECTIONANGLE(NCOOR,XNORM1,XNORM2,DE,POINT)
        IMPLICIT real*8 (A-H,O-Z)
        DIMENSION DE(3,3),XNORM1(3),XNORM2(3),XD(3),YD(3),ZD(3),
     *  POINT1(3),POINT2(3),POINT3(3),POINT(3)
6       FORMAT(2X,8I8)
7       FORMAT(2X,6E12.5)
 
        DO I=1,3
          DO J=1,3
            DE(I,J)=0.0
          ENDDO
        ENDDO
        DO NR=1,NCOOR
          POINT1(NR)=POINT(NR)+XNORM1(NR)
          POINT2(NR)=POINT(NR)+XNORM2(NR)
          POINT3(NR)=(POINT1(NR)+POINT2(NR))/2.
        ENDDO
        DO NR=1,NCOOR
          ZD(NR)=POINT3(NR)-POINT(NR)
        ENDDO
        ZLEN=0.0
        DO NR=1,NCOOR
          ZLEN=ZLEN+ZD(NR)**2.
        ENDDO
        ZLEN=SQRT(ZLEN)
        DO NR=1,NCOOR
          ZD(NR)=ZD(NR)/ZLEN
        ENDDO
        IF (NCOOR.EQ.2) THEN
          DE(1,1)= ZD(1)
          DE(2,1)= ZD(2)
          DE(1,2)=-ZD(2)
          DE(2,2)= ZD(1)
        ENDIF
        IF (NCOOR.EQ.3) THEN
          DO NR=1,NCOOR
            DE(NR,1)=ZD(NR)
            ZD(NR)=XNORM1(NR)
            XD(NR)=XNORM2(NR)
          ENDDO
          CALL GETE3(ZD,XD,YD)
          DO NR=1,NCOOR
            ZD(NR)=DE(NR,1)
          ENDDO
          CALL GETE3(YD,ZD,XD)
          XLEN=0.0
          YLEN=0.0
          DO NR=1,NCOOR
            XLEN=XLEN+XD(NR)**2.
            YLEN=YLEN+YD(NR)**2.
          ENDDO
          XLEN=SQRT(XLEN)
          YLEN=SQRT(YLEN)
          DO NR=1,NCOOR
            XD(NR)=XD(NR)/XLEN
            YD(NR)=YD(NR)/YLEN
          ENDDO
          DO NR=1,NCOOR
            DE(NR,2)=XD(NR)
            DE(NR,3)=YD(NR)
          ENDDO
        ENDIF
        RETURN
        END
 
C -------- 4 : GETE3
        subroutine gete3(e1,e2,e3)
        implicit real*8(a-h,o-z)
        dimension e1(3),e2(3),e3(3)
        e3(1) = e1(2)*e2(3)-e1(3)*e2(2)
        e3(2) = e1(3)*e2(1)-e1(1)*e2(3)
        e3(3) = e1(1)*e2(2)-e1(2)*e2(1)
        return
        end
 
        subroutine getec(e1,e2,e3)
        implicit real*8(a-h,o-z)
        dimension e1(3),e2(3),e3(3)
        e = 1.e-10
        de1 = 0.0
        de2 = 0.0
        de3 = 0.0
        do 300 i=1,3
        de1 = de1 + e1(i)*e1(i)
        de2 = de2 + e2(i)*e2(i)
        de3 = de3 + e3(i)*e3(i)
300     continue
        if (de1.lt.e .or. de2.lt.e .or. de3.lt.e) then
        write(*,*) ' Improperty case in getec =='
        write(*,*) 'de = ',de1,de2,de3
        stop 2222
        endif
        do 400 i=1,3
        e1(i) = e1(i)/sqrt(de1)
        e2(i) = e2(i)/sqrt(de2)
        e3(i) = e3(i)/sqrt(de3)
400     continue
        det3 = e1(1)*e2(2)*e3(3)+e2(1)*e3(2)*e1(3)+e1(2)*e2(3)*e3(1)
     &       - e1(3)*e2(2)*e3(1)-e2(3)*e3(2)*e1(1)-e1(2)*e2(1)*e3(3)
        if (det3.lt.0.0) then
        do 500 i=1,3
        d = e1(i)
        e1(i) = e2(i)
        e2(i) = d
500     continue
        endif
        return
        end

	
       subroutine ch(char1,ext,name)
       character char1*8,ext*3,name*12
       do 100 i=1,12
       name(i:i)=' '
100    continue
       lef=0
       do 200 i =1,8
       if(char1(i:i).eq.' ') goto 10
       lef=lef+1
       name(i:i)=char1(i:i)
200    continue
10     lef = lef+1
       name(lef:lef)='.'
       do 300 i = 1,3
       j = lef+i
       name(j:j) = ext(i:i)
300    continue
c       write(*,*) lef,',',name
       return
       end



       subroutine dualnode(num,nnode,kelem,knode,node,idelem,inode,
     + idnode,ndual,nbnodes,knode3)
       implicit real*8(a-h,o-z)
       logical flg,iflg
       dimension node(kelem),inode(knode),idnode(knode),nbnodes(knode3)
       dimension idelem(num)
       dimension mul(20)
c
c      find the boundary nodes which shared by more than one subdomain
c
       do i=1,knode
       inode(i) = -100
       idnode(i) = 1
       end do
       do i=1,knode3
       nbnodes(i) = -200
       end do
c
       nb=0
       nne = nnode-1
       do i=1,num
       do j=1,nne
       noij=node((i-1)*nnode+j)
       jb = idelem(i)
       ib = inode(noij)
       if(ib.eq.-100) then
       inode(noij)=jb
       else
       if(ib.ne.jb) then
       idnode(noij)=idnode(noij) + 1
       nb=nb+1
       nbnodes((nb-1)*6+1) = noij
       nbnodes((nb-1)*6+2) = ib
       nbnodes((nb-1)*6+4) = jb
       nbnodes((nb-1)*6+6) = 1
       end if
       end if
       end do
       end do
c
       do i=1,knode
       if(idnode(i).ge.3) then
       noyes=0
       nmul=0
       do j = 1,nb
       if(nbnodes((j-1)*6+1).eq.i) then
       noyes = 0
       do k=1,nmul
       if(mul(k).eq.nbnodes((j-1)*6+4)) noyes=1
       end do
       if(noyes.eq.0) then
       nmul=nmul+1
       mul(nmul)=nbnodes((j-1)*6+4)
       end if
       if(noyes.eq.1) then
       nbnodes(j*6) = 0
       end if
       end if             !if(nbnodes((j-1)*6+1).eq.i)
       end do             !do j=1,nb
       end if             !end if(idnode(i).ge.3)
2000   continue
       end do             !i=1,knode
c
       ndual = 0
       do i=1,nb
       if(nbnodes((i-1)*6+2).eq.nbnodes((i-1)*6+4)) then
       write(*,*) 'fatal error when get dual nodal parts'
       stop 1111
       end if
       if(nbnodes(i*6).eq.1) then
       ndual = ndual + 1
       nbnodes((ndual-1)*6+1) = nbnodes((i-1)*6+1)
       nbnodes((ndual-1)*6+2) = nbnodes((i-1)*6+2)
       nbnodes((ndual-1)*6+3) = nbnodes((i-1)*6+3)
       nbnodes((ndual-1)*6+4) = nbnodes((i-1)*6+4)
       nbnodes((ndual-1)*6+5) = nbnodes((i-1)*6+5)
       nbnodes((ndual-1)*6+6) = nbnodes((i-1)*6+6)
       end if
       end do

       write(*,*) ndual,'  nodal parts found..........'

c
       open(80,file='multipliers',form='formatted',status='unknown')
       write(80,*) ndual,6
       do i=1,ndual
       write(80,1000)  (nbnodes((i-1)*6+j),j=1,6)
       end do
       close(80)
1000   format(6i10)
c
       return
       end

