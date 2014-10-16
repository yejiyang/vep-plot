      SUBROUTINE GETDFS(PE,PV,YITA,DT,STRESS,STRAIN,D,F)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION D(6,6),F(6),STRESS(6),STRAIN(6),DE(6,6),
     *DE_CONTRA(6,6),DV_CONTRA(6,6),D_CONTRA(6,6),FTEMP(6)
      FACT=PE/(1.d0+PV)/(1.d0-2.d0*PV)
      DO I=1,6
      F(I)=0.d0
      FTEMP(I)=0.d0
      DO J=1,6
      DE(I,J)=0.d0
      DV_CONTRA(I,J)=0.d0
      ENDDO
      ENDDO
      DE(1,1)=FACT*(1.d0-PV)
      DE(1,2)=FACT*PV
      DE(1,3)=FACT*PV
      DE(2,1)=FACT*PV
      DE(2,2)=FACT*(1.d0-PV)
      DE(2,3)=FACT*PV
      DE(3,1)=FACT*PV
      DE(3,2)=FACT*PV
      DE(3,3)=FACT*(1.d0-PV)
      DE(4,4)=FACT*(0.5d0-PV)
      DE(5,5)=FACT*(0.5d0-PV)
      DE(6,6)=FACT*(0.5d0-PV)
      CALL GET_CONTRA(DE,DE_CONTRA)
      DV_CONTRA(1,1)=1.d0/3.d0/YITA
      DV_CONTRA(1,2)=-1.d0/6.d0/YITA
      DV_CONTRA(1,3)=-1.d0/6.d0/YITA
      DV_CONTRA(2,1)=-1.d0/6.d0/YITA
      DV_CONTRA(2,2)=1.d0/3.d0/YITA
      DV_CONTRA(2,3)=-1.d0/6.d0/YITA
      DV_CONTRA(3,1)=-1.d0/6.d0/YITA
      DV_CONTRA(3,2)=-1.d0/6.d0/YITA
      DV_CONTRA(3,3)=1.d0/3.d0/YITA
      DV_CONTRA(4,4)=1.d0/YITA
      DV_CONTRA(5,5)=1.d0/YITA
      DV_CONTRA(6,6)=1.d0/YITA
      DO I=1,6
      DO J=1,6
      D_CONTRA(I,J)=DE_CONTRA(I,J)/DT+DV_CONTRA(I,J)
      ENDDO
      ENDDO
      CALL GET_CONTRA(D_CONTRA,D)
      DO I=1,6
      DO J=1,6
      D(I,J)=D(I,J)/DT
      ENDDO
      ENDDO
      DO I=1,6
      DO J=1,6
      FTEMP(I)=FTEMP(I)+DE_CONTRA(I,J)*STRESS(J)
      ENDDO
c      FTEMP(I)=FTEMP(I)+STRAIN(I)
      ENDDO
      DO I=1,6
      DO J=1,6
      F(I)=F(I)+D(I,J)*FTEMP(J)
      ENDDO
      ENDDO
      RETURN
      END
 
      SUBROUTINE GET_CONTRA(D1,D2)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION D1(6,6),D2(6,6),DTEMP(6,6)
      DO I=1,6
      DO J=1,6
      D2(I,J)=0.d0
      DTEMP(I,J)=D1(I,J)
      ENDDO
      D2(I,I)=1.d0
      ENDDO
      DO I=1,6
      C=DTEMP(I,I)
      DO K=1,6
      DTEMP(I,K)=DTEMP(I,K)/C
      D2(I,K)=D2(I,K)/C
      ENDDO
      DO J=I+1,6
      C=DTEMP(J,I)
      DO K=1,6
      DTEMP(J,K)=DTEMP(J,K)-DTEMP(I,K)*C
      D2(J,K)=D2(J,K)-D2(I,K)*C
      ENDDO
      ENDDO
      ENDDO
      DO I=6,2,-1
      DO J=I-1,1,-1
      C=DTEMP(J,I)
      DO K=1,6
      DTEMP(J,K)=DTEMP(J,K)-DTEMP(I,K)*C
      D2(J,K)=D2(J,K)-D2(I,K)*C
      ENDDO
      ENDDO
      ENDDO
      RETURN
      END
 
      SUBROUTINE GETDE(PE,PV,STRESS,STRAIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION STRESS(6),STRAIN(6),DE(6,6),
     *DE_CONTRA(6,6),FTEMP(6)
      FACT=PE/(1.d0+PV)/(1.d0-2.d0*PV)
      DO I=1,6
      FTEMP(I)=0.d0
      DO J=1,6
      DE(I,J)=0.d0
      ENDDO
      ENDDO
      DE(1,1)=FACT*(1.d0-PV)
      DE(1,2)=FACT*PV
      DE(1,3)=FACT*PV
      DE(2,1)=FACT*PV
      DE(2,2)=FACT*(1.d0-PV)
      DE(2,3)=FACT*PV
      DE(3,1)=FACT*PV
      DE(3,2)=FACT*PV
      DE(3,3)=FACT*(1.d0-PV)
      DE(4,4)=FACT*(0.5d0-PV)
      DE(5,5)=FACT*(0.5d0-PV)
      DE(6,6)=FACT*(0.5d0-PV)
      CALL GET_CONTRA(DE,DE_CONTRA)
      DO I=1,6
      DO J=1,6
      FTEMP(I)=FTEMP(I)+DE_CONTRA(I,J)*STRESS(J)
      ENDDO
      ENDDO
      DO I=1,6
      STRAIN(I)=FTEMP(I)
      ENDDO
      RETURN
      END
