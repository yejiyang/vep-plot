      SUBROUTINE GETSTR(X,U,DE,DV,ddv,D,DP,KKK,F1,FSUB)
      implicit real*8 (a-h,o-z)
        DIMENSION H(7),U(7),DE(6),DV(6),D(6,6),DP(6,6),V(7),DU(6),X(*)
        DIMENSION P(21),DF(7),ddv(6),vc1(7),va(7),vb1(7),vb2(7)
        dimension D_CONTRA(6,6),dde(6)
      EXTERNAL F
      EXTERNAL FSUB
C     write(*,*) 'U =',U
C     write(*,*) 'DE =',DE
        N = 7
        NDF = 6
c        E = 1.d5
      E = x(2)*1.d-5
      call get_contra(D,D_CONTRA)
      do i=1,NDF
        dde(i)=0.d0
        do j=1,NDF
          dde(i)=dde(i)+D_CONTRA(i,j)*ddv(j)
        enddo
      enddo
      DO 40 I=1,NDF
      C = 0.0+ddv(i)
      DO 20 J=1,NDF
      C = C+D(I,J)*DE(J)
20    CONTINUE
      DU(I) = C
40    CONTINUE
c       de-delta strain compotent  du-elastic stress delta u-old stress
c        v-new stress   dv-delta stress
      DO 100 I=1,N
      V(I) = U(I)
      IF (I.LE.NDF) V(I) = V(I)+DU(I)
100   CONTINUE
      F1 = FSUB(X,V)
      vc1(N)=V(N)
      va(N)=V(N)
      vb1(N)=V(N)
      vb2(N)=V(N)
c       WRITE(*,*) 'U & X =','f1'
c       WRITE(*,*) (U(I),I=1,4),(X(I),I=1,2),f1
      IF (F1.LT.-E) THEN
      DO 150 I=1,NDF
      DV(I) = DU(I)
150   CONTINUE
c       write(*,*) 'elastic deformation'
      RETURN
      ENDIF
      F0 = FSUB(X,U)
      R = 0.0
      IF (F0.LT.-E*1.1) then
        R = F0/(F0-F1)
        beta0=0.d0
        beta1=R
        fc0=f0
        i=0
123     continue
        i=i+1  
        do j=1,NDF
          vc1(j)=u(j)+beta1*du(j)
        enddo
        fc1=fsub(x,vc1)
        beta2=beta1-fc1*(beta0-beta1)/(fc0-fc1)
        beta0=beta1
        beta1=beta2
        fc0=fc1
        if (i.le.5.and.dabs(beta1-beta0).ge.1.d-8) then
        goto 123
        endif
        R=beta2
      endif

c      if (R.lt.0.d0) R=0.d0
c      if (R.gt.1.d0) R=1.d0

c      write(*,*) R,beta2
c      WRITE(*,*) 'F0,F1,R ------------ ',F0,F1,R
c      IF(F1.GT.100.) THEN
c      WRITE(*,*) 'V================'
c        WRITE(*,*) (V(I),I=1,4)
c      WRITE(*,*) 'DU==============='
c        WRITE(*,*) (DU(I),I=1,3)
c      END IF
      DO 200 I=1,NDF
      V(I) = U(I)+R*DU(I)
200   CONTINUE
C       WRITE(*,*) 'V =',(V(I),I=1,4
      CALL GETDP(X,V,D,DP,KKK,FSUB)
      DO 400 I=1,NDF
      Va(I) = U(I)+DU(I)
      C = 0.0
      DO 300 J=1,NDF
      C = C+DP(I,J)*(DE(J)+dde(j))
300   CONTINUE
      DV(I) = C
400   CONTINUE
      DO 430 I=1,NDF
      IF (F0.GT.-E*1.1) THEN
      Vb1(I) = Va(I)-DV(I)
      vb2(i)=vb1(i)
      ELSE
      Vb1(I) = Va(I) - (1.-R)*DV(I)
      vb2(i)=vb1(i)
      ENDIF
430   CONTINUE
      do 234 ijk=1,5
      CALL GETDP(X,Vb2,D,DP,KKK,FSUB)
      DO 401 I=1,NDF
      C = 0.0
      DO 301 J=1,NDF
      C = C+DP(I,J)*(DE(J)+dde(j))
301   CONTINUE
      DV(I) = C
401   CONTINUE
      DO 431 I=1,NDF
      IF (F0.GT.-E*1.1) THEN
      Vb2(I) = Va(I)-DV(I)
      ELSE
      Vb2(I) = Va(I) - (1.-R)*DV(I)
      ENDIF
431   CONTINUE
      do i=1,6
        v(i)=(vb1(i)+vb2(i))/2
        vb2(i)=v(i)
      enddo
234   continue
      CALL GETDF(X,V,DF,FSUB)
      DO 450 I=1,NDF
      DV(I) =DF(I)
450   CONTINUE
      R = 1.0
      DD=0.1
      nliq=0
470   DO 500 I=1,N
      P(I) = V(I)
      IF (I.LT.N) P(I+N) = DV(I)
      P(I+N*2) = X(I)
500   CONTINUE
      CALL ROOT(R,DD,E,p,F,fsub,k)
      nliq=nliq+1
      if(k.eq.1.and.nliq.le.20) then
      do 550 i=1,ndf
      v(i)=v(i)-(1.-r)*dv(i)
550   continue
      call getdf(x,v,df,fsub)
      do 570 i=1,ndf
      dv(i)=df(i)
570   continue
      r = 1.0
      goto 470
      end if
      DO 600 I=1,NDF
      V(I) = V(I) - (1.-R)*DV(I)
      DV(I) = V(I)-U(I)
600   CONTINUE
c      do i=1,6
c      va(i)=(v(i)+va(i))/2
c      enddo
      CALL GETDP(X,V,D,DP,KKK,FSUB)    
      F1=FSUB(X,V)
C       WRITE(*,*) 'V +',(V(I),I=1,6)
      if (PRAG.gt.1.d5) then
       WRITE(*,*) 'R,PRAG =',R,F1,nliq
      endif
      if (f1.lt.f0) f1=f0
      RETURN
      END
 
      real*8 FUNCTION F(R,P,FSUB)
      implicit real*8 (a-h,o-z)
      DIMENSION X(7),U(7),P(*)
      EXTERNAL FSUB
        DO 100 I=1,7
        U(I) = P(I)-P(I+7)*(1.-R)
100     CONTINUE
        U(7) = P(7)
      DO 200 I=1,7
        X(I) = P(I+14)
200     CONTINUE
      F=FSUB(X,U)
c        WRITE(*,*) 'U -----',(U(I),I=1,3)
      RETURN
      END
 
c      SUBROUTINE GETK(DE,DV,D,DK)
c      implicit real*8 (a-h,o-z)
c        DIMENSION DE(3),DV(3),DEP(3),D(3,3)
c        DET =+D(1,1)*D(2,2)-D(2,1)*D(1,2)
c      dep(1)=de(1)-(dv(1)*d(2,2)-dv(2)*d(1,2))/det
c      dep(2)=de(2)-(-dv(1)*d(2,1)+dv(2)*d(1,1))/det
c      dk=dep(1)+dep(2)
c      END
 
      subroutine root(x1,d,e,p,f,fsub,k)
      implicit real*8 (a-h,o-z)
      dimension p(*)
      EXTERNAL FSUB
      EXTERNAL F
c      implicit double precision (a-h,o-z)
C      write (*,1000)
1000  format (1x,12hinput x1 d e)
c     k=0
1099  continue
      k=0
      f1 = f(x1,p,fsub)
      if (abs(f1).lt.e) goto 999
      if (f1.gt.0.0) then
      x2 = x1 - d
      else
      x2 = x1 + d
      endif
      f2 = f(x2,p,fsub)
c        WRITE(*,*) 'x1,d,x2,e ='
c      WRITE (*,'(1x,4f16.6)') x1,d,x2,e
c      write(*,*) 'f1,f2 =',f1,f2
c      write(*,*) 'p==',p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8)
c       write(*,*) 'p==',p(11),p(12),p(13),p(14)
      n12 = 0
10    if (f1*f2 .le. 0.0) go to 30
      if (f1*f1 .lt. f2*f2) go to 20
      x1 = x2 + 1.1*(x2-x1)
      f1 = f(x1,p,fsub)
c      write(*,*) 'x1,f1 =',x1,f1
      if (n12.ge.1) then
      x=x2
      go to 85
      end if
      n12 = 1
      goto 10
20    x2 = x1 + 1.1*(x1-x2)
      f2 = f(x2,p,fsub)
c      write(*,*) 'x2,f2 =',x2,f2
      if (n12.le.-1) then
      x=x1
      go to 85
      end if
      n12 = -1
      go to 10
30    if (x1 .ge. x2) go to 40
      xm1 = x1
      xm2 = x2
      go to 50
40    xm1 = x2
      xm2 = x1
c       write(*,*) 'x1,x2,f1,f2 =',x1,x2,f1,f2
50    continue
      if (f2.gt.f1) then
      x3 = (x1*f2-x2*f1)/((f2-f1)+1.d-20)
      else
      x3 = (x1*f2-x2*f1)/((f2-f1)-1.d-20)
      endif
      if (x3 .ge. xm1) go to 60
      x3 = xm1
60    if (x3 .le. xm2) go to 70
      x3 = xm2
70    f3 = f(x3,p,fsub)
c       write(*,*) 'x1,x2,x3 =',x1,x2,x3
c       write(*,*) 'f1,f2,f3 =',f1,f2,f3
      if (abs(f3) .lt. e*0.1) go to 90
      if (abs(x3).gt.1e-6.and.abs((x3-x1)/x3).lt.1.0e-6) go to 90
      if (abs(x3).gt.1e-6.and.abs((x3-x2)/x3).lt.1.0e-6) go to 90
      if (f2*f3 .ge. 0.0) then
      f2 = f3
      x2 = x3
      else
      f1 = f3
      x1 = x3
      end if
      go to 50
85    continue
c      write (*,3000)
      if (dabs(x1-x2).lt.e*0.1) then   
      k=1
      x1=x
c       if(d.gt.1.e17)return
      return
      endif
      x1=x
c       d=d*1.1
c       write(*,*)'goto 1099'
       goto 1099
990    continue
c      write (*,2000) x1, f(x1,p,fsub)
c      write(*,*) 'p==',p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8)
c      write(*,*) 'p==',p(9),p(10),p(11),p(12)
      return
3000  format(5x,18hfail,nonmonotone f)
90    x1 = x3
c      write(*,*) 'nliq',nliq,f(x1,p,fsub)
999   continue
      if (dabs(x1).gt.1.d20) x1=-1.d20
c      write(*,*) x1, f(x1,p,fsub)
2000  format(1x,3hx =,f16.6,3x,3hf =,f16.6)
      return
      end
