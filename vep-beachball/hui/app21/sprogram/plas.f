      SUBROUTINE GETDP(X,U,D,DP,KKK,yield)
      implicit real*8 (a-h,o-z)
        DIMENSION U(7),D(6,6),DP(6,6),X(*),DF(7),DDF(6)
        dimension dg(7),ddg(6),xg(7)
      EXTERNAL YIELD
c        write(*,*)'(df(i),i=1,4)'
c        write(*,*)(df(i),i=1,4)
        N = 7
        NDF = 6
      do i=1,N
      xg(i)=x(i)
      enddo
      xg(1)=0.d0

      CALL GETDF(X,U,DF,yield)
      DO 200 I=1,NDF
      C = 0.0
      DO 100 J=1,NDF
      C = C+D(I,J)*DF(J)
100     CONTINUE
      DDF(I) = C
200     CONTINUE
      CALL GETDF(xg,U,dg,yield)
      DO I=1,NDF
      C = 0.0
      DO J=1,NDF
      C = C+D(I,J)*dg(J)
      enddo
      ddg(I) = C
      enddo
      A = 0.0
      AM = 0.0
      DO 400 I=1,NDF
      A = A+DF(I)*ddg(I)
cc    am hasn't been modified yet. it should be corrected later
cc      IF (KKK.EQ.1) AM = AM+U(I)*DF(I)
cc      IF (KKK.EQ.2 .AND. I.LE.3) AM = AM+DF(I)
cc      IF (KKK.EQ.3) AM = AM+DF(I)*DF(I)
400     CONTINUE
cc      IF (KKK.EQ.3) AM = SQRT(AM)
cc      A = A - AM*DF(N)
c       write(*,*) 'a=', a
c       if(abs(a).lt.1e-12)a=1.
      DO 600 I=1,NDF
      DO 500 J=1,NDF
c        write(*,*)'i,j,dp(i,j)'
c        write(*,*)i,j,dp(i,j)
        DP(I,J) = DDF(I)*ddg(J)/A
500     CONTINUE
600     CONTINUE
      RETURN
      END
 
      SUBROUTINE GETDF(X,U,DF,yield)
      implicit real*8 (a-h,o-z)
        DIMENSION H(7),U(7),X(*),DF(7),V(7)
      EXTERNAL YIELD
        N = 7
      DO 10 K = 1,N
c        H(K) = 0.1
c        H(K) = 10000.
        H(K) = x(2)*1.d-6
      V(K) = U(K)
10      CONTINUE
      DO 50 K = 1,N
      V(K) = U(K) + H(K)
      F2 = yield(X,V)
      V(K) = U(K) - H(K)
      F1 = yield(X,V)
      DF(K) = (F2-F1)/H(K)*0.5
      V(K) = U(K)
50      CONTINUE
      RETURN
      END
 
      real*8 FUNCTION PRAGER(X,U)
      implicit real*8 (a-h,o-z)
        DIMENSION U(7),S(6),x(4)
      tgo=x(1)
      q0=x(4)
cc      c=x(2)+q0*U(7)
c     get rid of strengthening
      c=x(2)
      PV= X(3)

        D = (U(1)+U(2)+U(3))/3.
          
      S(1) = U(1) - D
      S(2) = U(2) - D
      S(3) = U(3) - D
        S(4) = U(4)
        S(5) = U(5)
        S(6) = U(6)
        UK = U(7)
        SJ = (S(1)**2+S(2)**2+S(3)**2)/2.+S(4)**2
        SJ = SJ+S(5)**2+S(6)**2
         cosf=1./dsqrt(1.+tgo*tgo)
         sinf=dsqrt(1-cosf*cosf)
         alpha=2*sinf/(dsqrt(3.d0)*(3+sinf))
         ck=6*c*cosf/(dsqrt(3.d0)*(3+sinf))
c     PRAGER = DSQRT(SJ) + ALPHA*D*3 - CK

c     use the mean stress of last loading step
      PRAGER = DSQRT(SJ) + ALPHA*UK*3 - CK


      RETURN
      END
 
      real*8 FUNCTION MISES(X,U)
      implicit real*8 (a-h,o-z)
      DIMENSION U(7),S(6),x(4)
      tgo=x(1)
      q0=x(4)
      c=x(2)+q0*U(7)
      PV= X(3)
c        write(*,*)'tgo,c,v=',tgo,c,pv
        D = (U(1)+U(2)+U(3))/3.
      S(1) = U(1) - D
      S(2) = U(2) - D
      S(3) = U(3) - D
        S(4) = U(4)
        S(5) = U(5)
        S(6) = U(6)
        UK = U(7)
        SJ = (S(1)**2+S(2)**2+S(3)**2)/2.+S(4)**2
        SJ = SJ+S(5)**2+S(6)**2
      MISES = DSQRT(SJ) - C
      RETURN
      END

c      real*8 FUNCTION SEIM(X,U)
c      implicit real*8 (a-h,o-z)
c        DIMENSION U(4),X(4)
c        UK = U(4)
c        a=0.1
c      d=x(1)
c        c=x(2)
c        q1=x(3)
c        q2=x(4)
c        qx=U(1)*q1+U(3)*q2
c        qy=U(3)*q1+U(2)*q2
c        UN=qx*q1+qy*q2
c        UT=sqrt(qx*qx+qy*qy-UN*UN)
c        SEIM = SQRT(UT**2+a**2*c**2)+d*UN-c
c      RETURN
c      END
 
