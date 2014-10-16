 
      SUBROUTINE SMIT(K,N,M,R,Y,Z,NC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(K,M),Y(N,M),Z(K,K),R1(3),NC(1)
C .... K : number of dimension in global space
C .... N : number of dimension in local space
C .... M : number of nodes
C .... R : coordinate values of nodes in global coordinates
C .... Z : new coordinate axes for the global space
C ........ obtained by Schmidt's orthogonalization method
C .... Y : coordinate values of nodes in local coordinates
C ........ obtained by projecting R into Z
C .... Y = Z(t)*R i.e. Y(J,I) = Z(L,J)*R(L,I)  (L=1,K)
C ........                       (J=1,N & I=1,M)
      DO 1 I=1,K
      L=NC(1)
1     R1(I)=R(I,L)
      DO 3 J=1,M
      DO 3 I=1,K
3     R(I,J)=R(I,J)-R1(I)
c      WRITE(*,*) 'R ='
c      DO 1100 I=1,M
c      WRITE(*,*) (R(J,I),J=1,K)
c1100   CONTINUE
      DO 20 I=1,N
      I1=NC(I+1)
      DO 5 L=1,K
5     Z(L,I)=R(L,I1)
C     WRITE(*,*) 'Z ='
C     WRITE(*,*) (Z(J,I),J=1,K)
      DO 10 J=1,I-1
      C=0.0
      DO 7 L=1,K
7     C=C+R(L,I1)*Z(L,J)
      DO 9 L=1,K
9     Z(L,I)=Z(L,I)-C*Z(L,J)
10    CONTINUE
      C=0.0
      DO 15 L=1,K
15    C=C+Z(L,I)**2
      C=SQRT(C)
      DO 17 L=1,K
17    Z(L,I)=Z(L,I)/C
20    CONTINUE
 
      IF (K.EQ.N) GOTO 60
      IF (N.EQ.1) THEN
      Z(1,2) = -Z(2,1)
      Z(2,2) = Z(1,1)
      C = 0.0
      DO 30 I=1,K
30    C = C+Z(I,2)*Z(I,2)
      C=SQRT(C)
      IF (C.LT.1.0E-5) THEN
      Z(1,2) = 1.0
      Z(2,2) = 0.0
      ELSE
      DO 35 I=1,K
35    Z(I,2) = Z(I,2)/C
      ENDIF
      ENDIF
      IF (N.EQ.2) THEN
      Z(1,3) = Z(2,1)*Z(3,2)-Z(3,1)*Z(2,2)
      Z(2,3) = Z(3,1)*Z(1,2)-Z(1,1)*Z(3,2)
      Z(3,3) = Z(1,1)*Z(2,2)-Z(2,1)*Z(1,2)
      ENDIF
60    CONTINUE
      DO 90 I=1,M
      DO 80 J=1,N
      C=0.0
      DO 70 L=1,K
70    C=C+R(L,I)*Z(L,J)
      Y(J,I)=C
80    CONTINUE
90    CONTINUE
      DO 100 I=1,K
      DO 100 J=1,I-1
      ZIJ = Z(I,J)
      Z(I,J) = Z(J,I)
      Z(J,I) = ZIJ
100   CONTINUE
 
c      WRITE(*,*) 'Z ='
c      DO 200 I=1,K
c      WRITE(*,*) (Z(J,I),J=1,K)
c200   CONTINUE
c      WRITE(*,*) 'Y ='
c      DO 300 I=1,M
c      WRITE(*,*) (Y(J,I),J=1,N)
c300   CONTINUE
      END
 
 
      SUBROUTINE TKT (NT,NA,T,A,TAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION  T(NA,NT) , A(NA,NA) , TA(50) , TAT(NT,NT)
C.......COMPUTE TAT = T(t)*A*T
      DO 50 I = 1 , NT
      DO 20 K = 1 , NA
      CC = 0.0
      DO 10 L = 1 , NA
   10 CC = CC+T(L , I)*A(L , K)
   20 TA(K) = CC
      DO 40 J = 1 , NT
      DD = 0.0
      DO 30 K = 1 , NA
   30  DD = DD + TA(K)*T(K , J)
      TAT(I ,J) = DD
   40 CONTINUE
   50 CONTINUE
      RETURN
      END
 
      SUBROUTINE TMT (NT,NA,T,F,TF)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION  T(NA,NT) , F(NA), TF(NT)
C.......CUMPUTE TF = T(t)*F*T
      DO 20 I = 1 , NT
      CC = 0.0
      DO 10 L = 1 , NA
   10 CC = CC+T(L,I)*F(L)*T(L,I)
   20 TF(I) = CC
      RETURN
      END

      SUBROUTINE NTMT (NT,NA,T,F,TF)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION  T(NA,NT) , F(NA), TF(NT)
C.......CUMPUTE TF = F
      DO 10 L = 1 , NA
   10 TF(L) = F(L)
      RETURN
      END
       
      SUBROUTINE TL (NT,NA,T,F,TF)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION  T(NA,NT) , F(NA), TF(NT)
C.......CUMPUTE TF = T(t)*F
      DO 20 I = 1 , NT
      CC = 0.0
      DO 10 L = 1 , NA
   10 CC = CC+T(L , I)*F(L)
   20 TF(I) = CC
      RETURN
      END
      SUBROUTINE CCTL(NGVAR,NLVAR,EGS,EGM,EGL,ELS,ELM,ELL)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION EGS(NGVAR,NGVAR),EGM(NGVAR),EGL(NGVAR),
     *         ELS(NLVAR,NLVAR),ELM(NLVAR),ELL(NLVAR)
      LOGICAL FAIL
      CALL INVER(NGVAR,NLVAR,EGS,ELS,FAIL)
      CALL SOLVE(NGVAR,NLVAR,EGS,EGL,ELL)
      DO 10 I=1,NLVAR
10    ELM(I) = EGM(I+NGVAR-NLVAR)
      RETURN
      END
 
      SUBROUTINE CCTD(NGVAR,NLVAR,EGS,EGM,EGL,ELS,ELM,ELL)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION EGS(NGVAR,NGVAR),EGM(NGVAR,NGVAR),EGL(NGVAR),
     *         ELS(NLVAR,NLVAR),ELM(NLVAR,NLVAR),ELL(NLVAR)
      LOGICAL FAIL
      CALL INVER(NGVAR,NLVAR,EGS,ELS,FAIL)
      CALL SOLVE(NGVAR,NLVAR,EGS,EGL,ELL)
      CALL INVER(NGVAR,NLVAR,EGM,ELM,FAIL)
      RETURN
      END
 
      SUBROUTINE INVER(N,M,A,D,FAIL)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N,N),D(M,M)
      LOGICAL FAIL
      FAIL = .FALSE.
C     WRITE (*,*) 'N,M =',N,M
CC    K = N-M
      DO 1000 I = 1,N
      IF (I.LE.N-M) THEN
      I1 = I
      ELSE
      I1 = N-M
      ENDIF
       DO 900 J = 1,I1
        R0 = A(I,J)
100     IF (I .EQ. 1) GO TO 300
        J1 = J-1
        DO 200 L = 1,J1
200      R0 = R0 - A(I,L) * A(J,L)
300     IF (J .NE. I) GO TO 500
      IF (R0 .LT. 1E-19) THEN
        FAIL = .TRUE.
        WRITE (*,2000) I,R0
2000    FORMAT (1X,4X,17HA is singular  I=,I4,9H  Lii**2=,E15.7)
        RETURN
      ELSE
        A(I,I) = SQRT(R0)
        GO TO 1000
      ENDIF
500     A(I,J) = R0 / A(J,J)
900   CONTINUE
1000  CONTINUE
 
      DO 1500 I=1,M
      DO 1400 J=1,M
      R0 = A(N-M+I,N-M+J)
      DO 1100 L=1,N-M
      R0 = R0-A(N-M+I,L)*A(N-M+J,L)
1100  CONTINUE
      D(I,J) = R0
1400   CONTINUE
1500  CONTINUE
 
      RETURN
      END
 
      SUBROUTINE SOLVE(N,M,A,F,G)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N,N),F(N),G(N)
C     WRITE(*,*) 'N,M =',N,M
      DO 200 I=1,N-M
       S0 = F(I)
       IF (I .EQ. 1) GO TO 200
       I1 = I - 1
       DO 100 J = 1,I1
100     S0 = S0 - A(I,J) * F(J)
200    F(I) = S0 / A(I,I)
      DO 400 I = 1,M
      S0 = F(N-M+I)
       DO 300 J = 1,N-M
300   S0 = S0 - A(N-M+I,J)*F(J)
      G(I) = S0
400    CONTINUE
      RETURN
      END
      SUBROUTINE QMESH(NX,NY,X,Y)
C .... REAL*4 FUNCTION QQINT(X,Y,N1,N2,N3,N4,U)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),Y(*),XB(401),YB(401)
      N1 = NX+1
      N2 = NY+1
      N3 = N1
      N4 = N2
      NTOL = (NX+NY)*2
      DO 100 N=1,NTOL
      XB(N) = X(N)
      YB(N) = Y(N)
100   CONTINUE
      XB(NTOL+1)  = XB(1)
      YB(NTOL+1)  = YB(1)
      DO 300 J=1,NY+1
      DO 200 I=1,NX+1
      XM = 1.0*(I-1)/NX
      YM = 1.0*(J-1)/NY
      M = (J-1)*(NX+1)+I
      X(M) = XB(I)*(1.0-YM)+XB(NX*2+NY-I+2)*YM
     &     + XB(NTOL-J+2)*(1.0-XM)+XB(NX+J)*XM
      Y(M) = YB(I)*(1.0-YM)+YB(NX*2+NY-I+2)*YM
     &     + YB(NTOL-J+2)*(1.0-XM)+YB(NX+J)*XM
      X(M) = X(M)-XB(1)*(1.-XM)*(1.-YM)-XB(NX+1)*XM*(1.-YM)
     &      -XB(NX+NY+1)*XM*YM-XB(NTOL-NY+1)*(1.-XM)*YM
      Y(M) = Y(M)-YB(1)*(1.-XM)*(1.-YM)-YB(NX+1)*XM*(1.-YM)
     &      -YB(NX+NY+1)*XM*YM-YB(NTOL-NY+1)*(1.-XM)*YM
C     WRITE(*,*) 'I,J,M,XM,YM =',I,J,M,XM,YM
C     X(M) = QQINT(XM,YM,N1,N2,N3,N4,XB)
C     Y(M) = QQINT(XM,YM,N1,N2,N3,N4,YB)
C     WRITE(*,*) 'X,Y =',X(M),Y(M)
200   CONTINUE
300   CONTINUE
      RETURN
      END
 
      SUBROUTINE TMESH(NX,X,Y)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),Y(*),XB(301),YB(301)
      NTOL = NX*3
      DO 100 N=1,NTOL
      XB(N) = X(N)
      YB(N) = Y(N)
100   CONTINUE
      M = (NX+2)*(NX+1)/2
      X(M) = XB(NX*2+1)
      Y(M) = YB(NX*2+1)
C     WRITE(*,*) 'M,X,Y ..........',M,X(M),Y(M)
      M = 0
      DO 300 J=1,NX
      DO 200 I=1,NX+2-J
      XM = 1.0*(I-1)/NX
      YM = 1.0*(J-1)/NX
      ZM = (1.0-XM-YM)*0.5
      MM = NX*3+2-I+J
      M2 = MM/2
      M1 = MM - M2
C     WRITE(*,*) 'J,I,XM,YM,ZM =',J,I,XM,YM,ZM
C     WRITE(*,*) 'M1,M2,MM =',M1,M2,MM
      M = M+1
      IF (I.EQ.1 .AND. J.EQ.1) GOTO 200
      IF (I.EQ.NX+1 .AND. J.EQ.1) GOTO 200
      X(M) = XB(I)*XM*ZM+XB(3*NX-J+2)*YM*ZM+(XB(M1)+XB(M2))*XM*YM/2
      X(M) = X(M)/(XM*YM+YM*ZM+ZM*XM)
      Y(M) = YB(I)*XM*ZM+YB(3*NX-J+2)*YM*ZM+(YB(M1)+YB(M2))*XM*YM/2
      Y(M) = Y(M)/(XM*YM+YM*ZM+ZM*XM)
C     WRITE(*,*) 'M,X,Y ..........',M,X(M),Y(M)
200   CONTINUE
300   CONTINUE
      RETURN
      END
 
c      REAL*4 FUNCTION QQINT(X,Y,N1,N2,N3,N4,U)
c      IMPLICIT REAL*4 (A-H,O-Z)
c      DIMENSION U(1),V(100)
c      WRITE(*,*) 'X,Y =',X,Y
c      M1=1
c      M2=N1
c      M3=M2+N2-1
c      M4=M3+N3-1
c      M =M4+N4-2
c      B = U(M1)*(1.-X)*(1.-Y)+U(M2)*X*(1.-Y)
c     *  +U(M3)*X*Y + U(M4)*(1.-X)*Y
cc      WRITE(*,*) 'M1,M2,M3,M,B=',M1,M2,M3,M,B
cc      WRITE(*,*) 'U =',(U(I),I=1,M)
c      C1 = FLINT(X,N1,U(M1))*(1.-Y)
c      C2 = FLINT(Y,N2,U(M2))*X
c      DO 10 I=1,N3
c      V(I) = U(M4-I+1)
c10    CONTINUE
c      C3 = FLINT(X,N3,V)*Y
c      V(1) = U(1)
c      DO 20 I=2,N4
c      V(I) = U(M-I+2)
c20    CONTINUE
c      C4 = FLINT(Y,N4,V)*(1.-X)
c      QQINT = C1 + C2 + C3 + C4 - B
cc      WRITE (*,*) 'C1,C2,C3,C4,B =',C1,C2,C3,C4,B
c      RETURN
c      END
 
      SUBROUTINE PLINE(N0,N1,M,X,Y)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),Y(*),XL(200),YL(200),s(200)
      if (M.lt.2) then
      write(*,*) 'M ........ ',M
      write(*,*) 'NUMBER OF PLINE NODES CANNOT < 2'
      endif
      NL = 0
      DO 100 N=N0,N1
      NL = NL+1
      XL(NL) = X(N)
      YL(NL) = Y(N)
100   CONTINUE
c      write(*,*) 'NL,N0,N1 = ',NL,N0,N1,'  x & y ='
c      write(*,*) (x(i),i=n0,n1)
c      write(*,*) (y(i),i=n0,n1)
c      write(*,*) 'NL =',NL,'  M = ',M
      call gets(NL,XL,YL,s)
c      write(*,*) 'NL =',NL
      DO 200 I=1,M
      C = s(NL)*(I-1)/(M-1)
      X(N0+I-1) = FLINT(C,NL,XL,s)
      Y(N0+I-1) = FLINT(C,NL,YL,s)
200   CONTINUE
      N1 = N0+M-1
c      write(*,*) 'NL,N0,N1 = ',NL,N0,N1,'  x & y ='
c      write(*,*) (x(i),i=n0,n1)
c      write(*,*) (y(i),i=n0,n1)
      RETURN
      END
      
      REAL*8 FUNCTION FLINT(Y,N,U,X)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION U(N),X(200)
      F = 0.0
      DO 50 I=1,N
      C = U(I)
      DO 30 J=1,N
      IF (J.NE.I) C = C*(Y-X(J))/(X(I)-X(J))
30    CONTINUE
      F = F + C
50    CONTINUE
      FLINT = F
C     WRITE(*,*) 'N,Y,F = ',N,Y,F
C     WRITE(*,*) (X(I),I=1,N)
C     WRITE(*,*) (U(I),I=1,N)
      RETURN
      END
 
      subroutine gets1(n,x,y,s)
      implicit REAL*8 (a-h,o-z)
      dimension s(*),x(*),y(*)
c     write(*,*) 'n = ',n
      s0=sqrt((x(n)-x(1))**2+(y(n)-y(1))**2)*1.0d-6
      s(1)=0.0
      do 100 k=2,n
      s(k)=sqrt((x(k)-x(k-1))**2+(y(k)-y(k-1))**2)
100   continue
      m = 1
      do 200 k=2,n
      if (s(k).gt.s0) then
      m = m+1
      s(m) = s(k)
      x(m) = x(k)
      y(m) = y(k)
      endif
200   continue
      do 300 k=2,m
      s(k) = s(k)+s(k-1)
300   continue
c     WRITE(*,*) 'm =',m,'  s = '
c     write(*,*) (s(i),i=1,m)
      n = m
      return
      end
 
      subroutine getxy(n,x,y,s)
      implicit REAL*8 (a-h,o-z)
      dimension s(*),x(*),y(*)
c     write(*,*) '-------- getxy n = ',n
      s0=(x(n)-x(1))*1.0d-6
      s(1)=0.0
      do 100 k=2,n
      s(k)=x(k)-x(k-1)
100   continue
      m = 1
      do 200 k=2,n
      if (s(k)**2.gt.s0**2) then
      m = m+1
      s(m) = s(k)
      x(m) = x(k)
      y(m) = y(k)
      endif
200   continue
      do 300 k=2,m
      s(k) = s(k)+s(k-1)
300   continue
c     WRITE(*,*) 'm =',m,'  s = '
c     write(*,*) (s(i),i=1,m)
      n = m
      return
      end
      
      subroutine xplin(n0,n1,m,x,y)
      implicit REAL*8 (a-h,o-z)
      dimension tl(200),xl(200),yl(200),tt(1000),x(*),y(*)
      do 500 k=n0,n1
      nl = k-n0+1
      xl(nl) = x(k)
      yl(nl) = y(k)
500   continue
c     WRITE(*,*) 'm,nl =',m,nl,'  xl,yl ='
c     WRITE(*,*) (xl(i),yl(i),i=1,nl)
      WRITE(*,*) 'm,nl =',m,nl,'  xl,yl ='
      call getxy(nl,xl,yl,tl)
      WRITE(*,*) (xl(i),i=1,nl)
      WRITE(*,*) (yl(i),i=1,nl)
       write(*,*) 'NL =',NL
      if (nl.le.2) then
      dx = (xl(2)-xl(1))/(m-1)
      dy = (yl(2)-yl(1))/(m-1)
      do 600 k=1,m-1
      x(k+n0) = x(n0)+dx*k
      y(k+n0) = y(n0)+dy*k
600   continue
      else
      dt=tl(nl)/(m-1)
      tt(1)=0.0
      do 800 i=2,m
      tt(i)=tt(i-1)+dt
800   continue
      write(*,*) 'tt =',(tt(i),i=1,m)
      call splin1(nl,m,tl,xl,tt,x(n0))
      call splin1(nl,m,tl,yl,tt,y(n0))
      endif
      N1 = N0+M-1
c      write(*,*) 'NL,N0,N1 = ',NL,N0,N1,'  x & y ='
c      write(*,*) (x(i),i=n0,n1)
c      write(*,*) (y(i),i=n0,n1)
      return
      end
 
      subroutine splin(n0,n1,m,x,y)
      implicit REAL*8 (a-h,o-z)
      dimension tl(200),xl(200),yl(200),tt(1000),x(*),y(*)
      do 500 k=n0,n1
      nl = k-n0+1
      xl(nl) = x(k)
      yl(nl) = y(k)
500   continue
c     WRITE(*,*) 'm,nl =',m,nl,'  xl,yl ='
c     WRITE(*,*) (xl(i),yl(i),i=1,nl)
      WRITE(*,*) 'm,nl =',m,nl,'  xl,yl ='
      call gets(nl,xl,yl,tl)
      WRITE(*,*) (xl(i),i=1,nl)
      WRITE(*,*) (yl(i),i=1,nl)
       write(*,*) 'NL =',NL
      if (nl.le.2) then
      dx = (xl(2)-xl(1))/(m-1)
      dy = (yl(2)-yl(1))/(m-1)
      do 600 k=1,m-1
      x(k+n0) = x(n0)+dx*k
      y(k+n0) = y(n0)+dy*k
600   continue
      else
      dt=tl(nl)/(m-1)
      tt(1)=0.0
      do 800 i=2,m
      tt(i)=tt(i-1)+dt
800   continue
      write(*,*) 'tt =',(tt(i),i=1,m)
      call splin1(nl,m,tl,xl,tt,x(n0))
      call splin1(nl,m,tl,yl,tt,y(n0))
      endif
      N1 = N0+M-1
c      write(*,*) 'NL,N0,N1 = ',NL,N0,N1,'  x & y ='
c      write(*,*) (x(i),i=n0,n1)
c      write(*,*) (y(i),i=n0,n1)
      return
      end
 
      subroutine splin1(nl,m,xl,yl,xx,yy)
      implicit REAL*8 (a-h,o-z)
      dimension xl(*),yl(*),dy(200),
     *            h(200),xx(*),yy(*)
c        double precision x,y,dy,ddy,h,h0,h1,a,b
c     WRITE(*,*) 'nl,m =',nl,m,'  xl,yl ='
c     WRITE(*,*) (xl(i),i=1,nl)
c     WRITE(*,*) (yl(i),i=1,nl)
c     WRITE(*,*) 'xx =',(xx(i),i=1,m)
      dy(1)=-0.5
      h0=xl(2)-xl(1)
      h(1)=3.0*(yl(2)-yl(1))/(2.0*h0)
      do 10 j=2,nl-1
      h1=xl(j+1)-xl(j)
      a=h0/(h0+h1)
      b=(1.0-a)*(yl(j)-yl(j-1))/h0
      b=3.0*(b+a*(yl(j+1)-yl(j))/h1)
      dy(j)=-a/(2.0+(1.0-a)*dy(j-1))
      h(j)=(b-(1.0-a)*h(j-1))
      h(j)=h(j)/(2.0+(1.0-a)*dy(j-1))
      h0=h1
10    continue
      dy(nl)=(3.0*(yl(nl)-yl(nl-1))/h1-h(nl-1))/(2.0+dy(nl-1))
        do 20 j=nl-1,1,-1
      dy(j)=dy(j)*dy(j+1)+h(j)
20    continue
      do 30 j=1,nl-1
      h(j)=xl(j+1)-xl(j)
30    continue
      do 70 j=1,m
       if (xx(j).ge.xl(nl)) then
         i=nl-1
       else
         i=1
60       if (xx(j).gt.xl(i+1)) then
            i=i+1
          goto 60
         endif
       endif
      h1=(xl(i+1)-xx(j))/h(i)
      yy(j)=(3.0*h1*h1-2.0*h1*h1*h1)*yl(i)
      yy(j)=yy(j)+h(i)*(h1*h1-h1*h1*h1)*dy(i)
      h1=(xx(j)-xl(i))/h(i)
      yy(j)=yy(j)+(3.0*h1*h1-2.0*h1*h1*h1)*yl(i+1)
      yy(j)=yy(j)-h(i)*(h1*h1-h1*h1*h1)*dy(i+1)
70    continue
c     WRITE(*,*) 'yy =',(yy(i),i=1,m)
      return
      end
 
      subroutine exps1(n,a,b,s)
      implicit real*8 (a-h,o-z)
      dimension s(100)
      if (n.lt.2) then
      write(*,*) 'Error number of nodes < 2 '
      stop
      return
      endif
      c = 1.0
      d = 0.0
      e = b*(n-1)
      if (e.lt.1.0) then
      c = 1.0 - e
      d = b
      endif
      do 100 i=1,n
      f = 1.0*(i-1)/(n-1)
c     s(i) = c*(i-1)**a/(n-1)**a+d*(i-1)
      s(i) = c*f**a+d*(i-1)
c        write(*,*) s(i)
100     continue
      return
      end
 
      subroutine exps2(n,a,b,s)
      implicit real*8 (a-h,o-z)
      dimension s(100)
      call exps1(n,a,b,s)
      m = n/2
      do 100 i=1,m
      c = 1.-s(i)
      s(i) = 1.-s(n-i+1)
      s(n-i+1) = c
100   continue
      kk = n - n/2*2
      if (kk.gt.0) s(m+1) = 1.0-s(m+1)
      return
      end
 
      subroutine exps3(n,a,bb,s)
      implicit real*8 (a-h,o-z)
      dimension s(100)
      m = n/2+1
      b = bb*2.0
c     write(*,*) 'n,m = ',n,m
      call exps1(m,a,b,s)
      kk = n - n/2*2
      if (kk.eq.0) then
      k = m-1
      else
      k = m
      endif
c     call exps2(m,a,b,s(k))
      do 100 i=1,k
      if (kk.eq.0) s(i) = s(i)*n/(n-1)
      s(i) = s(i)/2.0
      s(n-i+1) = 1-s(i)
100   continue
      return
      end
 
      subroutine spline(n0,n1,m,x,y,s)
      implicit REAL*8 (a-h,o-z)
      dimension tl(200),xl(200),yl(200),tt(1000),x(*),y(*),s(*)
      WRITE(*,*) 'n0,n1,m =',n0,n1,m,'  s(1)_s(m) ='
      WRITE(*,'(1x,5f13.4)') (s(i),i=1,m)
      do 500 k=n0,n1
      nl = k-n0+1
      xl(nl) = x(k)
      yl(nl) = y(k)
500   continue
c     WRITE(*,*) 'm,nl =',m,nl,'  xl,yl ='
c     WRITE(*,*) (xl(i),yl(i),i=1,nl)
      WRITE(*,*) 'm,nl =',m,nl,'  xl,yl ='
      call gets(nl,xl,yl,tl)
      WRITE(*,*) (xl(i),i=1,nl)
      WRITE(*,*) (yl(i),i=1,nl)
       write(*,*) 'NL =',NL,'  n0 = ',n0
      if (nl.le.2) then
c     dx = (xl(2)-xl(1))/(m-1)
c     dy = (yl(2)-yl(1))/(m-1)
      dx = xl(2)-xl(1)
      dy = yl(2)-yl(1)
      do 600 k=1,m-1
      x(k+n0) = x(n0)+dx*s(k+1)
      y(k+n0) = y(n0)+dy*s(k+1)
600   continue
      else
c     dt=tl(nl)/(m-1)
      dt=tl(nl)
      tt(1)=0.0
      do 800 k=2,m
c     tt(k)=tt(k-1)+dt
      tt(k)=dt*s(k)
800   continue
      write(*,*) 'tt =',(tt(i),i=1,m)
      call splin1(nl,m,tl,xl,tt,x(n0))
      call splin1(nl,m,tl,yl,tt,y(n0))
      endif
      N1 = N0+M-1
c      write(*,*) 'NL,N0,N1 = ',NL,N0,N1,'  x & y ='
c      write(*,*) (x(i),i=n0,n1)
c      write(*,*) (y(i),i=n0,n1)
      return
      end
 
      SUBROUTINE SHAP(NREFC,NCOOR,NVAR,SHPR,SHPC,CR,
     *                   DORD,TOLC,TOLR)
      IMPLICIT REAL*8 (A-H,O-Z)
C     IMPLICIT INTEGER*2 (I-N)
      INTEGER TOLC,TOLR,DORD
      DIMENSION SHPR(NVAR,TOLR),SHPC(NVAR,TOLC),CR(NREFC,NCOOR)
C     WRITE(*,*) 'NREFC,NCOOR,NVAR,DORD,TOLR,TOLC='
C     WRITE(*,9) NREFC,NCOOR,NVAR,DORD,TOLR,TOLC
C     WRITE(*,*) 'SHPR ='
C     DO 21 I=1,NVAR
C21   WRITE(*,8) (SHPR(I,J),J=1,TOLR)
      DO 100 I=1,NVAR
      LR=1
      LC=1
      SHPC(I,LC)=SHPR(I,LR)
      DO 20 J=1,NCOOR
      C=0.0
      DO 10 K=1,NREFC
10    C = C + SHPR(I,LR+K)*CR(K,J)
20    SHPC(I,LC+J) = C
      LC = LC + NCOOR
      LR = LR + NREFC
      IF (DORD.EQ.1) GOTO 100
      DO 60 J1=1,NCOOR
      DO 50 J2=J1,NCOOR
      J = (2*NCOOR-J1)*(J1-1)/2+J2
      C = 0.0
      DO 40 K1=1,NREFC
      DO 30 K2=K1,NREFC
      K = (2*NREFC-K1)*(K1-1)/2+K2
      C = C + SHPR(I,LR+K)*CR(K1,J1)*CR(K2,J2)
      IF (K1.LT.K2) C = C + SHPR(I,LR+K)*CR(K2,J1)*CR(K1,J2)
30    CONTINUE
40    CONTINUE
50    SHPC(I,LC+J) = C
60    CONTINUE
100   CONTINUE
C     WRITE(*,*) 'SHPC ='
C     DO 22 I=1,NVAR
C22   WRITE(*,8) (SHPC(I,J),J=1,TOLC)
8     FORMAT(1X,8F9.2)
9     FORMAT(1X,15I4)
      END
 
      SUBROUTINE SHAPn(NREFC,NCOOR,NVAR,SHPR,SHPC,CR,
     *                   DORD,TOLC,TOLR)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER TOLC,TOLR,DORD
      DIMENSION SHPR(NVAR,TOLR),SHPC(NVAR,TOLC),CR(NREFC,NCOOR)
C     WRITE(*,*) 'NREFC,NCOOR,NVAR,DORD,TOLR,TOLC='
C     WRITE(*,9) NREFC,NCOOR,NVAR,DORD,TOLR,TOLC
C     WRITE(*,*) 'SHPR ='
C     DO 21 I=1,NVAR
C21   WRITE(*,8) (SHPR(I,J),J=1,TOLR)
      DO 100 I=1,NVAR
      LR=1
      LC=1
      SHPC(I,LC)=SHPR(I,LR)
      DO 20 J=1,NCOOR
      C=0.0
      DO 10 K=1,NREFC
10    C = C + SHPR(I,LR+K)*CR(K,J)
20    SHPC(I,LC+J) = C
      LC = LC + NCOOR
      LR = LR + NREFC
      IF (DORD.EQ.1) GOTO 100
      DO 60 J1=1,NCOOR
      DO 50 J2=J1,NCOOR
      J = (2*NCOOR-J1)*(J1-1)/2+J2
      C = 0.0
      DO 40 K1=1,NREFC
      DO 30 K2=K1,NREFC
      K = (2*NREFC-K1)*(K1-1)/2+K2
      C = C + SHPR(I,LR+K)*CR(K1,J1)*CR(K2,J2)
      IF (K1.LT.K2) C = C + SHPR(I,LR+K)*CR(K2,J1)*CR(K1,J2)
30    CONTINUE
40    CONTINUE
50    SHPC(I,LC+J) = C
60    CONTINUE
100   CONTINUE
C     WRITE(*,*) 'SHPC ='
C     DO 22 I=1,NVAR
C22   WRITE(*,8) (SHPC(I,J),J=1,TOLC)
8     FORMAT(1X,6e13.4)
9     FORMAT(1X,15I5)
      return
      END
 
      SUBROUTINE SHAPc(NREFC,NCOOR,NVAR,SHPR,SHPC,CR,
     *                   DORD,TOLC,TOLR)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER TOLC,TOLR,DORD
      DIMENSION SHPR(NVAR,TOLR),SHPC(NVAR,TOLC),CR(NREFC,NCOOR)
C     WRITE(*,*) 'NREFC,NCOOR,NVAR,DORD,TOLR,TOLC='
C     WRITE(*,9) NREFC,NCOOR,NVAR,DORD,TOLR,TOLC
C     WRITE(*,*) 'SHPR ='
C     DO 21 I=1,NVAR
C21   WRITE(*,8) (SHPR(I,J),J=1,TOLR)
      DO 100 I=1,NVAR
      LR=0
      LC=0
      DO 20 J=1,NCOOR
      C=0.0
      DO 10 K=1,NREFC
10    C = C + SHPR(I,LR+K)*CR(K,J)
20    SHPC(I,LC+J) = C
      LC = LC + NCOOR
      LR = LR + NREFC
      IF (DORD.EQ.1) GOTO 100
      DO 60 J1=1,NCOOR
      DO 50 J2=J1,NCOOR
      J = (2*NCOOR-J1)*(J1-1)/2+J2
      C = 0.0
      DO 40 K1=1,NREFC
      DO 30 K2=K1,NREFC
      K = (2*NREFC-K1)*(K1-1)/2+K2
      C = C + SHPR(I,LR+K)*CR(K1,J1)*CR(K2,J2)
      IF (K1.LT.K2) C = C + SHPR(I,LR+K)*CR(K2,J1)*CR(K1,J2)
30    CONTINUE
40    CONTINUE
50    SHPC(I,LC+J) = C
60    CONTINUE
100   CONTINUE
C     WRITE(*,*) 'SHPC ='
C     DO 22 I=1,NVAR
C22   WRITE(*,8) (SHPC(I,J),J=1,TOLC)
8     FORMAT(1X,6E13.4)
9     FORMAT(1X,15I5)
      return
      END

c..................................................................
c     the following subroutines for generating triangular mesh
c       according to circle boundary nodal points
c     trimesh generates nodal points of mesh
c     trinod generates nodal No. of each triangle
c     tricb generates boundary edges according to
c       continuous nodal No. of circle boundary
c     tricbn generates boundary edges according to
c       discontinuous nodal No. of circle boundary
c..................................................................
c     trimesh(nums,numnod,numedge,mc,x,y)
c     trinod(numedge,nt,x,y,nod1,nod2,nod3)
c     tricb(numedge,ns,ne)
c     tricbn(numedge,numnod,node)
c..................................................................
c     nums ....... starting edge No. - 1
c     numnod ..... number of nodes
c     numedge .... number of edges
c     mc ......... number of iteration for modifying nodal points
c     nt ......... number of triangles
c     ns ......... first nodal point
c     ne ......... last nodal point
c..................................................................
c     x,y ............... coordinate arrays of nodal points
c     nod1,nod2,nod3 .... nodal No. arrays for each triangle
c     node .............. circle boundary nodal No. array
c..................................................................
 
      subroutine trimesh(nums,numnod,nume,mc,x,y)
c .... generate new nodal points for triagle elements
c .... store all nodal No. of edges in arrays nod1 & nod2
      implicit real*8 (a-h,o-z)
      common /edge/ nod1(8000),nod2(8000)
      dimension x(*),y(*),nv(100)
      COMMON /TMDATA/ index(8000)
      integer edge1,edge2
 
      numb = numnod
cc    write(*,*) 'numb =',numb,'  x,y ='
cc    write(*,'(1x,6f12.4)') (x(n),n=1,numnod)
cc    write(*,'(1x,6f12.4)') (y(n),n=1,numnod)
cc    write(*,*) 'nume =',nume,'  nod1,nod2 ='
cc    write(*,'(1x,15i5)') (nod1(n),n=1,nume)
cc    write(*,'(1x,15i5)') (nod2(n),n=1,nume)
 
      dmin = 1.e20
      dmax = 0.0
      do 50 n=1,nume
      index(n) = 0
      i = nod1(n)
      j = nod2(n)
      d = (x(i)-x(j))**2+(y(i)-y(j))**2
      if (d.lt.dmin) dmin = d
      if (d.gt.dmax) dmax = d
50    continue
      dmin = sqrt(dmin)
      dmax = sqrt(dmax)
      write(*,*) 'dmin,dmax =',dmin,dmax
      write(*,*) 'index ='
      write(*,'(1x,15i5)') (index(n),n=1,nume)
 
c ....  m is the middle point of n1 & n2
c ....  m_h is the normal vector of n1_n2
c ....  h is the relative high with the bottom n1_n2
cccc  nums = 0
      h = sqrt(3.0)
cccc  h = 1.d0
100   continue
      nums = nums+1
c     write(*,*) 'nums, index =',nums,index(nums)
      if (index(nums).gt.0) goto 2400
c     write(*,*) 'nums ,nume =',nums,nume
      n1 = nod1(nums)
      n2 = nod2(nums)
c     write(*,*) 'n1,n2 =',n1,n2
      xx1 = x(n1)
      yy1 = y(n1)
      xx2 = x(n2)
      yy2 = y(n2)
      xm = (xx1+xx2)/2
      ym = (yy1+yy2)/2
      a = xm-xx1
      b = ym-yy1
      dd = (a*a+b*b)*4
      d = sqrt(dd)
      hh = h
      if (d.lt.dmin) hh = h*dmin/d
      if (d.gt.dmax) hh = h*dmax/d
      xh = xm - b*hh
      yh = ym + a*hh
c     write(*,*) 'xh,yh =',xh,yh
      e = d*hh*1.0e-4
      ee = e*1.e-30
 
      nn = 0
      smin = dmax*dmax
      do 140 i=1,numnod
cccc  do 140 n=nums+1,nume
cccc  i = nod2(n)
      xx = x(i)
      yy = y(i)
      s = (xx-xh)**2+(yy-yh)**2
      if (s.lt.smin) then
      j = i
      smin = s
      endif
140   continue
      smin = sqrt(smin)
c     write(*,*) 'j,smin =',j,smin
cccc  if (smin.lt.dmin*0.5) then
      if (smin.lt.dmin*0.5 .or. smin.lt.d*0.5) then
      nn = j
c     write(*,*) 'nn,smin,d =',nn,smin,d
      xh = x(nn)
      yh = y(nn)
c     write(*,*) 'xh,yh =',xh,yh
      index(nums) = 2
      endif
 
cccc  do 1000 n=nums+1,nume
      do 1000 n=1,nume
      if (n.eq.nums) goto 1000
      k1 = nod1(n)
      k2 = nod2(n)
      if (k1.eq.nn .or. k2.eq.nn) goto 1000
      x1 = x(k1)
      x2 = x(k2)
      y1 = y(k1)
      y2 = y(k2)
      dm = dettri(x1,x2,xm,y1,y2,ym)
      dh = dettri(x1,x2,xh,y1,y2,yh)
      if (dm*dh.gt.-ee) goto 1000
      d1 = dettri(xm,xh,x1,ym,yh,y1)
      d2 = dettri(xm,xh,x2,ym,yh,y2)
      if (d1*d2.gt.-ee) goto 1000
c     write(*,*) 'dm,dh =',dm,dh
c     write(*,*) 'd1,d2 =',d1,d2
c     write(*,*) 'n,k1,k2 =',n,k1,k2
c     write(*,*) 'x1,x2,y1,y2 =',x1,x2,y1,y2
      if (k1.eq.n1 .or. k1.eq.n2) then
      nn = k2
      goto 500
      endif
      if (k2.eq.n1 .or. k2.eq.n2) then
      nn = k1
      goto 500
      endif
      if (d1*d1.lt.d2*d2) then
      nn = k1
      else
      nn = k2
      endif
500   continue
      xh = x(nn)
      yh = y(nn)
      index(nums) = 4
1000  continue
 
      do 400 n=1,numnod
      if (n.eq.n1 .or. n.eq.n2 .or. n.eq.nn) goto 400
      xn = x(n)
      yn = y(n)
      if (dettri(xn,xx1,xx2,yn,yy1,yy2).lt.-e) goto 400
      if (dettri(xn,xx2,xh,yn,yy2,yh).lt.-e) goto 400
      if (dettri(xn,xh,xx1,yn,yh,yy1).lt.-e) goto 400
      nn = n
      xh = xn
      yh = yn
      index(nums) = 3
c     write(*,*) 'n1,nn,n2,e =',n1,nn,n2,e
c     write(*,*) 'nn,xh,yh ============',nn,xh,yh
400   continue
 
c     write(*,*) 'nn,nums,index =',nn,nums,index(nums)
      if (nn.ne.0) goto 2100
      index(nums) = 1
      numnod = numnod+1
      x(numnod) = xh
      y(numnod) = yh
c     write(*,*) 'numnod,x,y =',numnod,x(numnod),y(numnod)
      nume = nume+1
      index(nume) = 0
      nod1(nume) = nod1(nums)
      nod2(nume) = numnod
      nod1(nume+1) = numnod
      nod2(nume+1) = nod2(nums)
      nume = nume+1
      index(nume) = 0
cccc  nums = nums+1
      goto 2400
 
2100  continue
      edge1 = 1
      edge2 = 1
c     write(*,*) 'n1,nn,n2 =',n1,nn,n2
cccc  do 2200 n=nums+1,nume
      do 2200 n=1,nume
c     write(*,*) 'n,nod1,nod2 =',n,nod1(n),nod2(n)
      if ((n1.eq.nod1(n) .and. nn.eq.nod2(n)) .or.
     &     (n1.eq.nod2(n) .and. nn.eq.nod1(n))) then
      index(n) = 5
c     write(*,*) 'n ..........',n
      edge1 = 0
      endif
      if ((n2.eq.nod1(n) .and. nn.eq.nod2(n)) .or.
     &     (n2.eq.nod2(n) .and. nn.eq.nod1(n))) then
      index(n) = 5
c     write(*,*) 'n ..........',n
      edge2 = 0
      endif
c     write(*,*) 'edge1,edge2 =',edge1,edge2
2200  continue
c     write(*,*) 'edge1,edge2 =',edge1,edge2
      if (edge1.eq.1) then
      nume = nume+1
      index(nume) = 0
      nod1(nume) = n1
      nod2(nume) = nn
      endif
      if (edge2.eq.1) then
      nume = nume+1
      index(nume) = 0
      nod1(nume) = nn
      nod2(nume) = n2
      endif
 
2400  continue
 
      if (nume.gt.nums) goto 100
 
      call delnod(numb,numnod,nume,nod1,nod2,x,y,nv)
 
c     mc = 20
      do 2800 k=1,mc
      do 2700 n=numb+1,numnod
      m = 0
      do 2600 i=1,nume
      if (n.eq.nod1(i)) then
      m = m+1
      nv(m) = nod2(i)
      goto 2600
      endif
      if (n.eq.nod2(i)) then
      m = m+1
      nv(m) = nod1(i)
      endif
2600  continue
      xc = 0.0
      yc = 0.0
      do 2650 i=1,m
      j = nv(i)
      xc = xc+x(j)
      yc = yc+y(j)
2650  continue
      x(n) = xc/m
      y(n) = yc/m
2700  continue
2800  continue
 
c     write(*,*) 'numnod =',numnod
c     write(*,'(1x,6f12.4)') (x(n),n=1,numnod)
c     write(*,'(1x,6f12.4)') (y(n),n=1,numnod)
c     write(*,*) 'nume =',nume,'  nod1,nod2 ='
c     write(*,'(1x,15i5)') (nod1(n),n=1,nume)
c     write(*,'(1x,15i5)') (nod2(n),n=1,nume)
cc    write(*,*) '   n  nod1 nod2 index '
cc    do 3000 n=1,nume
cc    write(*,'(1x,4i5)') n,nod1(n),nod2(n),index(n)
cc3000continue
 
      return
      end
 
      subroutine delnod(numb,numnod,nume,nod1,nod2,x,y,nv)
      implicit real*8 (a-h,o-z)
      dimension nod1(nume),nod2(nume),x(numnod),y(numnod),nv(*)
250     continue
      do 270  n=numb+1,numnod
      m = 0
      do 260  i=1,nume
      if (n.eq.nod1(i)) then
      m = m+1
      nv(m) = i
      goto 260
      endif
      if (n.eq.nod2(i)) then
      m = m+1
      nv(m) = i
      endif
      if (m.gt.3) goto 270
260     continue
      if (m.le.3) goto 300
270     continue
      return
300     continue
      do 310  i=n+1,numnod
      x(i) = x(i+1)
      y(i) = y(i+1)
310     continue
      nn = 0
      do 320  i=1,nume
      do 315 j=1,m
      if (i.eq.nv(j)) goto 320
315     continue
      nn = nn+1
      nod1(nn) = nod1(i)
      nod2(nn) = nod2(i)
      if (nod1(nn).gt.n) nod1(nn) = nod1(nn)-1
      if (nod2(nn).gt.n) nod2(nn) = nod2(nn)-1
320     continue
      numnod = numnod-1
      nume = nn
      goto 250
      end
 
      subroutine trinod(nume,nt,x,y,nod1,nod2,nod3)
c .... generate nodal No. of each triangle
c .... according to nodal No. of edges
      implicit real*8 (a-h,o-z)
      common /edge/ n1(8000),n2(8000)
      dimension nod1(*),nod2(*),nod3(*)
      dimension x(*),y(*)
 
cc    write(*,*) 'nume =',nume,'  n1,n2 ='
cc    write(*,'(1x,15i5)') (n1(n),n=1,nume)
cc    write(*,'(1x,15i5)') (n2(n),n=1,nume)
 
      nt = 0
      do 500 n=1,nume
      i = n1(n)
      j = n2(n)
      do 300 m=n+1,nume
      m1 = n1(m)
      m2 = n2(m)
      if (m1.eq.i .or. m2.eq.i) then
      if (m1.eq.i) then
      m12 = m2
      else
      m12 = m1
      endif
      do 100 k=n+1,nume
      k1 = n1(k)
      k2 = n2(k)
      if ((k1.eq.j .and. k2.eq.m12) .or.
     &     (k2.eq.j .and. k1.eq.m12)) then
      x1 = x(i)
      x2 = x(j)
      x3 = x(m12)
      y1 = y(i)
      y2 = y(j)
      y3 = y(m12)
      d = dettri(x1,x2,x3,y1,y2,y3)
      nt = nt+1
      if (d.gt.0.0d0) then
      nod1(nt) = i
      nod2(nt) = j
      else
      nod1(nt) = j
      nod2(nt) = i
      endif
      nod3(nt) = m12
      endif
100   continue
      endif
300   continue
500   continue
 
cc    write(*,*) 'nt =',nt,'    nod1,nod2,nod3 ='
cc    do 1000 n=1,nt
cc    write(*,'(1x,3i5)') nod1(n),nod2(n),nod3(n)
cc1000continue
 
      return
      end
 
      real*8 function dettri(x1,x2,x3,y1,y2,y3)
      implicit real*8 (a-h,o-z)
      dettri = x1*y2+x2*y3+x3*y1
     &       -y1*x2-y2*x3-y3*x1
c     write(*,*) 'x1,x2,x3,y1,y2,y3 ='
c     write(*,'(1x,3f12.4)') x1,x2,x3
c     write(*,'(1x,3f12.4)') y1,y2,y3
c     write(*,*) 'dettri =',dettri
      return
      end
 
      subroutine tricb(n,ns,ne)
c .... generate circle boundary edges
c .... according to circle boundary nodes
c .... ns first nodal No.
c .... ne last nodal No.
c .... n starting edge No. - 1
      implicit real*8 (a-h,o-z)
      common /edge/ nod1(8000),nod2(8000)
      do 100 i=ns,ne
      n = n+1
      nod1(n) = i
      nod2(n) = i+1
100   continue
      nod2(n) = ns
      return
      end
 
      subroutine tricbn(n,numnod,node)
c .... generate circle boundary edges
c .... according to circle boundary nodel No.
c .... numnod number of nodes
c .... node nodal No. array
c .... n starting edge No. - 1
      implicit real*8 (a-h,o-z)
      common /edge/ nod1(8000),nod2(8000)
      dimension node(numnod)
      do 100 i=1,numnod
      n = n+1
      nod1(n) = node(i)
      if (i.lt.numnod) nod2(n) = node(i+1)
100   continue
      nod2(n) = node(1)
      return
      end
 
 
      subroutine gqm4(nx,ny,x,y)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      integer n
      dimension px(4),py(4),p(4)
      m = 4
      do 50 i=1,m
      px(i) = x(i)
      py(i) = y(i)
   50 continue
      dx=2.0/nx
      dy=2.0/ny
      do  4 i=1,ny+1
      do  3 j=1,nx+1
      n=(i-1)*(nx+1)+j
c      if (n.le.0) goto  3
      x(n)=dx*(j-1)-1.0
      y(n)=dy*(i-1)-1.0
   3  continue
   4  continue
      do 200 i=1,n
      xi=x(i)
      yi=y(i)
      x(i) = 0.0
      y(i) = 0.0
      p(1)=(1.-xi)*(1.-yi)/4
      p(2)=(1.+xi)*(1.-yi)/4
      p(3)=(1.+xi)*(1.+yi)/4
      p(4)=(1.-xi)*(1.+yi)/4
      do 100 j=1,m
      x(i) = x(i)+p(j)*px(j)
      y(i) = y(i)+p(j)*py(j)
  100 continue
  200 continue
      return
      end
 
      subroutine gtnod4(nx,ny,nod1,nod2,nod3)
      implicit real*8 (a-h,o-z)
      dimension nod1(*),nod2(*),nod3(*)
      integer n,nod1,nod2,nod3
      n = -1
      do  6 j=1,ny
      do  5  i=1,nx
      n=n+2
c      if (n.le.0) goto  5
      nod1(n)= (j-1)*(nx+1)+i
      nod2(n)=(j-1)*(nx+1)+i+1
      nod3(n)=j*(nx+1)+i
   5  continue
   6  continue
      n = 0
      do  8 j=1,ny
      do  7  i=1,nx
      n=n+2
c      if (n.le.0) goto  7
      nod1(n)= j*(nx+1)+i+1
      nod2(n)=j*(nx+1)+i
      nod3(n)=(j-1)*(nx+1)+i+1
   7  continue
   8  continue
      return
      end
 
      subroutine gqm8(nx,ny,x,y)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      integer n
      dimension px(8),py(8),p(8)
      m = 8
      do 50 i=1,m
      px(i) = x(i)
      py(i) = y(i)
   50 continue
      dx=2.0/nx
      dy=2.0/ny
      do 10 i=1,ny+1
      do  9 j=1,nx+1
      n=(i-1)*(nx+1)+j
c      if (n.le.0) goto  9
      x(n)=dx*(j-1)-1.0
      y(n)=dy*(i-1)-1.0
   9  continue
  10  continue
      do 200 i=1,n
      xi=x(i)
      yi=y(i)
      x(i) = 0.0
      y(i) = 0.0
      p(1)=(1.-xi)*(1.-yi)*(-xi-yi-1.)/4
      p(3)=(1.+xi)*(1.-yi)*(+xi-yi-1.)/4
      p(5)=(1.+xi)*(1.+yi)*(+xi+yi-1.)/4
      p(7)=(1.-xi)*(1.+yi)*(-xi+yi-1.)/4
      p(2)=(1.-yi)*(1.-xi*xi)/2
      p(6)=(1.+yi)*(1.-xi*xi)/2
      p(4)=(1.+xi)*(1.-yi*yi)/2
      p(8)=(1.-xi)*(1.-yi*yi)/2
      do 100 j=1,m
      x(i) = x(i)+p(j)*px(j)
      y(i) = y(i)+p(j)*py(j)
  100 continue
  200 continue
      return
      end
 
      subroutine gtm3(nx,x,y)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      integer n
      dimension px(3),py(3),p(3)
      m = 3
      do 50 i=1,m
      px(i) = x(i)
      py(i) = y(i)
   50 continue
      dx=1.0/nx
      n = 0
      do 300 i=1,nx+1
      do 200 j=1,nx+2-i
      n = n+1
      s2 = dx*(j-1)
      s3 = dx*(i-1)
      s1 = 1.0-s2-s3
      p(1)=s1
      p(2)=s2
      p(3)=s3
      x(n) = 0.0
      y(n) = 0.0
      do 100 k=1,m
      x(n) = x(n)+p(k)*px(k)
      y(n) = y(n)+p(k)*py(k)
  100 continue
  200 continue
  300 continue
      return
      end
 
      subroutine gtnod3(nx,nod1,nod2,nod3)
      implicit real*8 (a-h,o-z)
      dimension nod1(*),nod2(*),nod3(*)
      integer n,nod1,nod2,nod3
      nod = 0
      n = 0
      do 300 i=1,nx
      ni = nx+2-i
      do 200 j=1,nx+1-i
      n=n+1
      nod1(n)= nod+j
      nod2(n)=nod+j+1
      nod3(n)=nod+j+ni
      if (j.eq.ni-1) goto 200
      n=n+1
      nod1(n)= nod+j+1+ni
      nod2(n)=nod+j+ni
      nod3(n)=nod+j+1
  200 continue
      nod = nod+ni
  300 continue
      return
      end
 
      subroutine gtm6(nx,x,y)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      integer n
      dimension px(6),py(6),p(6)
      m = 6
      do 50 i=1,m
      px(i) = x(i)
      py(i) = y(i)
   50 continue
      dx=1.0/nx
      n = 0
      do 300 i=1,nx+1
      do 200 j=1,nx+2-i
      n = n+1
      s2 = dx*(j-1)
      s3 = dx*(i-1)
      s1 = 1.0-s2-s3
      p(1)=s1*(2.*s1-1.)
      p(2)=4.*s1*s2
      p(3)=s2*(2.*s2-1.)
      p(4)=4.*s2*s3
      p(5)=s3*(2.*s3-1.)
      p(6)=4.*s3*s1
      x(n) = 0.0
      y(n) = 0.0
      do 100 k=1,m
      x(n) = x(n)+p(k)*px(k)
      y(n) = y(n)+p(k)*py(k)
  100 continue
  200 continue
  300 continue
      return
      end
 
      subroutine boun4(n,nl,ib)
      implicit real*8 (a-h,o-z)
      dimension ib(*)
      do 100 i=1,n+1
      goto (1,2,3,4), nl
1     ib(i) = i
      goto 100
2     ib(i) = (n+1)*i
      goto 100
3     ib(i) = (n+1)*n+i
      goto 100
4     ib(i) = (n+1)*(i-1)+1
100   continue
      return
      end
 
      subroutine boun3(n,nl,ib)
      implicit real*8 (a-h,o-z)
      dimension ib(*)
      do 100 i=1,n+1
      goto (1,2,3), nl
1     ib(i) = i
      goto 100
2     ib(i) = (2*n+3-i)*i/2
      goto 100
3     ib(i) = ((n+2)*(n+1)-(i+1)*i)/2+1
100   continue
      return
      end
      subroutine merge(numnod,m,x,y,np)
      implicit real*8 (a-h,o-z)
      dimension np(numnod),x(numnod),y(numnod)
c     write(*,*) 'x,y ='
c     write(*,*) x
c     write(*,*) y
      e = 1.e-5
      xmin = x(1)
      xmax = x(1)
      ymin = y(1)
      ymax = y(1)
      do 100 n=1,numnod
      if (x(n).lt.xmin) xmin = x(n)
      if (x(n).gt.xmax) xmax = x(n)
      if (y(n).lt.ymin) ymin = y(n)
      if (y(n).gt.ymax) ymax = y(n)
100   continue
      ex = (xmax-xmin)*e
      ex = ex*ex
      ey = (ymax-ymin)*e
      ey = ey*ey
      write(*,*) 'ex,ey =',ex,ey
      np(1) = 1
      do 300 n=2,numnod
      np(n) = n
      do 200 m=1,n-1
      dx = x(m)-x(n)
      dy = y(m)-y(n)
      if (dx*dx.lt.ex .and. dy*dy.lt.ey) then
c     write(*,*) 'n,m =',n,m
c     write(*,*) 'dx*dx,dy*dy =',dx*dx,dy*dy
c     write(*,*) 'xn,yn,xm,ym =',x(n),y(n),x(m),y(m)
      np(n) = m
      goto 300
      endif
200   continue
300   continue
cc    write(*,*) 'np =',np
      m = 0
      do 400 n=1,numnod
      k = np(n)
      if(k.eq.n) then
      m = m+1
      np(n) = m
      x(m) = x(k)
      y(m) = y(k)
      else
      np(n) = np(k)
      endif
400   continue
cc    write(*,*) 'np =',np
      return
      end
 
      subroutine mid(numnp,np,id)
      implicit real*8 (a-h,o-z)
      dimension np(numnp),id(numnp),nc(10)
      WRITE(*,*) 'NP & ID ='
      WRITE(*,'(1X,15I4)') NP
      WRITE(*,'(1X,15I4)') ID
      do 300 n=1,numnp
      k = np(n)
      i = 1
      nc(1) = n
      do 100 m=n+1,numnp
      if (np(m).ne.k) goto 100
      i = i+1
      nc(i) = m
100   continue
c     write(*,*) 'i =',i,'  nc =',(nc(j),j=1,i)
      if (i.le.1) goto 300
      do 200 j=1,i
      jc = nc(j)
      if (id(jc).lt.1) then
      do 150 l=1,i
      lc = nc(l)
      id(lc) = id(jc)
150   continue
      goto 300
      endif
200   continue
300   continue
      do 400 n=1,numnp
      k = np(n)
      id(k) = id(n)
400   continue
      return
      end
 
      subroutine getname(name,IT)
      implicit real*8 (a-h,o-z)
      character name*12,ch3*3
c     IF (IT.LT.10) WRITE(UNIT=CH3,FMT='(I1)') IT
c     IF (IT.GE.10) WRITE(UNIT=CH3,FMT='(I2)') IT
c     IF (IT.GE.100) WRITE(UNIT=CH3,FMT='(I3)') IT
      call getext(it,ch3)
c     write(*,*) 'name =',name,'++++ CH3 =',CH3
      do 10 i=1,12
      if (name(i:i).eq.' ') then
      j=i
      goto 20
      endif
10    continue
20    continue
      if (j.gt.9) then
      write(*,*) 'Error, plot filename too long .......',name
      write(*,*) ' the length of name must be less or equal 8 character'
      stop 111
      endif
c     read(*,'(a3)') ch3
      name(j:j)='.'
      name(j+1:j+4)=ch3
c     write(*,*) 'plot filename = ',name
      return
      end
 
      subroutine getext(ii,ch3)
      implicit real*8 (a-h,o-z)
      character ch3*3
      it = ii
      ch3 = '   '
      k = 0
      if (ii.ge.100) then
      n = it/100
      k = k+1
      call getchar(n,k,ch3)
      it = it - n*100
      endif
      if (ii.ge.10) then
      n = it/10
      k = k+1
      call getchar(n,k,ch3)
      it = it - n*10
      endif
      n = it
      k = k+1
      call getchar(n,k,ch3)
      return
      end
 
      subroutine getchar(n,k,ch3)
      implicit real*8 (a-h,o-z)
      character ch3*3
      if (n.eq.0) ch3(k:k) = '0'
      if (n.eq.1) ch3(k:k) = '1'
      if (n.eq.2) ch3(k:k) = '2'
      if (n.eq.3) ch3(k:k) = '3'
      if (n.eq.4) ch3(k:k) = '4'
      if (n.eq.5) ch3(k:k) = '5'
      if (n.eq.6) ch3(k:k) = '6'
      if (n.eq.7) ch3(k:k) = '7'
      if (n.eq.8) ch3(k:k) = '8'
      if (n.eq.9) ch3(k:k) = '9'
      return
      end

      subroutine gets(n,x,y,s)
      implicit REAL*8 (a-h,o-z)
      dimension s(*),x(*),y(*)
c     write(*,*) 'n = ',n
      s0=sqrt((x(n)-x(1))**2+(y(n)-y(1))**2)*1.0d-6
      s(1)=0.0
      do 100 k=2,n
      s(k)=sqrt((x(k)-x(k-1))**2+(y(k)-y(k-1))**2)
100   continue
      m = 1
      do 200 k=2,n
      if (s(k).gt.s0) then
      m = m+1
      s(m) = s(k)
      x(m) = x(k)
      y(m) = y(k)
      endif
200   continue
      do 300 k=2,m
      s(k) = s(k)+s(k-1)
300   continue
c     WRITE(*,*) 'm =',m,'  s = '
c     write(*,*) (s(i),i=1,m)
      n = m
      return
      end
 
 
      subroutine dshap(f,x,shap,nrefc,nvar,ndord)
      implicit real*8 (a-h,o-z)
      dimension x(3),shap(nvar,*),dfdx(10)
      external f
c      write(*,*) 'dshap sub nvar ====',nvar
      do 200 n=1,nvar
c      write(*,*) 'dshap sub n ====',n
      call sdfdx(f,x,nrefc,dfdx,ndord,n,m)
      do 100 i=1,m
      shap(n,i) = dfdx(i)
100   continue
200   continue
      return
      end
 
      subroutine dcoef(f,x,cc,dc,nrefc,nvar,ndord)
      implicit real*8 (a-h,o-z)
      dimension x(3),cc(*),dc(nvar,*),dfdx(10)
      external f
c      write(*,*) 'dcoef sub nvar ====',nvar
      do 200 n=1,nvar
c      write(*,*) 'dcoef sub n ====',n
      call sdfdx(f,x,nrefc,dfdx,ndord,n,m)
      cc(n) = dfdx(1)
      do 100 i=1,m-1
      dc(n,i) = dfdx(i+1)
100   continue
200   continue
      return
      end
 
      subroutine sdfdx(f,x,n,dfdx,ndord,iv,m)
      implicit real*8 (a-h,o-z)
      common /coord/ coor(3),coora(27,3)
      dimension x(3),y(3),fx(27),dfdx(*)
      external f
      h =0.02
      h2 = h*2
c      write(*,*) 'h =',h,'     x =',x
      m = 0
      do 300 i=-1,1
      y(3) = x(3)+h*i
      do 200 j=-1,1
      y(2) = x(2)+h*j
      do 100 k=-1,1
      y(1) = x(1)+h*k
      m = m+1
      if ((n.eq.1) .and. (i*j.ne.0)) goto 100
      if ((n.eq.2) .and. (i.ne.0)) goto 100
      if (ndord.eq.1) then
      kk = k*k+j*j+i*i
      if (kk.gt.1) goto 100
      endif
      do 50 l=1,3
      coor(l)=coora(m,l)
50    continue
      fx(m) = f(y,iv)
c      write(*,'(1x,i3,4f10.4)') m,y,fx(m)
100   continue
200   continue
300   continue
400   continue
c      write(*,*) 'iv,n,m =',iv,n,m,'     fx ='
c      write(*,'(9f8.3)') (fx(i),i=1,m)
      mm = m/2+1
      dfdx(1) = fx(mm)
      m = n+1
      do 500 i=1,n
      j = 3**(i-1)
      dfdx(i+1) = (fx(mm+j)-fx(mm-j))/h2
500   continue
      if (ndord.eq.1) goto 1000
      do 900 i=1,n
      do 800 j=i,n
      i1 = 3**(i-1)
      j1 = 3**(j-1)
      m = m+1
      if (i.eq.j) then
      dfdx(m) = (fx(mm+i1)+fx(mm-i1)-fx(mm)*2)/h/h
      else
      dfdx(m) = (fx(mm+j1+i1)+fx(mm-j1-i1)
     &          -fx(mm+j1-1)-fx(mm-i1-1))/h2/h2
      endif
c     write(*,*) 'i,j,i1,j1,m,dfdx =',i,j,i1,j1,m,dfdx(m)
800   continue
900   continue
1000  continue
c      write(*,*) 'n,m =',n,m,'      dfdx ='
c      write(*,*) (dfdx(i),i=1,m)
      return
      end
 
      subroutine dcoor(f,x,cc,dc,nrefc,nvar,ndord)
      implicit real*8 (a-h,o-z)
      common /coord/ coor(3),coora(27,3)
      dimension x(3),cc(*),dc(nvar,*),dfdx(10)
      external f
c      write(*,*) 'dcoef sub nvar ====',nvar
      do 200 n=1,nvar
c      write(*,*) 'dcoef sub n ====',n
      call sdcdx(f,x,nrefc,dfdx,coora(1,n),ndord,n,m)
      cc(n) = dfdx(1)
      do 100 i=1,m-1
      dc(n,i) = dfdx(i+1)
100   continue
200   continue
c      write(*,*) 'coora =',coora
      return
      end
 
      subroutine sdcdx(f,x,n,dfdx,fx,ndord,iv,m)
      implicit real*8 (a-h,o-z)
      dimension x(3),y(3),fx(27),dfdx(*)
      external f
      h =0.02
      h2 = h*2
c      write(*,*) 'h =',h,'     x =',x
      m = 0
      do 300 i=-1,1
      y(3) = x(3)+h*i
      do 200 j=-1,1
      y(2) = x(2)+h*j
      do 100 k=-1,1
      y(1) = x(1)+h*k
      m = m+1
      if ((n.eq.1) .and. (i*j.ne.0)) goto 100
      if ((n.eq.2) .and. (i.ne.0)) goto 100
      if (ndord.eq.1) then
      kk = k*k+j*j+i*i
      if (kk.gt.1) goto 100
      endif
      fx(m) = f(y,iv)
c      write(*,'(1x,i3,4f10.4)') m,y,fx(m)
100   continue
200   continue
300   continue
400   continue
c      write(*,*) 'iv,n,m =',iv,n,m,'     fx ='
c      write(*,'(9f8.3)') (fx(i),i=1,m)
      mm = m/2+1
      dfdx(1) = fx(mm)
      m = n+1
      do 500 i=1,n
      j = 3**(i-1)
      dfdx(i+1) = (fx(mm+j)-fx(mm-j))/h2
500   continue
      if (ndord.eq.1) goto 1000
      do 900 i=1,n
      do 800 j=i,n
      i1 = 3**(i-1)
      j1 = 3**(j-1)
      m = m+1
      if (i.eq.j) then
      dfdx(m) = (fx(mm+i1)+fx(mm-i1)-fx(mm)*2)/h/h
      else
      dfdx(m) = (fx(mm+j1+i1)+fx(mm-j1-i1)
     &          -fx(mm+j1-1)-fx(mm-i1-1))/h2/h2
      endif
c     write(*,*) 'i,j,i1,j1,m,dfdx =',i,j,i1,j1,m,dfdx(m)
800   continue
900   continue
1000  continue
c      write(*,*) 'n,m =',n,m,'      dfdx ='
c      write(*,*) (dfdx(i),i=1,m)
      return
      end
 
      SUBROUTINE MSTRESS6(KDGOF,ESTR,EMSTR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ESTR(KDGOF),EMSTR(KDGOF)
      a1=estr(1)+estr(2)+estr(3)
      a2=estr(1)*estr(2)+estr(2)*estr(3)+estr(3)*estr(1)
      a2=a2-estr(4)**2-estr(5)**2-estr(6)**2
      a3=estr(1)*estr(2)*estr(3)
      a3=a3+2*estr(4)*estr(5)*estr(6)
      a3=a3-estr(1)*estr(4)**2-estr(2)*estr(5)**2
      a3=a3-estr(3)*estr(6)**2
      emstr(1)=a1 
      emstr(2)=a2
      emstr(3)=a3 
      pi=4*datan(1.d0)
      p=a2-a1*a1/3.
      q=a1*a2/3.-2.*a1**3/27.-a3
      dd=(q*0.5)**2+(p/3.)**3
      if(dd.gt.0)  dd=-dd
      if(q.eq.0) then
      ang=pi/2
      else
      ang=datan2(2.*dsqrt(-dd),-q)
      end if
      if(p.gt.0) p=-p
      aa=2.*dsqrt(-p/3.)
      emstr(4)=aa*dcos(ang/3.)
      emstr(5)=aa*dcos((2*pi+ang)/3.)
      emstr(6)=aa*dcos((2*pi-ang)/3.)
      emstr(4)=emstr(4)+a1/3
      emstr(5)=emstr(5)+a1/3
      emstr(6)=emstr(6)+a1/3
      if (emstr(4).lt.emstr(5)) then
      etemp=emstr(4)
      emstr(4)=emstr(5)
      emstr(5)=etemp
      end if
      if (emstr(5).lt.emstr(6)) then
      etemp=emstr(5)
      emstr(5)=emstr(6)
      emstr(6)=etemp
      end if
      if (emstr(4).lt.emstr(5)) then
      etemp=emstr(4)
      emstr(4)=emstr(5)
      emstr(5)=etemp
      end if 
      return
      END

      SUBROUTINE MSTRESS3(KDGOF,ESTR,EMSTR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ESTR(KDGOF),EMSTR(KDGOF)
      emstr(1)=estr(1)+estr(2)
      coefb=-estr(1)-estr(2)
      coefc=estr(1)*estr(2)-estr(3)**2
      emstr(2)=((-coefb)+dsqrt(coefb**2-4*coefc))/2
      emstr(3)=((-coefb)-dsqrt(coefb**2-4*coefc))/2
      if (emstr(2).lt.emstr(3)) then
      etemp=emstr(2)
      emstr(2)=emstr(3)
      emstr(3)=etemp
      endif
      return
      END

      SUBROUTINE MSTRESS4(KDGOF,ESTR,EMSTR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ESTR(KDGOF),EMSTR(KDGOF)
      emstr(1)=estr(1)+estr(2)+estr(3)
      emstr(2)=estr(2)
      coefb=-estr(1)-estr(3)
      coefc=estr(1)*estr(3)-estr(4)**2
      emstr(3)=((-coefb)+dsqrt(coefb**2-4*coefc))/2
      emstr(4)=((-coefb)-dsqrt(coefb**2-4*coefc))/2
      if (emstr(2).lt.emstr(3)) then
      etemp=emstr(2)
      emstr(2)=emstr(3)
      emstr(3)=etemp
      endif
      if (emstr(3).lt.emstr(4)) then
      etemp=emstr(3)
      emstr(3)=emstr(4)
      emstr(4)=etemp
      endif
      if (emstr(2).lt.emstr(3)) then
      etemp=emstr(2)
      emstr(2)=emstr(3)
      emstr(3)=etemp
      endif
      return
      END
