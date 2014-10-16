      implicit real*8 (a-h,o-z)
      TMAX=1.D1
      DT=1.D0
      TIME=0.D0
      IT=1
      OPEN(1,FILE='time',FORM='UNFORMATTED')
      WRITE(1) TMAX,DT,TIME,IT
      CLOSE(1)
      END

