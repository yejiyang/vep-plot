        subroutine openf(iu,if,istatus)
      implicit real*8 (a-h,o-z)
        common /d2cd/ initi(100),initr(100),iu2fi(100),iu2fr(100),
     &  iu2fi0(100),iu2fr0(100),lasti,lastr
        if (istatus.eq.0) then
        if (initi(if).eq.0) initi(if) = lasti
        if (initr(if).eq.0) initr(if) = lastr
        endif
        iu2fi(iu) = initi(if)
        iu2fr(iu) = initr(if)
        iu2fi0(iu) = initi(if)
        iu2fr0(iu) = initr(if)
c      if (if.eq.13) then
c        write(*,*) 'initi(if),initr(if) =='
c        write(*,*)  initi(if),initr(if)
c      endif
   !     write(*,*) ' 11111       istatus,lasti =',istatus,lasti
   !     write(*,*) ' 22222       iu,if,initi(if) =',iu,if,initi(if)
   !    call findrunk(myid)
   !    if (myid.eq.1) then
   !    print *,(  initi(i),i=1,15)
   !    endif
        return
        end
 
        subroutine rewindf(iu)
      implicit real*8 (a-h,o-z)
        common /d2cd/ initi(100),initr(100),iu2fi(100),iu2fr(100),
     &  iu2fi0(100),iu2fr0(100),lasti,lastr
        iu2fi(iu) = iu2fi0(iu)
        iu2fr(iu) = iu2fr0(iu)
        return
        end
 
        subroutine initrwf(iu,iui,iur)
      implicit real*8 (a-h,o-z)
        common /d2cd/ initi(100),initr(100),iu2fi(100),iu2fr(100),
     &  iu2fi0(100),iu2fr0(100),lasti,lastr
        iui = iu2fi(iu)
        iur = iu2fr(iu)
   !     write(*,*) '   333333     iui,iur =',iui,iur
        return
        end
 
        subroutine endrwf(iu,iui,iur)
      implicit real*8 (a-h,o-z)
        common /d2cd/ initi(100),initr(100),iu2fi(100),iu2fr(100),
     &  iu2fi0(100),iu2fr0(100),lasti,lastr
        iu2fi(iu) = iui
        iu2fr(iu) = iur
        if (iui.gt.lasti) lasti = iui
        if (iur.gt.lastr) lastr = iur
   !     write(*,*) 'iu,iui,lasti,iu2fi(iu) = ',iu,iui,lasti,iu2fi(iu)
        return
        end
 
        subroutine initpool
      implicit real*8 (a-h,o-z)
        common /d2cd/ initi(100),initr(100),iu2fi(100),iu2fr(100),
     &  iu2fi0(100),iu2fr0(100),lasti,lastr
        lasti = 0
        lastr = 0
        do 100 i=1,100
        initi(i)=0
        initr(i)=0
100     continue
        return
        end
