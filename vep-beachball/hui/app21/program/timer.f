c*********************************************************************c
c FILE: timer.f                                                       c
c Author: Huai Zhang                                                  c
c DESCRIPTION:                                                        c
c*********************************************************************c
c                                                                     c
c       case itimer :                                                 c
c       1 i/o timer                                                   c  
c       2 timer for element computing                                 c
c       3 timer for reverse of gross A metrix                         c
c       4 timer for all communications                                c
c       5 timer for solver(iteration)time consuming                   c  
c       6 timer for usr time consuming                                c  
c       7 timer for system time consuming                             c
c       8 timer for all time consuming                                c
c       9 total time consuming                                        c
c       default error                                                 c
c                                                                     c
c       case istatus:                                                 c
c       1 timer begin                                                 c
c       2 timer end                                                   c
c       default error                                                 c
c                                                                     c
c*********************************************************************c
       subroutine timer(itimer,istatus)
       implicit real*8 (a-h,o-z)
       common /timer1/timio(3),timges(3),timar(3),timcmn(3),timitr(3)
       common /timer2/timusr(3),timsys(3),timall(3),timtot(3)

       dimension nwtimer(2),systimer(2)
       dimension nwtimer2(2),systimer2(2)
c
c      get wall timer to record present time 
c
        call fwtiming(nwtimer)
c
c      caculate the real wall time
c
       wtimer = nwtimer(1) + nwtimer(2) * 1.0e-6
c
c       call chgsec(wtimer,nday,nhour,nmin,sec) 
c       write(*,*) 'current wall time is=========',wtimer
c
c
c      get system timer to record present system and user time consuming
c
        call futiming(systimer)
c
c      caculte the user time consuming (user time and system time)
c      
c      
       usertimer = systimer(1)
       systemtimer = systimer(2)
c
c       write(*,*) 'user timer is=========',usertimer
c       write(*,*) 'system timer is=========',systemtimer
c
        select case(itimer) 
        case(1) 
          select case (istatus)
          case(1) 
          call fwtiming(nwtimer)
          wtimer = nwtimer(1) + nwtimer(2) * 1.0e-6
          timio(2) = wtimer
          case(2)
          call fwtiming(nwtimer)
          wtimer = nwtimer(1) + nwtimer(2) * 1.0e-6
          timio(3) = wtimer
          dwtimer  = timio(3) - timio(2)
          if (dwtimer.lt.0.d0) then
          write (*,*) 'fatal error when get timer of I/O'
          call endjob(ierr)
          end if
          timio(1) = timio(1) + dwtimer
          case default
          write(*,*) 'wrong input parameters for I/O timer'
          call endjob(ierr)
          end select !istatus
        case(2)
          select case (istatus)
          case(1) 
          call fwtiming(nwtimer)
          wtimer = nwtimer(1) + nwtimer(2) * 1.0e-6
          timges(2) = wtimer
          case(2)
          call fwtiming(nwtimer)
          wtimer = nwtimer(1) + nwtimer(2) * 1.0e-6
          timges(3) = wtimer
          dwtimer  = timges(3) - timges(2)
          if (dwtimer.lt.0.d0) then
          write (*,*) 'fatal error in timer of element computing'
          call endjob(ierr)
          end if
          timges(1) = timges(1) + dwtimer
          case default
          write(*,*) 'wrong input parameters for element compute timer'
          call endjob(ierr)
          end select !istatus
        case(3)
          select case (istatus)
          case(1) 
          call fwtiming(nwtimer)
          wtimer = nwtimer(1) + nwtimer(2) * 1.0e-6
          timar(2) = wtimer
          case(2)
          call fwtiming(nwtimer)
          wtimer = nwtimer(1) + nwtimer(2) * 1.0e-6
          timar(3) = wtimer
          dwtimer  = timar(3) - timar(2)
          if (dwtimer.lt.0.d0) then
          write (*,*) 'fatal error in timer of A reverse'
          call endjob(ierr)
          end if
          timar(1) = timar(1) + dwtimer
          case default
          write(*,*) 'wrong input parameters for A reverse timer'
          call endjob(ierr)
          end select !istatus
        case(4)
          select case (istatus)
          case(1) 
          call fwtiming(nwtimer)
          wtimer = nwtimer(1) + nwtimer(2) * 1.0e-6
          timcmn(2) = wtimer
          case(2)
          call fwtiming(nwtimer)
          wtimer = nwtimer(1) + nwtimer(2) * 1.0e-6
          timcmn(3) = wtimer
          dwtimer  = timcmn(3) - timcmn(2)
          if (dwtimer.lt.0.d0) then
          write (*,*) 'fatal error in timer of communication'
          call endjob(ierr)
          end if
          timcmn(1) = timcmn(1) + dwtimer
          case default
          write(*,*) 'wrong input parameters for communication timer'
          call endjob(ierr)
          end select !istatus
        case(5)
          select case (istatus)
          case(1) 
          call fwtiming(nwtimer)
          wtimer = nwtimer(1) + nwtimer(2) * 1.0e-6
          timitr(2) = wtimer
          case(2)
          call fwtiming(nwtimer)
          wtimer = nwtimer(1) + nwtimer(2) * 1.0e-6
          timitr(3) = wtimer
          dwtimer  = timitr(3) - timitr(2)
          if (dwtimer.lt.0.d0) then
          write (*,*) 'fatal error in timer of iteration solver'
          call endjob(ierr)
          end if
          timitr(1) = timitr(1) + dwtimer
          case default
          write(*,*) 'wrong input parameters for solver timer'
          call endjob(ierr)
          end select !istatus
        case(6)
          select case (istatus)
          case(1) 
          call futiming(systimer)
          usertimer = systimer(1)
          timusr(1) = timusr(1) + usertimer
          case default
          write(*,*) 'wrong input parameters for user timer'
          call endjob(ierr)
          end select !istatus
        case(7)
          select case (istatus)
          case(1) 
          call futiming(systimer)
          systemtimer = systimer(2)
          timsys(1) = timsys(1) + systemtimer
          case default
          write(*,*) 'wrong input parameters for system timer'
          call endjob(ierr)
          end select !istatus
       case(8)
          select case (istatus)
          case(1) 
          timall(1) = timio(1)+timges(1)+timar(1)+timcmn(1)+timitr(1)
          timall(2) = timusr(1)
          timall(3) = timsys(1)
          case default
          write(*,*) 'wrong input parameters for timer summering'
          call endjob(ierr)
          end select !istatus
       case(9)
          select case (istatus)
          case(1) 
          call fwtiming(nwtimer)
          wtimer = nwtimer(1) + nwtimer(2) * 1.0e-6
          timtot(2) = wtimer
          case(2)
          call fwtiming(nwtimer)
          wtimer = nwtimer(1) + nwtimer(2) * 1.0e-6
          timtot(3) = wtimer
          dwtimer  = timtot(3) - timtot(2)
          if (dwtimer.lt.0.d0) then
          write (*,*) 'fatal error in timer of total time get'
          call endjob(ierr)
          end if
          timtot(1) = timtot(1) + dwtimer
          case default
          write(*,*) 'wrong input parameters for total time get'
          call endjob(ierr)
          end select !istatus
       case default
          write(*,*) 'Wrong input timer parameters or timer selection'
          call endjob(ierr)
       end select   
        return
        end

       subroutine initimer(ierr)
       implicit real*8 (a-h,o-z) 
       common /timer1/timio(3),timges(3),timar(3),timcmn(3),timitr(3)
       common /timer2/timusr(3),timsys(3),timall(3),timtot(3)
c
c      start initial the timers
c       
       do i = 1,3
       timio(i) = 0.d0
       timges(i) = 0.d0
       timar(i) = 0.d0
       timcmn(i) = 0.d0
       timitr(i) = 0.d0
       timusr(i) = 0.d0
       timsys(i) = 0.d0
       timall(i) = 0.d0
       timtot(i) = 0.d0
       end do 
       ierr = 0
       return
       end
       

       subroutine chgsecond(s,nday,nhour,nmin,sec)
       implicit real*8 (a-h,o-z) 
c
        sec = 0.d0
        nday = 0
        nhour  = 0
        nmin = 0
c
        nday    =  S/86400.d0
        nhour   =  (S - nday*86400.d0)/3600.d0
        sec    =  (S - nday*86400.d0 - nhour*3600)/60.d0
        nmin    =  sec
        sec    =  (sec - nmin) * 60.d0
C
        return
        end

c*********************************************************************c
c FILE: timereport.f                                                  c
c Author: Huai Zhang                                                  c
c DESCRIPTION:                                                        c
c*********************************************************************c
c   Summerizing and report tht time consumming for each distributed   c
c   parallel computer node and all the nodes add togetherc.           c
c   Data communication is needed when write the report in the Master  c
c   node.                                                             c
c                                                                     c
c*********************************************************************c

       subroutine timereport(myid,myblk,numblk,ioport,isreport)
       implicit real*8 (a-h,o-z)
       common /timer1/timio(3),timges(3),timar(3),timcmn(3),timitr(3)
       common /timer2/timusr(3),timsys(3),timall(3),timtot(3)
       dimension timera(9,514),timeriblk(9)
       dimension timermax(9),timermin(9),balance(9)
       dimension nodmax(9),nodmin(9)
       character*80 path
c
       master = 0
       do i = 1,9
       nodmax(i) = 0
       nodmin(i) = 0
       timermax(i) = 0.d0
       timermin(i) = 0.d0
       balance(i) = 0.d0
       timeriblk(i) = 0.d0
       end do
c
c
c      data error check
c
       if((myid.lt.master).or.(myid.gt.numblk-1)) then
       write(*,*) 'Node rank input error when timer report'
       end if
       if((myblk.lt.master).or.(myblk.gt.numblk)) then
       write(*,*) 'Block number input error when timer report'
       end if
c       
c      first collect all the timer data from all the slave nodes
c
c
c      slaves send all the timer data to master node 
c

       if((myid.gt.master).and.(myid.le.numblk-1)) then
       timeriblk(1) = timio(1)
       timeriblk(2) = timges(1)
       timeriblk(3) = timar(1)
       timeriblk(4) = timcmn(1)
       timeriblk(5) = timitr(1)
       timeriblk(6) = timusr(1)
       timeriblk(7) = timsys(1)
       timeriblk(8) = timall(1)
       timeriblk(9) = timtot(1)
       call sendar(master,myblk,timeriblk,9)
       end if
c
c      master node collect all the timer data from slave nodes
c
       if(myid.eq.master) then
       do iblk = 1,numblk-1
       call recvar(master,iblk,timeriblk,9 )
       do i = 1,9
       timera(i,iblk)=timeriblk(i)
       end do
       end do 
c        
c      write out the timer report        
c        
       select case(isreport)
       case default
       write(*,*) ' Error input of isreport, 0 or 1, check please    ' 
       case(0)
       write(*,*) 'You have disable the timer report function'
       write(*,*) 'This function is very important for your parallel'//
     +          ' computing analysys, and if you want to turn on it,'//
     +            'just isreport=1, please, good luck '
          
       case(1) 
       nn = ioport
       write(nn,*) '**************************************************'
       write(nn,*) '* Time consuming of all nodes list as following: *'
       write(nn,*) '**************************************************'
       write(nn,*) '* Discriptions of table title                    *'
       write(nn,*) '* I/O Time consuming of I/O                      *'
       write(nn,*) '* ges Time consuming of element computing        *'
       write(nn,*) '* ar  Time consuming of reverse of gross A metrix*'
       write(nn,*) '* cmn Time consuming of all data communications  *'
       write(nn,*) '* itr Time consuming of parallel solvers         *'
       write(nn,*) '* usr Time consuming of usr time consuming       *'
       write(nn,*) '* sys Time consuming of system time consuming    *'
       write(nn,*) '* all Time consuming of total time consuming     *'
       write(nn,*) '* tot total time consuming of each node          *'
       write(nn,*) '**************************************************'
       timeriblk(1) = timio(1)
       timeriblk(2) = timges(1)
       timeriblk(3) = timar(1)
       timeriblk(4) = timcmn(1)
       timeriblk(5) = timitr(1)
       timeriblk(6) = timusr(1)
       timeriblk(7) = timsys(1)
       timeriblk(8) = timall(1)
       timeriblk(9) = timtot(1)

       write(nn,*) '*           Master node                          *'
       write(nn,*) 'rank   I/O     ges     ar     cmn    itr     usr  '
     +             //'   sys      all     tot'
       write(nn,9) master,(timeriblk(i),i = 1,9)
       write(nn,*) '*           Slave nodes                          *'
       write(nn,*) 'iblk   I/O     ges     ar     cmn    itr     usr  '
     +             //'   sys      all     tot'
       do iblk = 1,numblk-1
       write(nn,9) iblk,(timera(i,iblk),i = 1,9)
       end do
       write(nn,*) '**************************************************'
c
c      Now, analysis the loadblance amoung slave nodes and master node
c
       do i = 1,9
       balance(i) = 0.d0
       timermax(i) = timera(i,1)
       nodmax(i) = 1
       timermin(i) = timera(i,1)
       nodmin(i) = 1
       do iblk = 2,numblk-1
       if(timera(i,iblk).gt.timermax(i)) then
       timermax(i) = timera(i,iblk)
       nodmax(i) = iblk
       end if
       if(timera(i,iblk).lt.timermin(i)) then
       timermin(i) = timera(i,iblk)
       nodmin(i) = iblk
       end if
       end do
       if(timermin(i).lt.1.0e-4) then
       balance(i) = 0.d0
       else
       balance(i) = timermax(i)/timermin(i)
       end if
       end do
       if(timeriblk(8).gt.timermax(8)) then
       timermax(8) = timeriblk(8)
       nodmax(8) = master
       end if
       if(timeriblk(8).lt.timermin(8)) then
       timermin(8) = timeriblk(8)
       nodmin(8) = master
       end if
       if(timeriblk(9).gt.timermax(9)) then
       timermax(9) = timeriblk(9)
       nodmax(9) = master
       end if
       if(timeriblk(9).lt.timermin(9)) then
       timermin(9) = timeriblk(9)
       nodmin(9) = master
       end if
c
       write(nn,*)   
       write(nn,*) '**************************************************'
       write(nn,*) '*******   Analysys load imbalance    *************'
       write(nn,*) '**************************************************'
       write(nn,*) ' Type  maxtime  rank   mintime   rank  imbalance  '
       write(nn,8) ' I/O ', timermax(1), nodmax(1), timermin(1),
     +               nodmin(1), balance(1)
       write(nn,8) ' ges ', timermax(2), nodmax(2), timermin(2),
     +               nodmin(2), balance(2)
       write(nn,8) ' ar  ', timermax(3), nodmax(3), timermin(3),
     +               nodmin(3), balance(3)
       write(nn,8) ' cmn ', timermax(4), nodmax(4), timermin(4),
     +               nodmin(4), balance(4)
       write(nn,8) ' itr ', timermax(5), nodmax(5), timermin(5),
     +               nodmin(5), balance(5)
       write(nn,8) ' usr ', timermax(6), nodmax(6), timermin(6),
     +               nodmin(6), balance(6)
       write(nn,8) ' sys ', timermax(7), nodmax(7), timermin(7),
     +               nodmin(7), balance(7)
       write(nn,8) ' all ', timermax(8), nodmax(8), timermin(8),
     +               nodmin(8), balance(8)
       write(nn,8) ' TOT ', timermax(9), nodmax(9), timermin(9),
     +               nodmin(9), balance(9)
       write(nn,*) '**************************************************'
       write(nn,*) 'Total time elipse =========',timermax(9)
       write(nn,*) '**************************************************'
       end select
c
       do i = 1,80
       path(i:i) = ' '
       end do
c        call fgetcwd(lenth, path)
c        write(*,*) lenth, path
       end if
c
8      format(a,3x,e7.2,1x,i5,4x,e7.2,1x,i5,4x,e7.2)
9      format(i3,1x,9(1x,e7.2))        
       return
       end 
