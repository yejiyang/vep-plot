       program gid2prg
       implicit real*8 (a-h,o-z)
       parameter(maxb=3000,maxnode=35000,maxcell=600)
       dimension x(maxnode),y(maxnode),
     *  itri(4,maxnode*2),mat(maxnode*2)
       dimension xb(3000),yb(3000),ibse1(3000),ibse2(3000),
     * markbeg(3000),idx(3000),idy(3000),idm(3000),
     * no1(6000),no2(6000),nob(1000)
       dimension xx(maxnode),yy(maxnode),ibse3(6000),ibse4(6000),
     * nnnode(maxnode),ibs1(6000),ibs2(6000)
     
     
c...............................................................
         
c         open(21,file='test.flavia.msh',status='unknown',
c     *   form='formatted')
c         read(21,*)
c         read(21,*)
c         nnode=1210
c         ntri=1151
c         do i=1,nnode
c         read(21,*) k,x(i),y(i)
c         enddo
c         read(21,*)
c         read(21,*)
c         do i=1,ntri
c         read(21,*) k,itri(1,i),itri(2,i),itri(3,i),itri(4,i),mat(i)
c         enddo
c         close(21)

         open(21,file='data1',status='unknown',form='formatted')
         read(21,*)nnode
         do i=1,nnode
         read(21,*) k,x(i),y(i)
         enddo
         read(21,*)ntri
         do i=1,ntri
         read(21,*) k,itri(1,i),itri(2,i),itri(3,i),itri(4,i),mat(i)
         enddo
         close(21)
         write(*,*)'node = ',nnode,'    nqudr =',ntri



c...............................................................
         open(20,file='mesh.prj',form='formatted',status='unknown')
         open(21,file='mesh.cor',form='formatted',status='unknown')
         open(22,file='mesh.elm',form='formatted',status='unknown')
         write(20,*) 'mesh'
         write(20,*) 0
         write(20,*) 0
         write(20,*) 0
         write(21,*) nnode , 2
         do i=1,nnode
         write(21,1000) x(i), y(i)
         enddo
         write(22,2000) ntri, 5
        do i=1,ntri
        write(22,2000) itri(1,i),itri(2,i),itri(3,i),itri(4,i),mat(i)
         enddo
         close(20)
         close(21)
         close(22)
c
c        find the node on the boundaries
c

        iprint = 0
c        iprint = 1
        ierr = 1
         
c...............................................................
        open(10,file='t2.dat',form='formatted',status='unknown')
        read(10,*)  numb
        do i=1,numb
        read(10,*)ii, xb(i), yb(i)
        end do
        read(10,*) num0
        do i=1,num0
        read(10,*) 
        end do
        read(10,*) nbseg
        do i=1,nbseg
        read(10,*) ii,ibse1(i), ibse2(i), markbeg(i),
     #   idx(i),idy(i),idm(i)
        end do
        close(10)

        if(iprint.eq.1) then
        write(*,*)  'numb=', numb
        do i=1,numb
c        write(*,*) xb(i), yb(i)
        end do
        write(*,*) 'nbseg=',nbseg
        do i=1,nbseg
c        write(*,*) ibse1(i), ibse2(i), markbeg(i)
        end do
c        pause
        end if

        numnodb=nnode
        write(*,*)  'numnodb=', nnode
        
c       numb2: the total nodes on boundaries
c       numb3: the total lines on boundaries
        numb2=0
        numb3=0

        do k = 1,nbseg
        
        x0 = xb(ibse1(k))
        y0 = yb(ibse1(k))
        x1 = xb(ibse2(k))
        y1 = yb(ibse2(k))
        x0min=x0
        x0max=x0
        y0min=y0
        y0max=y0
        if(x0.gt.x1) then
        x0max = x0
        x0min = x1
        else
        x0max = x1
        x0min = x0
        end if
c
        if(y0.gt.y1) then
        y0max = y0
        y0min = y1
        else
        y0max = y1
        y0min = y0
        end if
        
        numb1=0
c       numb1: the total nodes on echo boundary 
        
	    do i = 1,numnodb
	    x2 = x(i)
	    y2 = y(i)
	    ddet = dabs((y1-y0)*(x2-x0)-(y2-y0)*(x1-x0))
	    if(ddet.lt.1.0e-4) then
          dx0=abs(x2-x0)
          dx1=abs(x2-x1)
          dy0=abs(y2-y0)
          dy1=abs(y2-y1)
          if(((dx0.le.1e-3).and.(dy0.le.1e-3)).or.
     *      ((dx1.le.1e-3).and.(dy1.le.1e-3))) then
              numb1=numb1+1
              nob(numb1) = i
          endif
             
          if((x2.ge.x0min).and.(x2.le.x0max)) then
          if((y2.ge.y0min).and.(y2.le.y0max)) then
          numb1=numb1+1
	      nob(numb1) = i
	    end if
	    end if 
	    end if
	    end do
	    
        do i=1,numb1
        do j=1,i
        x0=x(nob(i))
        x1=x(nob(j))        
        if (x0 .lt. x1) then
        ttt=nob(i)
        nob(i)=nob(j)
        nob(j)=ttt
        endif
        enddo
        enddo        
        
        do i=1,numb1-1
        no1(i)=i
        no2(i)=i+1
        enddo

        do i=1,numb1
        numb2=numb2+1
        xx(numb2)=x(nob(i))
        yy(numb2)=y(nob(i))
        enddo
        
        do i=1,numb1-1
        if (k .eq. 1) then
        ibse3(i+numb3)=no1(i)
        ibse4(i+numb3)=no2(i)
        else        
        ibse3(i+numb3)=no1(i)+numb2-numb1
        ibse4(i+numb3)=no2(i)+numb2-numb1
        endif
        enddo
        numb3=numb3+numb1-1
        
	    
c        iprint=1
        if(iprint.eq.1) then
	    write(*,*) k,numb1
	    do i=1,numb1
	    write(*,*) nob(i),x(nob(i)),y(nob(i))
	    enddo
        endif
        
        enddo
	    write(*,*) 'numb2=',numb2,'    numb3=',numb3
	    
	    do i=1,numb2
	    enddo
	    do i=1,numb3
	    enddo

c...............................................................
	    
c
c      merge the same dots
c      
      nnp=numb2
      numb=numb3
      call merge2(nnp,mm,xx,yy,nnnode)
      write(*,*) 'nnp=',nnp,'     mm=',mm

      do i=1,numb
      ibs1(i)=nnnode(ibse3(i))
      ibs2(i)=nnnode(ibse4(i))
      enddo
      
      do i=1,numb
      ibs1(i)=nnnode(ibse3(i))
      ibs2(i)=nnnode(ibse4(i))
      enddo

      k=0
      do i=1, numb      
      if (ibs1(i) .ne. ibs2(i)) then
      k=k+1      
      ibse3(k)=ibs1(i)
      ibse4(k)=ibs2(i)
      endif
      enddo      
      numb=k
      write(*,*) 'numb=',numb
      
c
c      write the result files
c      
      open(21,file='bd.prj',status='unknown',form='formatted')
      write(21,*)'bd'
      write(21,*)'0'
      write(21,*)'0'
      write(21,*)'0'
      close(21)

      open(21,file='bd.cor',status='unknown',form='formatted')
      write(21,*)mm,2
      do i=1,mm
      write(21,*)xx(i),yy(i)
      enddo
      close(21)
 
      open(21,file='bd.elm',status='unknown',form='formatted')
      write(21,*)numb,2
      do i=1,numb
      write(21,*)ibse3(i),ibse4(i)
      enddo
      close(21)
      


1000     format(2f15.5)
2000     format(14i10)      
1100     format(i5,20f9.4)
      end
      
      
       
      subroutine merge2(numnod,mm,x,y,np)
      implicit real*8 (a-h,o-z)
      double precision cn
      dimension np(numnod),x(numnod),y(numnod),index(2000)
c     write(*,*) 'x,y ='
c     write(*,*) x
c     write(*,*) y
      e = 1.e-5
      xmin = x(1)
      xmax = x(1)
      ymin = y(1)
      ymax = y(1)
      do 100 n=1,numnod
      np(n) = n
      if (x(n).lt.xmin) xmin = x(n)
      if (x(n).gt.xmax) xmax = x(n)
      if (y(n).lt.ymin) ymin = y(n)
      if (y(n).gt.ymax) ymax = y(n)
100   continue
      ex = (xmax-xmin)*e
      ey = (ymax-ymin)*e
      write(*,*) 'ex,ey =',ex,ey
      cn = numnod
      mk = dsqrt(cn)
      xe = (xmax-xmin)/mk
      do 550 ik=1,mk
      k = 0
      x0 = xmin+xe*(ik-1)-ex
      x1 = xmin+xe*ik+ex
      do 150 n=1,numnod
      if (x0 .gt. x(n)) goto 150
      if (x(n) .gt. x1) goto 150
      k = k+1
      index(k) = n
150   continue
      do 500 in=2,k
      n = index(in)
      do 400 im=1,in-1
      m = index(im)
      dx = x(m)-x(n)
      if (dx.lt.0.0) dx=-dx
      dy = y(m)-y(n)
      if (dy.lt.0.0) dy=-dy
      if (dx.lt.ex .and. dy.lt.ey) then
c     write(*,*) 'n,m =',n,m
c     write(*,*) 'dx,dy =',dx,dy
c     write(*,*) 'xn,yn,xm,ym =',x(n),y(n),x(m),y(m)
      np(n) = m
      goto 500
      endif
400   continue
500   continue
550   continue       
cc    write(*,*) 'np =',np
      m = 0
      do 600 n=1,numnod
      k = np(n)
      if (k.eq.n) then
      m = m+1
      np(n) = m
      x(m) = x(k)
      y(m) = y(k)
      else
      np(n) = np(k)
      endif
600   continue
      mm = m
cc    write(*,*) 'np =',np
      return
      end
