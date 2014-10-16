       program main
c       subroutine bdseg(ierr)
       implicit real*8 (a-h,o-z)
       dimension xb(3000),yb(3000),ibse1(3000),ibse2(3000),
     * markbeg(3000),idx(3000),idy(3000),idm(3000),
     * x(6000),y(6000),no1(6000),no2(6000),matb(6000),matn(6000)
       dimension elemx(30000),elemy(30000),index(6000)

c..................................................................
c        read(*,*) boundary line file
c.................................................................
          iprint = 0
c          Iprint = 1
	  ierr = 1

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
c
c
	  if(iprint.eq.1) then
	  write(*,*)  numb
        do i=1,numb
        write(*,*) xb(i), yb(i)
        end do
        write(*,*) nbseg
        do i=1,nbseg
        write(*,*) ibse1(i), ibse2(i), markbeg(i)
        end do
c        pause
	  end if

c.................................................................
         open(11,file='bd.cor',form='formatted',status='unknown')
         open(12,file='bd.elm',form='formatted',status='unknown')
         read(11,*) numnodb 
         do i=1,numnodb 
         read(11,*) x(i), y(i)
         end do
         read(12,*) nemb
         do i=1,nemb
         read(12,*) no1(i), no2(i)
         end do
         close(11)
         close(12)
         if(iprint.eq.1) then
         write(*,*) numnodb 
         do i=1,numnodb 
         write(*,*) x(i), y(i)
         end do
         write(*,*) nemb
         do i=1,nemb
         write(*,*) no1(i), no2(i)
         end do
c         pause
	   end if
c..................................................................
c        read(*,*) all mesh file and tran the No of nodes on
c        boundaries to all node No of mesh element
c.................................................................
         open(11,file='mesh.cor',form='formatted',status='unknown')
         read(11,*) nelem
         do i=1,nelem
         read(11,*) elemx(i),elemy(i)
         enddo
         close(11)
         do i=1,numnodb
         do j=1,nelem
         error=(x(i)-elemx(j))*(x(i)-elemx(j))
     +  +(y(i)-elemy(j))*(y(i)-elemy(j))
         if (error .lt. 1e-7) then
         index(i)=j
         exit
         endif
         enddo
         enddo

c..................................................................
c        begin boundary processing
c.................................................................
       open(21,file='all-boundaries',form='formatted',status='unknown')
       open(22,file='boundaries',form='formatted',status='unknown')
c
	 write(21,*) nbseg
       do k = 1,nbseg
c
	   do i = 1, numnodb
	   matn(i) = 0
	   end do

         do i = 1, nemb
	   matb(i) = 0
	   end do
c	   
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
c	  
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
              matn(i)=1
          endif
          if((x2.ge.x0min).and.(x2.le.x0max)) then
          if((y2.ge.y0min).and.(y2.le.y0max)) then
	  matn(i) = 1
	  end if
	  end if 
	  end if
	  end do

        do i=1,nemb
        if((matn(no1(i)).eq.1).and.(matn(no2(i)).eq.1)) then
	  matb(i) = 1
	  end if
	  end do
c
	  nb = 0
        do i=1,nemb
        if(matb(i).eq.1) then
	  nb = nb + 1
	  matn((nb-1)*3+1) = i
	  matn((nb-1)*3+2) = index(no1(i))
	  matn((nb-1)*3+3) = index(no2(i))
        end if
        end do
c
c       write out boundary information 
c
         write(21,*) k,markbeg(k),nb,idx(k),idy(k),idm(k)
         if(markbeg(k).eq.-1)
     #   write(22,*) k,markbeg(k),nb,idx(k),idy(k),idm(k)
	   do i = 1,nb
         write(21,*)i, (matn((i-1)*3+j),j=2,3)
         if(markbeg(k).eq.-1)
     #   write(22,*)i, (matn((i-1)*3+j),j=2,3)
         end do
c
c	 boundary line loop finished
c   
       end do !k loop
       close(21)
	 ierr = 0
c         return
	 end
	   
