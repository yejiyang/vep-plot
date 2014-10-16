      subroutine slavep(iblk)
      master = 0

      istop=-1
c      istop=0
      kend=1
      if (iblk.eq.master) then
      call Mazsendpart(istop)
c      print *,'send azsendpart    ok...................'
      call Mmsazrecvpart(istop)
c      print *,'msazrecvpart       ok...................'
      end if

      call Mazrecvpart(iblk,istop)
c      print *,'azrecvpart ok...................',iblk

      call Mstart(iblk)
c      print *,'start      ok...................',iblk
      call Mnzdmbsdiag(iblk)
c      print *,'nzdmbsdiag ok...................',iblk
1     continue
      Call Mmbft(iblk,istop,kend)
      Call Mbft(iblk,istop,kend)
c      print *,'bft        ok...................',iblk
      call Meddm(iblk)
c      print *,'eddm       ok...................',iblk
      call Mdmbs2dmsr(iblk)
c      print *,'dmbs2dmsr  ok...................',iblk
      call Mazsolv(iblk)
c      print *,'azoslv     ok...................',iblk
      call Muddm(iblk,istop,kend)
c      print *,'uddm       ok...................',iblk
      call Mesddm(iblk,istop,kend)
c      print *,'esddm      ok...................',iblk
      if (kend.ne.1) goto 1
c      call Mesddm(iblk,istop,kend)

      if (iblk.eq.master) then
      call Mnzrecvdip(istop)
c      print *,'nzrecvdip          ok...................'
      call Mnzrecvstr(istop)
c      print *,'nzrecvstr          ok...................'
      end if

      if (istop.le.0) goto 1

      end
