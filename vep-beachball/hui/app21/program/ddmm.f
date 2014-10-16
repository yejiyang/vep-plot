      subroutine masterp
c      call Mazpartition
      call Mazsendpart
      print *,'send azsendpart    ok...................'
      call Mmsazrecvpart
      print *,'msazrecvpart       ok...................'
      return
      end

      subroutine masterp1 
      call Mnzrecvdip
      print *,'nzrecvdip          ok...................'
c      return
      call Mnzrecvstr
      print *,'nzrecvstr          ok...................'
      return
      end
