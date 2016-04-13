c-----------------------------------------------------------------------
      module particle
      implicit none
      real*8 mx
      integer ch,istat

      contains

      real*8 function intgrnd(ener)
      real*8 ener,dshayield
      intgrnd=dshayield(mx,ener,ch,151,istat)
      return
      end function intgrnd

      end module particle
c-----------------------------------------------------------------------