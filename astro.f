c-----------------------------------------------------------------------
      module astro
      use constants
      use particle
      use cluster_params
      implicit none
c-----declare global variables:
      real*8 z                  !redshift
      real*8 bmu                !magnetic field in microGauss      
      real*8 rkpc               !radius for ROI in kpc
      real*8 :: bsyn0=0.0253d0  !GeV/s
      real*8 :: bic0=0.265d0    !GeV/s
c-----------------------------------------------------------------------
      contains
c-----------------------------------------------------------------------
c-----functions      
c-----------------------------------------------------------------------
c-----function that defines the distance as a function of redshift integral
      real*8 function distintegrand(ztemp)
      real*8 ztemp
      distintegrand=Mpccm*clight/(h0*sqrt(omegaM*(1+ztemp)**3+omegaL))
      return
      end function distintegrand
c-----------------------------------------------------------------------
c-----distance as a function of redshift:      
      real*8 function dist(ztemp)
      real*8 ztemp
c-----integration parameters
      integer ier,last,neval
      integer :: limit=5000
      integer :: iord(5000)
      real*8 :: alist(5000),blist(5000),elist(5000),rlist(5000)
      real*8 abserr,result
      real*8 :: epsabs=1.0d-20
      real*8 :: epsrel=1.0d-3
c-----function definition
      call dqagse(distintegrand,1.d-16,ztemp,
     &     epsabs,epsrel,limit,
     &     result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      dist=result
      return
      end function dist
c-----------------------------------------------------------------------
c-----equilibrium spectrum; requires intgrnd (mod: particle) and bloss
      real*8 function dndeeq(ener)
      real*8 ener
c-----integration parameters
      integer ier,last,neval
      integer :: limit=5000
      integer :: iord(5000)
      real*8 :: alist(5000),blist(5000),elist(5000),rlist(5000)
      real*8 abserr,result
      real*8 :: epsabs=1.0d-20
      real*8 :: epsrel=1.0d-3
c------function definition
      call dqagse(intgrnd,ener,mx,
     &     epsabs,epsrel,limit,
     &     result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      dndeeq=result/bloss(ener)
      return
      end function dndeeq
c-----------------------------------------------------------------------
c-----energy loss function: IC,synch,coulomb,brem losses 
c-----expression from Colafrancesco+2006
      real*8 function bloss(ener)
      real*8 ener
      bloss=(bic0*(1+z)**4+bmu**2*bsyn0)*ener**2
     &     +nele*6.13d0*(1.d0+dlog((ener/me)/nele)/75.d0)
     &     +nele*1.51d0*(dlog((ener/me))+0.36d0)*ener
      
      bloss=bloss*1.d-16
      return
      end function bloss
c-----------------------------------------------------------------------
c-----dark matter profile (NFW)
      real*8 function rhox(r)
      real*8 r
      rhox=rhos/((r/rs)*(1+r/rs)**2)
      return
      end
c-----------------------------------------------------------------------
      end module astro
c-----------------------------------------------------------------------