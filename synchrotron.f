c-----------------------------------------------------------------------
      module synchrotron
c-----------------------------------------------------------------------
      use astro
      use cluster_params
      implicit none
c-----set radio observation specific parameters
      real*8 Sin,Sout           !Sin=synch signal from DM, Sout=Radio UL (telescope-specific)
      real*8 :: nu=1.4!0.120        !observation frequency in GHz
c-----here add in how Sout is defined in terms of telescope parameters
      real*8 :: frms=1.d-5 !noise per beam in Jy
      real*8 :: thetaB=25.d0    !beam size in arcsec
      contains
c-----------------------------------------------------------------------
c-----synchrotron emission spectral function (from Colafrancesco+2006)
      real*8 function fff(x)
      real*8 x
      fff=1.25d0*x**(1.d0/3.d0)*dexp(-x)*(648.d0+x**2)**(1.d0/12.d0)
      return
      end function fff
c-----------------------------------------------------------------------
c-----synchrotron single electron power, averaged over isotropic pitch angles
      real*8 function psyn(ener)
      real*8 ener
c     external psynintegrand
c-----integration parameters
      integer ier,last,neval
      integer :: limit=5000
      integer :: iord(5000)
      real*8 :: alist(5000),blist(5000),elist(5000),rlist(5000)
      real*8 abserr,result
      real*8 :: epsabs=1.0d-20
      real*8 :: epsrel=1.0d-3
c-----function definition      
      call dqagse(psynintegrand,1.d-16,pi,
     &     epsabs,epsrel,limit,
     &     result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      
      psyn=result*2*dndeeq(ener) !factor of 2 for positrons + electrons
      
      contains
c-----synchrotron single electron power, integrand
      real*8 function psynintegrand(theta)
      real*8 theta,xxx,nu_em
      real*8 :: psyn0=1.46323d-25 !GeV/s/Hz
      real*8 :: x0=62.1881      !dimensionless constant
      nu_em=(1+z)*nu            !emitted frequency = (1+z)(observed frequency)
      xxx=x0*nu_em/(bmu*ener**2)
      psynintegrand=psyn0*bmu*0.5*(dsin(theta))**2
     &     *fff(xxx/dsin(theta))
      end function psynintegrand
      
      end function psyn
c-----------------------------------------------------------------------
      real*8 function jsyn(r)
      real*8 r
c-----integration parameters
      integer ier,last,neval
      integer :: limit=5000
      integer :: iord(5000)
      real*8 :: alist(5000),blist(5000),elist(5000),rlist(5000)
      real*8 abserr,result
      real*8 :: epsabs=1.0d-20
      real*8 :: epsrel=1.0d-3
c-----function definition starts here
c-----this is obnoxious, but assigning bmu here instead of as its own function dramatically speeds up this computation
      bmu=b0*(1+(r/rc)**2)**(-1.5*beta*eta)
      call dqagse(psyn,me,mx,
     &     epsabs,epsrel,limit,
     &     result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      
      jsyn=result
      return
      end function jsyn
c-----------------------------------------------------------------------
c-----integrand for los integral (DM density^2 times emissivity)
      real*8 function ssynintegrand(r)
      real*8 r,dist_z
      dist_z=dist(z)/(1+z)      !angular diameter distance
      ssynintegrand=4*pi/dist_z**2*r**2*rhox(r)**2*jsyn(r) !approximation; ok for small theta=radius/distance
      return
      end function ssynintegrand
c-----------------------------------------------------------------------
c-----dark matter integral along line of sight
      real*8 function ssyn(r)
      real*8 r
c-----integration parameters
      integer ier,last,neval
      integer :: limit=5000
      integer :: iord(5000)
      real*8 :: alist(5000),blist(5000),elist(5000),rlist(5000)
      real*8 abserr,result
      real*8 :: epsabs=1.0d-20
      real*8 :: epsrel=1.0d-3
c-----function definition
      call dqagse(ssynintegrand,1.d-16,r,
     &     epsabs,epsrel,limit,
     &     result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      
      ssyn=result
      
      return
      end function ssyn
c-----------------------------------------------------------------------
c-----minimum flux function (modified from eq. 10 in Cassano+2009)
      real*8 function min_flux(r)
      real*8 r
      real*8 thetaH,dist_z
      dist_z=dist(z)/(1+z)!angular diameter distance
      thetaH=r/dist_z*180.d0/pi*3600
      min_flux=4*log(2.0)*frms*(thetaH/thetaB)**2
      return
      end function min_flux
c-----------------------------------------------------------------------
c-----equivalent radius, where IC losses balance synch losses
c-----currently, this function is not used; rconst(ztemp) is used instead
      real*8 function req(ztemp)
      real*8 ztemp,beq,dist_z,r_b,req_temp
      dist_z=dist(z)/(1+z)      !angular diameter distance
      beq=sqrt(bic0/bsyn0)*(1+ztemp)**2 !magnetic field in muG where IC losses==synch losses
      r_b=0.5*dist_z*thetaB/(3600)*(pi/180) !beam size in cm
      if (beq<b0) then
         req_temp=rc*sqrt((b0/(sqrt(bic0/bsyn0)*(1+ztemp)**2))**(1/(1.5*beta*eta))-1)
         if (r_b<req_temp) then
            req=req_temp
         elseif (r_b>=req_temp) then
            req=r_b
         endif
      elseif (beq>=b0) then
         req_temp=rc*sqrt((1.01)**(1/(1.5*beta*eta))-1)
         if (r_b<req_temp) then
            req=req_temp
         elseif (r_b>=req_temp) then
            req=r_b
         endif
      endif
      return
      end function req
c-----------------------------------------------------------------------
c-----constant physical radius
c-----this function is currently used in the calculations
      real*8 function rconst(ztemp)
      real*8 ztemp,dist_z,r_b,r_temp
      dist_z=dist(z)/(1+z)      !angular diameter distance
      r_temp=rkpc*0.001*Mpccm   !kpc to Mpc to cm
      r_b=0.5*dist_z*thetaB/(3600)*(pi/180) !beam size in cm
      if (r_b<r_temp) then
         rconst=r_temp
      elseif (r_b>=r_temp) then
         rconst=r_b
      endif
      return
      end function rconst
c-----------------------------------------------------------------------
      end module synchrotron
c-----------------------------------------------------------------------