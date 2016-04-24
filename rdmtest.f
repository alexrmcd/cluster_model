c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----RadioDMcalc
c-----27 March 2015
c-----Calculates dark matter annihilation cross section, given some dark matter mass and annihilation channel, some total observed radio emission at some frequency
c-----To run:
c-----make RadioDMcalc
c-----./RadioDMcalc <channel number, 2 char> <redshift, 6 char>      
c-----example ./RadioDMcalc 25 0.0230
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----This line is 72 columns long--------------------------------------
c-----------------------------------------------------------------------
c-----Modules: declarations+definitions of global variables
c-----------------------------------------------------------------------
      module constants
      implicit none
      real*8 :: pi=3.14159265359
      real*8 :: Mpccm=3.085678d24 !Mpc to cm
      real*8 :: GeVJy=1.6d20    !GeV/s/Hz/cm^2 to Jy
      real*8 :: me=0.511d-3     !electron mass in GeV
c-----define constants for relationship between redshift and distance:
      real*8 :: clight=299792.458d0 !km/s
      real*8 :: h0=70.4d0       !km/s/Mpc
      real*8 :: omegaM=0.272d0
      real*8 :: omegaL=0.728d0
c-----------------------------------------------------------------------
      end module constants
c-----------------------------------------------------------------------
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
c-----------------------------------------------------------------------
      module cluster_params
      implicit none
c-----thermal electron density
      real*8 :: nele=1.d-3      !cm^-3
c-----magnetic field parameters
      real*8 :: beta=0.75       !beta from beta-model-fit 
      real*8 :: eta=0.5         !mag field scaling
c-----Coma; from Bonafede+2010
      real*8 :: b0=4.7        !central mag field value in microGauss
      real*8 :: rc=8.9793d23  !core radius in cm; 291kpc
c-----Perseus; from Chen+2007 and Taylor+2006
c      real*8 :: b0=25.0         !central mag field value in microGauss
c      real*8 :: rc=1.404d23     !core radius in cm; 46kpc
c-----dark matter NFW profile parameters; Coma specific; derived in mathematica script
      real*8 :: rhos=0.0399744  !dark matter characteristic density in GeV/cm^3
      real*8 :: rs=1.245518d24  !dark matter scale radius in cm
c-----------------------------------------------------------------------
      end module cluster_params
c-----------------------------------------------------------------------
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
      bloss= (bic0*(1+z)**4+bmu**2*bsyn0)*ener**2
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
      fff= 1.25d0*x**(1.d0/3.d0)*dexp(-x)*(648.d0+x**2)**(1.d0/12.d0)
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
      
      psyn=result*2*dndeeq(ener) !factor of 2 for positrons + electrons  tempremoved
      
      contains
c-----synchrotron single electron power, integrand
      real*8 function psynintegrand(theta)
      real*8 theta,xxx,nu_em
      real*8 :: psyn0=1.46323d-25 !GeV/s/Hz
      real*8 :: x0=62.1881      !dimensionless constant
      nu_em=(1+z)*nu            !emitted frequency = (1+z)(observed frequency)
      xxx=x0*nu_em/(bmu*ener**2)
      psynintegrand= psyn0*bmu*0.5*(dsin(theta))**2
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
c-----------------------------------------------------------------------
c-----Main Program
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      program radioDMcalc
      
      use constants
      use cluster_params
      use astro
      use synchrotron
      use particle

      implicit none

c-----Other variables, used only outside of functions
      real*8 sv,rmax,rmax_kpc,norm
c-----iterators, conversion factors, constants
      integer i,n
      real*8 j,ul,ll
c-----command line input variable names
      character(len=2) ch_char
      character(len=5) z_char
      character(len=4) mx_char
      character(len=3) r_char
      character(len=40) fname
c-----call darksusy initialization routine      
      call dsinit

c-----------------------------------------------------------------------
c---reference: dark matter annihilation channel numbers
c-----------------------------------------------------------------------
*** Ch No  Particles 
*** -----  ---------  
*** 13     W+ W-    
*** 15     e+ e-    
*** 17     mu+ mu-  
*** 19     tau+ tau-
*** 25     b b-bar  
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----set channel,redshift as command-line arguments
      call getarg(1,ch_char)
      read(ch_char,*) ch
      call getarg(2,z_char)
      read(z_char,*) z
c     call getarg(3,mx_char)
c     read(mx_char,*) mx
c      fname='mx_sv_r300kpc_ch'//ch_char//'_z'//z_char//'.txt'
      fname='test'//ch_char//'_z'//z_char//'.txt'
c     fname='nu_Ssyn_mx'//mx_char//'GeV_ch'
c     &     //ch_char//'_z'//z_char//'.txt'
c     fname='r_rhox_jsyn_ssyn_mx'//mx_char//'_ch'
c     &     //ch_char//'_z'//z_char//'.txt'
c-----open file for saving
      open(unit=21,file=fname)
c-----------------------------------------------------------------------
c-----THIS IS WHERE THE ACTION IS      
c-----------------------------------------------------------------------
c-----whatever you're actually trying to calculate goes here
c-----spectra in this case
      n=100
      rkpc=415                  !ROI radius in kpc
      rmax=rconst(z)
c      write(21,*)'#nu=',nu,'GHz | frms=',frms,'Jy/beam | B0=',b0,'muG | rmax=',rmax_kpc,'kpc | UL= ',min_flux(rmax),'Jy'
c      write(21,*) 'mass (GeV) | < sigma v > (cm^3 s^-1)'
      do i=0,n
         ll=5.d0                !lower limit on mx
         ul=1000.d0             !upper limit on mx
         mx=ll*(dexp((dlog(ul)-dlog(ll))/dble(n)*dble(i)))
         Sin=ssyn(rmax)*GeVJy   !ssyn also implicitly depends on DM mass, channel, redshift, magnetic field, etc
         Sout=min_flux(rmax) !min_flux also implicitly depends on beam size, redshift
         sv=8*pi*mx**2*(Sout/Sin) !units: cm^3 s^-1
         write(*,*)mx,sv
         write(21,103)mx,sv
      enddo
      
 103  format(20(1x,e20.10))     !formatting: ew.d, where w-d>7
      
      stop
      
      end program
c-----------------------------------------------------------------------
c-----ACTION ENDS HERE
c-----------------------------------------------------------------------      

