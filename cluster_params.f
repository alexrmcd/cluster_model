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