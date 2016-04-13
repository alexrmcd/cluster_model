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