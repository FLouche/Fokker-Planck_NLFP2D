!***************************************************
!
! Computation of several important physical quantities
!
!**************************************************
!
subroutine consts

use shared_plasma
use shared_grid
use shared_beam

use shared_timer
use shared_rf

use coulomb_log_mod

implicit none


double precision pi,twopi,gamma0,cinit,eonm
double precision :: teff

double precision rjtj,rr,tauel,rj
double precision cte0
double precision :: z1=0.d0
double precision arg,dennod
double PRECISION :: degtorad
double precision :: vcr, ecr
double precision :: lnae,lnab

integer i,j,ib

common/mathcons/pi,twopi

data eonm/9.57847729d7/
!data gamma0/3.5859906d0/ ! This is Lcoul*e^4/(4pi Eps0^2 mp^2) 
data gamma0/2.390775d-1/ ! This is e^4/(4pi Eps0^2 mp^2) 

!===================================================================
!
! Preliminary Computations
! ------------------------
!   -> Ions Densities 
!      --------------

nb(1)=ne
npart=ne*xpart
dennod=ne
do ib=2,nbulk
	nb(ib)=ne*xb(ib-1)
	dennod=dennod-zb(ib-1)*nb(ib)
enddo
dennod=dennod/za
!write(*,*) dennod,npart

!   -> Thermal velocities, gamma's, 
!             relatives masses (cfr. NRL )
!      -----------------------------------

do ib=1,nbulk

	if(ib == 1) then

		vt(ib)=4.19d5*dsqrt(t(ib))
        
        call coulomb_log_ae(ne, za, t(1), lnae)

		cte0=gamma0*lnae*(za/aa)**2

!  ---> Gamma for the electrons
!       -----------------------------

		gammab(ib)=cte0*nb(ib)
		maonmb(ib)=aa*1.8362d3
				      
		        else

		vt(ib)=9.79d3*dsqrt(t(ib)/ab(ib-1))
        
        call coulomb_log_ab(za, aa, t(2), npart, zb(ib-1), ab(ib-1), t(ib), nb(ib), lnab)
        
        cte0=gamma0*lnab*(za/aa)**2 ! strictly speaking this formula is only valid for thermal ions
        ! we use it for the steady-state solution

		gammab(ib)=cte0*nb(ib)*zb(ib-1)**2
		maonmb(ib)=aa/ab(ib-1)
		z1=z1+xb(ib-1)*zb(ib-1)**2*maonmb(ib)
        
        gammaa = cte0*npart*za**2

	endif
		  
enddo

vcr=9.d-2*(z1/aa)**.33333333d0*vt(1)

!  Vcr is the critical velocity where an equal amount of energy is transferred
!   from heated ions to background ions and electrons 
!     cfr Gaffey - JPP (1976) 16(2), pp. 146-169
!
!===================================================================
!
! Characteristic Relaxation Times (see Wesson 1997 pp 68-69)
! -------------------------------
!
!     -> Electron Collision Time
!        -----------------------
!
tauel=1.09d16*(t(1)/1000.d0)**1.5d0/npart/za**2/16.d0

!     -> Heat exchange time between ions and electrons
!        ---------------------------------------------
tauie=maonmb(1)/2.d0*tauel
!      write(6,*) ' '
!      write(6,*) 'Heat Exchange Time between Ions and Electrons = '
!    +       ,tauie,'  [s]'
spit=6.27d8*aa*t(1)**1.5/za**2/(ne*1.d-6)/15.d0

write(77,77) spit
77 format('Spitzer Slowing-down Time for resonating ion species = ', &
          D12.5,' [s]')

!===================================================================
!
! Initial distribution function parameters
! ----------------------------------------
!
   
allocate(fstix(nperp,npar))
   
!     -> Effective temperature  (Stix NF 15 737)
!        ---------------------
!
rj=0.d0
rjtj=0.d0
do j=1,nbulk-1
	rr=xb(j)*zb(j)**2*vt(1)/vt(j+1)
	rj=rj+rr
	rjtj=rjtj+rr*t(1)/t(j+1)
enddo

teff=t(1)*(1.d0+rj)/(1.d0+rjtj)

!     -> Effective thermal velocity  
!        --------------------------

vteff=9.79d3*dsqrt(teff/aa)
write(77,80) vteff
80     format('Effective thermal velocity of heated ions',D12.5)
write(77,81) teff*1.d-3
81 format('Effective Stix temperature of heated ions',D12.5,' keV')

!     -> Constant appearing in front of the exponential
!        ----------------------------------------------

cinit=npart/(2.d0*pi*vteff**2)**1.5d0

open(40, file=TRIM(outfile('fstix.dat')),status='unknown')	
do i=1,nperp
    do j=1,npar
        arg=(vpar(j)**2+vperp(i)**2)/(2.d0*vteff**2)
        fstix(i,j)=cinit*dexp(-arg)
        write(40,*) vperp(i),vpar(j),fstix(i,j)
    enddo
enddo
close(40)
!
!===================================================================

! Beam data's

degtorad = pi/180.d0

if (isource == -1) then
    beam_angle = beam_angle_deg*degtorad
    
    taum = 1.d0/taus
    
else
    
    taum = 0.d0
    
endif

!===================================================================

! RF Field data's

if (irf == -1) then

! --> RF fields rotating components

ceplus=dconjg(eplus)
cemin=dconjg(emin)

! RF pulsation 

omega=2.d0*pi*frek

! Phase velocity

vph = omega/kpar

! Cyclotron frequency of heated species

omc = 9.58d7*za/aa*b0

eta = kperp/omc

! Parallel resonant velocity

vres=(omega-nharm*omc)/kpar
write(*,*) 'Cyclotron frequency: ',omc/2.d0/pi/1d6, 'MHz'
write(*,*) 'Resonant velocity: ',vres, 'm/s'

! Constant appearing in front of the RF terms

rfcte = 0.125d0/kpar*pi*(za*eonm/aa)**2/vph**2

endif

vcr=9.d-2*(z1/aa)**.33333333d0*vt(1)
ecr=aa/2.d0/(9.79d0)**2*vcr**2*1.d-9

!  Vcr is the critical velocity where an equal amount of energy is transferred
!   from heated ions to background ions and electrons 
!     cfr Gaffey - JPP (1976) 16(2), pp. 146-169


write(*,*) 'Critical velocity: ', vcr, ' m/s'


end subroutine consts
