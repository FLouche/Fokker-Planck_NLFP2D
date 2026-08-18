!***************************************************
!
! Computation of relevant physical quantities
!
!  26/05/2026: updated with consistent evaluation of the 
!              Coulomb logarithm with NRL formulas (F. Louche)
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
use mod_ncint

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
! characteristic collision times, evaluated after vteff below; named distinctly
! from tauel/tauie above, which are the Wesson relaxation times
double precision :: tcoll_ii, tcoll_ie, tcoll_ei
double precision :: sum_nz2, lnaa

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

		! The linearised operator annihilates exp(-v^2/2 vth_a^2) only when
		! vt(b)^2 / maonmb(b) == vth_a^2 = (9.79d3)^2 * T/aa.  On the ion
		! branch below the same 9.79d3 appears on both sides and cancels
		! identically.  On the electron branch it did not: the literal 4.19d5
		! implies (4.19d5/9.79d3)^2 = 1831.7 for the mass ratio, while
		! maonmb(1) uses 1836.2 -- a 0.243% inconsistency, which left the
		! electron operator's equilibrium a Maxwellian at 0.24%-shifted
		! temperature instead of Te (measured: continuous residual L*M/M ~ 2e-3
		! for electrons, 1e-16 for deuterons).
		!
		! Deriving sqrt(e/me) from the ion constant and the mass ratio makes
		! the identity exact by construction.  Numerically this is 4.19503d5
		! against the previous 4.19d5, a +0.12% change to the electron thermal
		! velocity only; every other use of 9.79d3 in the code is untouched.
		vt(ib)=9.79d3*dsqrt(1.8362d3)*dsqrt(t(ib))
        
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
        
!        gammaa = cte0*npart*za**2

	endif
		  
enddo

vcr=9.d-2*(z1/aa)**.33333333d0*vt(1)

!  Vcr is the critical velocity where an equal amount of energy is transferred
!   from heated ions to background ions and electrons 
!     cfr Gaffey - JPP (1976) 16(2), pp. 146-169

write(*,66) vcr
66 format('Critical velocity = ',D12.5,' [m/s]')
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

write(*,77) spit
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
write(*,80) vteff
80     format('Effective thermal velocity of heated ions',D12.5)
write(*,81) teff*1.d-3
81 format('Effective Stix temperature of heated ions',D12.5,' keV')

!===================================================================
!
! Characteristic Collision Times
! ------------------------------
!
!  Placed here rather than in the relaxation-time section above because they
!  use vteff, the thermal velocity at the Stix effective temperature at which
!  the resonant species is initialised, which is only known at this point.
!
!  These are the standard Wesson (1997, p. 69) expressions evaluated at the
!  initial state, with the code's own NRL Coulomb logarithms rather than the
!  fixed values used by tauel/spit above.  They are reference timescales, not
!  the instantaneous rates the operator applies: the Fokker-Planck coefficients
!  scale as gammab/v**3, but that form is the fast-test-particle limit and is
!  not valid for ions colliding on the much faster electrons, so it must not be
!  used to build an ion-electron time.
!
!  "Ion" here is the resonant species, which is the one the code evolves.
!  tauel and tauie above are deliberately left untouched: tauie sets the
!  reference rate of the steady-state convergence diagnostic.

!     -> Electron-ion: electron collision time
!        -------------------------------------
!  tau_e = 1.09d16 * Te[keV]**1.5 / ( sum_i n_i Z_i**2 * lnae ), the sum running
!  over every ion species present, which is what the electrons collide against.
sum_nz2=npart*za**2
do ib=2,nbulk
    sum_nz2=sum_nz2+nb(ib)*zb(ib-1)**2
enddo
if (sum_nz2 > 0.d0 .and. lnae > 0.d0) then
    tcoll_ei=1.09d16*(t(1)*1.d-3)**1.5d0/(sum_nz2*lnae)
else
    tcoll_ei=0.d0
endif

!     -> Ion-ion: like-particle collisions of the resonant species
!        ----------------------------------------------------------
!  tau_i = 6.60d17 * sqrt(ma/mp) * Ti[keV]**1.5 / ( na * Za**4 * lnaa ),
!  evaluated at the Stix effective temperature at which it is initialised.
call coulomb_log_ab(za, aa, teff, npart, za, aa, teff, npart, lnaa)
if (npart > 0.d0 .and. lnaa > 0.d0) then
    tcoll_ii=6.60d17*dsqrt(aa)*(teff*1.d-3)**1.5d0/(npart*za**4*lnaa)
else
    tcoll_ii=0.d0
endif

!     -> Ion-electron: energy equipartition time
!        ----------------------------------------
!  tau_ie = (ma / 2 me) * tau_e
tcoll_ie=maonmb(1)/2.d0*tcoll_ei

write(*,82) tcoll_ii
82 format('Collision time ion-ion      (resonant, like-particle) = ', &
          D12.5,' [s]')
write(*,83) tcoll_ie
83 format('Collision time ion-electron (equipartition)           = ', &
          D12.5,' [s]')
write(*,84) tcoll_ei
84 format('Collision time electron-ion (electron collision time) = ', &
          D12.5,' [s]')

!     -> Constant appearing in front of the exponential
!        ----------------------------------------------

cinit=npart/(2.d0*pi*vteff**2)**1.5d0

do i=1,nperp
    do j=1,npar
        arg=(vpar(j)**2+vperp(i)**2)/(2.d0*vteff**2)
        fstix(i,j)=cinit*dexp(-arg)
    enddo
enddo

! Renormalise fstix to exactly npart on the numerical grid.
! The analytical cinit gives the correct continuous integral, but the
! discretised quadrature may differ by a small amount.  Rescaling here ensures
! time_density returns exactly npart at t=0 for any grid configuration.
block
    double precision :: fstix_dens
    double precision :: fint_tmp(nperp,npar)
    integer :: ii, jj
    do ii = 1, nperp
        do jj = 1, npar
            fint_tmp(ii,jj) = fstix(ii,jj) * jacob(ii,jj)
        enddo
    enddo
    call ncint_2d(fint_tmp, fstix_dens)
    if (fstix_dens > 0.d0) fstix = fstix * npart / fstix_dens
end block

! fstix.dat is written AFTER the renormalisation, so the file holds the same
! array the solver uses.  It used to be written before, which meant the file
! carried the continuum-normalised Maxwellian while every computed solution was
! grid-normalised through fout = fout*npart/dens_tmp.  Any comparison against
! the file therefore inherited the quadrature's own error as a spurious
! amplitude offset -- measured at 8.2e-5 on a 241x241 uniform grid, and
! converging only at the order of the quadrature rather than of the scheme.
open(40, file=TRIM(outfile('fstix.dat')),status='unknown')
do i=1,nperp
    do j=1,npar
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
