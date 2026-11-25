!==============================================================================
! cold_dispersion.f90 -- cold-plasma (Stix) dispersion relation for the RF
! parameters of an FP2D_QLRF_NL namelist.
!
! Usage:
!     cold_dispersion.exe  <FP2D namelist>  [<file with a &COLD group>]
!
! Reads from the &INPUT group of the FP2D namelist (all other entries are
! ignored, so any FP2D namelist works unchanged):
!     casename, ne, nbulk, xb(:), ab(:), zb(:)   electrons + bulk ions
!     xpart, aa, za                              minority (the FP2D species)
!     b0, frek, kpar, kperp, eplus, emin, nharm  RF
! Densities as in consts.f90: n_bulk(k) = ne*xb(k), n_minority = ne*xpart.
! Masses: A * m_p, as in the FP2D code (aa = 1 for H, 2 for D).
!
! Optional group &COLD (in the same file after &INPUT, or in a second file):
!     nu_omega  = 0       collision frequency / omega, added to every species
!                         (omega -> omega + i nu); regularises the resonances
!     nscan     = 0       0: point calculation at B = b0 only;
!                         >1: also a radial scan, B(R) = b0*r0/R, n uniform
!     r0        = 2.96    major radius where B = b0 (m)      [JET]
!     rmin      = 2.0     scan range (m)
!     rmax      = 4.0
!     nharm_max = 3       highest cyclotron harmonic listed in the layer table
!
! Physics. Wave fields ~ exp(i k.r - i omega t), B along z, k = (k_perp,0,k_par).
! Stix elements (sum over species s, Omega_s = Z_s e B / m_s signed):
!     R = 1 - sum wp^2 / (w (w~ + Omega)),   L = 1 - sum wp^2 / (w (w~ - Omega))
!     P = 1 - sum wp^2 / (w w~),   S = (R+L)/2,   D = (R-L)/2,   w~ = w + i nu
! Wave equation  n x (n x E) + K.E = 0, i.e. with n_par = c k_par / omega,
!     | S - n_par^2   -iD             n_par n_perp |
!     | iD            S - n^2         0            | E = 0
!     | n_par n_perp  0               P - n_perp^2 |
! whose determinant gives the biquadratic in x = n_perp^2:
!     S x^2 - [ (S - n_par^2)(S + P) - D^2 ] x + P (R - n_par^2)(L - n_par^2) = 0
! Two roots: fast wave (smaller |x|) and slow / ion Bernstein-like branch
! (larger |x|). Polarisation of each root (rows 2 and 3):
!     E_y/E_x = -i D / (S - n^2),   E_z/E_x = -n_par n_perp / (P - n_perp^2)
!     E+ = (E_x + i E_y)/sqrt2  (left-hand, rotates with the ions)
!     E- = (E_x - i E_y)/sqrt2  (right-hand)
!     E-/E+ = (L - n^2) / (R - n^2)
! E+ and E- are those of the FP2D QL operator (E+ multiplies J_{N-1}).
!
! At the minority layer a cold minority makes L singular, so the polarisation
! is also given for the dielectric WITHOUT the minority ("bulk only"), which is
! the usual estimate of E-/E+ seen by a trace minority.
!
! Output: stdout report; with nscan > 1 also cold_dispersion-<casename>.txt.
!==============================================================================
module cold_kinds
   implicit none
   integer, parameter :: dp = kind(1.d0)
   real(dp), parameter :: pi   = 3.141592653589793d0
   real(dp), parameter :: qe   = 1.602176634d-19
   real(dp), parameter :: me   = 9.1093837015d-31
   real(dp), parameter :: mp   = 1.67262192369d-27
   real(dp), parameter :: eps0 = 8.8541878128d-12
   real(dp), parameter :: clight = 2.99792458d8
end module cold_kinds

!------------------------------------------------------------------------------
! Extract one namelist group from a file, keeping only selected variables, so
! that a namelist with many other entries can be read with a short NAMELIST.
!------------------------------------------------------------------------------
module nml_extract
   use cold_kinds
   implicit none
contains

   ! Whole file as one string, '!' comments removed (outside quotes).
   subroutine read_text(fname, text, ok)
      character(*), intent(in) :: fname
      character(:), allocatable, intent(out) :: text
      logical, intent(out) :: ok
      integer :: u, ios, i
      character(len=4096) :: line
      character :: q
      text = ''
      open(newunit=u, file=fname, status='old', action='read', iostat=ios)
      ok = (ios == 0)
      if (.not. ok) return
      do
         read(u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         q = ' '
         do i = 1, len_trim(line)
            if (q /= ' ') then
               if (line(i:i) == q) q = ' '
            else if (line(i:i) == '''' .or. line(i:i) == '"') then
               q = line(i:i)
            else if (line(i:i) == '!') then
               line(i:) = ' '
               exit
            end if
         end do
         text = text // ' ' // trim(line)
      end do
      close(u)
   end subroutine read_text

   pure function lower(s) result(r)
      character(*), intent(in) :: s
      character(len(s)) :: r
      integer :: i, k
      r = s
      do i = 1, len(s)
         k = iachar(s(i:i))
         if (k >= 65 .and. k <= 90) r(i:i) = achar(k + 32)
      end do
   end function lower

   pure logical function is_name_char(c)
      character, intent(in) :: c
      is_name_char = (c >= 'a' .and. c <= 'z') .or. (c >= 'A' .and. c <= 'Z') .or. &
                     (c >= '0' .and. c <= '9') .or. c == '_'
   end function is_name_char

   ! Build '&<group> name=value ... /' from <text>, keeping the variables whose
   ! (lower-case) names appear in the blank-separated list <wanted>
   ! (wanted = '*' keeps everything). found = .false. if the group is absent.
   subroutine extract_group(text, group, wanted, buf, found)
      character(*), intent(in) :: text, group, wanted
      character(:), allocatable, intent(out) :: buf
      logical, intent(out) :: found
      character(:), allocatable :: lt
      integer :: i0, i1, i, j, k, depth, na, nstart(400)
      character :: q
      character(len=64) :: nm

      lt = lower(text)
      buf = ''
      found = .false.
      ! locate '&group' followed by a non-name character
      i0 = 0
      i = 1
      do
         k = index(lt(i:), '&' // lower(trim(group)))
         if (k == 0) return
         j = i + k - 1 + len_trim(group) + 1          ! first char after the name
         if (j > len(lt)) then
            i0 = j; exit
         end if
         if (.not. is_name_char(lt(j:j))) then
            i0 = j; exit
         end if
         i = j
      end do
      found = .true.
      ! scan the body: end at '/' (or '&end') outside quotes and parentheses,
      ! record the start of every 'name =' / 'name(..) =' assignment
      q = ' '; depth = 0; na = 0; i1 = len(lt)
      do i = i0, len(lt)
         if (q /= ' ') then
            if (lt(i:i) == q) q = ' '
            cycle
         end if
         select case (lt(i:i))
         case ('''', '"'); q = lt(i:i)
         case ('(');       depth = depth + 1
         case (')');       depth = depth - 1
         case ('/')
            if (depth == 0) then
               i1 = i - 1; exit
            end if
         case ('&')
            if (depth == 0) then
               i1 = i - 1; exit
            end if
         case ('=')
            if (depth == 0) then
               j = i - 1
               do while (j > i0 .and. lt(j:j) == ' ')
                  j = j - 1
               end do
               if (lt(j:j) == ')') then              ! array element: name(...)=
                  do while (j > i0 .and. lt(j:j) /= '(')
                     j = j - 1
                  end do
                  j = j - 1
               end if
               do while (j > i0 .and. is_name_char(lt(j-1:j-1)))
                  j = j - 1
               end do
               na = na + 1
               nstart(na) = j
            end if
         end select
      end do
      buf = '&' // trim(group)
      do k = 1, na
         j = nstart(k)
         i = merge(nstart(k+1) - 1, i1, k < na)
         ! variable name (up to '(' , ' ' or '=')
         nm = ''
         do i0 = j, i
            if (.not. is_name_char(lt(i0:i0))) exit
         end do
         nm = lt(j:i0-1)
         if (trim(wanted) == '*' .or. index(' ' // wanted // ' ', ' ' // trim(nm) // ' ') > 0) then
            buf = buf // ' ' // trim(adjustl(text(j:i)))
            ! make sure each kept assignment ends with a separator
            if (buf(len(buf):len(buf)) /= ',') buf = buf // ','
         end if
      end do
      buf = buf // ' /'
   end subroutine extract_group

end module nml_extract

!------------------------------------------------------------------------------
module cold_plasma
   use cold_kinds
   implicit none

   type species_t
      character(len=12) :: name
      real(dp) :: z, a, n, m          ! charge number, mass number, density, mass
   end type species_t

contains

   ! Stix R, L, P (S, D derived) for omega, B. Species with use=.false. skipped.
   subroutine stix(sp, use, omega, b, nu_omega, r, l, p, s, d)
      type(species_t), intent(in) :: sp(:)
      logical, intent(in) :: use(:)
      real(dp), intent(in) :: omega, b, nu_omega
      complex(dp), intent(out) :: r, l, p, s, d
      complex(dp) :: wt
      real(dp) :: wp2, oc
      integer :: k
      wt = cmplx(omega, nu_omega*omega, dp)
      r = (1.d0, 0.d0); l = r; p = r
      do k = 1, size(sp)
         if (.not. use(k)) cycle
         wp2 = sp(k)%n * (sp(k)%z*qe)**2 / (eps0*sp(k)%m)
         oc  = sp(k)%z*qe*b / sp(k)%m
         r = r - wp2 / (omega*(wt + oc))
         l = l - wp2 / (omega*(wt - oc))
         p = p - wp2 / (omega*wt)
      end do
      s = 0.5d0*(r + l)
      d = 0.5d0*(r - l)
   end subroutine stix

   ! Roots x = n_perp^2 of  a x^2 + b x + c = 0 (cancellation-free form),
   ! returned as x(1) = smaller |x| (fast wave), x(2) = larger (slow branch).
   subroutine nperp2_roots(s, d, p, r, l, npar, x)
      complex(dp), intent(in) :: s, d, p, r, l
      real(dp), intent(in) :: npar
      complex(dp), intent(out) :: x(2)
      complex(dp) :: a, b, c, sq, q, t
      real(dp) :: n2
      n2 = npar**2
      a = s
      b = -((s - n2)*(s + p) - d*d)
      c = p*(r - n2)*(l - n2)
      sq = sqrt(b*b - 4.d0*a*c)
      if (real(conjg(b)*sq) < 0.d0) sq = -sq
      q = -0.5d0*(b + sq)
      if (abs(q) == 0.d0) then
         x = (0.d0, 0.d0); return
      end if
      x(1) = c/q
      if (abs(a) > 0.d0) then
         x(2) = q/a
      else
         x(2) = cmplx(huge(1.d0), 0.d0, dp)
      end if
      if (abs(x(1)) > abs(x(2))) then
         t = x(1); x(1) = x(2); x(2) = t
      end if
   end subroutine nperp2_roots

   ! Field ratios for one root; nperp = principal sqrt of x.
   subroutine polarisation(s, d, p, r, l, npar, x, emep, eyex, ezex, resid)
      complex(dp), intent(in) :: s, d, p, r, l, x
      real(dp), intent(in) :: npar
      complex(dp), intent(out) :: emep, eyex, ezex
      real(dp), intent(out) :: resid
      complex(dp) :: n2, np, m(3,3), det, sc
      n2 = x + npar**2
      np = root_branch(x)
      eyex = -(0.d0, 1.d0)*d/(s - n2)
      ezex = -npar*np/(p - x)
      emep = (l - n2)/(r - n2)
      ! check: determinant of the wave matrix, relative to the size of its terms
      m(1,:) = [s - npar**2, -(0.d0,1.d0)*d, npar*np]
      m(2,:) = [(0.d0,1.d0)*d, s - n2, (0.d0,0.d0)]
      m(3,:) = [npar*np, (0.d0,0.d0), p - x]
      det = m(1,1)*(m(2,2)*m(3,3) - m(2,3)*m(3,2)) - m(1,2)*(m(2,1)*m(3,3) - m(2,3)*m(3,1)) &
          + m(1,3)*(m(2,1)*m(3,2) - m(2,2)*m(3,1))
      sc = abs(m(1,1)*m(2,2)*m(3,3)) + abs(m(1,2)*m(2,1)*m(3,3)) + abs(m(1,3)*m(3,1)*m(2,2))
      resid = abs(det)/max(abs(sc), tiny(1.d0))
   end subroutine polarisation

   ! sqrt with the branch Re >= 0, and Im >= 0 when Re = 0 (evanescent: decays
   ! towards +x, the direction of propagation chosen for the propagating roots)
   pure complex(dp) function root_branch(x) result(z)
      complex(dp), intent(in) :: x
      z = sqrt(x)
      if (real(z) < 0.d0 .or. (real(z) == 0.d0 .and. aimag(z) < 0.d0)) z = -z
   end function root_branch

end module cold_plasma

!==============================================================================
program cold_dispersion
   use cold_kinds
   use nml_extract
   use cold_plasma
   implicit none

   ! --- &INPUT variables (subset of the FP2D namelist), FP2D defaults -----
   character(len=256) :: casename = ''
   integer  :: nbulk = 2, nharm = 1
   real(dp) :: ne = 0.d0, xpart = 0.d0, aa = 1.d0, za = 1.d0
   real(dp) :: ab(9) = 0.d0, zb(9) = 0.d0, xb(9) = 0.d0
   real(dp) :: b0 = 0.d0, frek = 0.d0, kpar = 0.d0
   complex(dp) :: kperp = (0.d0,0.d0), eplus = (0.d0,0.d0), emin = (0.d0,0.d0)
   namelist /input/ casename, nbulk, nharm, ne, xpart, aa, za, ab, zb, xb, &
                    b0, frek, kpar, kperp, eplus, emin
   ! --- &COLD --------------------------------------------------------------
   real(dp) :: nu_omega = 0.d0, r0 = 2.96d0, rmin = 2.0d0, rmax = 4.0d0
   integer  :: nscan = 0, nharm_max = 3
   namelist /cold/ nu_omega, r0, rmin, rmax, nscan, nharm_max

   character(len=512) :: fnml, fcold, fout
   character(:), allocatable :: text, buf
   logical :: ok, found
   integer :: ios, nsp, k, i, iu, nh
   type(species_t), allocatable :: sp(:)
   logical, allocatable :: all_sp(:), no_min(:)
   real(dp) :: omega, npar, b, rr, oc, dens_check
   complex(dp) :: r, l, p, s, d, x(2), xb2(2), rb, lb, pb, sb, db
   character(len=400) :: msg

   if (command_argument_count() < 1) then
      write(*,'(a)') 'usage: cold_dispersion.exe <FP2D namelist> [<file with &COLD>]'
      stop 1
   end if
   call get_command_argument(1, fnml)
   fcold = fnml
   if (command_argument_count() >= 2) call get_command_argument(2, fcold)

   ! --- read &INPUT ---------------------------------------------------------
   call read_text(fnml, text, ok)
   if (.not. ok) then
      write(*,'(2a)') 'cannot open ', trim(fnml); stop 1
   end if
   call extract_group(text, 'input', &
        'casename nbulk nharm ne xpart aa za ab zb xb b0 frek kpar kperp eplus emin', buf, found)
   if (.not. found) then
      write(*,'(2a)') 'no &INPUT group in ', trim(fnml); stop 1
   end if
   read(buf, nml=input, iostat=ios, iomsg=msg)
   if (ios /= 0) then
      write(*,'(a)') 'error reading &INPUT: ' // trim(msg)
      write(*,'(a)') 'extracted text: ' // buf
      stop 1
   end if
   ! --- read &COLD (optional) -----------------------------------------------
   call read_text(fcold, text, ok)
   if (ok) then
      call extract_group(text, 'cold', '*', buf, found)
      if (found) then
         read(buf, nml=cold, iostat=ios, iomsg=msg)
         if (ios /= 0) then
            write(*,'(a)') 'error reading &COLD: ' // trim(msg); stop 1
         end if
      end if
   else
      write(*,'(2a)') 'cannot open ', trim(fcold); stop 1
   end if
   if (ne <= 0.d0 .or. b0 <= 0.d0 .or. frek <= 0.d0) then
      write(*,'(a)') 'ne, b0 and frek must be set and positive'; stop 1
   end if
   if (nbulk < 1 .or. nbulk > 10) then
      write(*,'(a)') 'nbulk must be 1..10'; stop 1
   end if

   ! --- species: electrons, bulk ions (nbulk-1, as in consts.f90), minority -
   nsp = 1 + (nbulk - 1) + 1
   allocate(sp(nsp), all_sp(nsp), no_min(nsp))
   sp(1) = species_t('electrons', -1.d0, me/mp, ne, me)
   do k = 1, nbulk - 1
      write(sp(1+k)%name, '(a,i0)') 'bulk ion ', k
      sp(1+k)%z = zb(k); sp(1+k)%a = ab(k); sp(1+k)%n = ne*xb(k); sp(1+k)%m = ab(k)*mp
   end do
   sp(nsp) = species_t('minority', za, aa, ne*xpart, aa*mp)
   all_sp = .true.
   no_min = .true.; no_min(nsp) = .false.
   do k = 2, nsp
      if (sp(k)%a <= 0.d0) then
         write(*,'(a,i0)') 'zero mass for species ', k; stop 1
      end if
   end do

   omega = 2.d0*pi*frek
   npar  = clight*kpar/omega

   ! --- report: plasma ------------------------------------------------------
   write(*,'(a)') repeat('=', 78)
   write(*,'(a)') ' Cold-plasma dispersion   case: ' // trim(casename)
   write(*,'(a)') '   namelist: ' // trim(fnml)
   write(*,'(a)') repeat('=', 78)
   write(*,'(a,f8.4,a,f9.4,a,es11.4,a)') ' B0 = ', b0, ' T   f = ', frek/1.d6, ' MHz   omega = ', omega, ' rad/s'
   write(*,'(a,f8.3,a,f9.4,a,es10.3)') ' k_par = ', kpar, ' 1/m   n_par = ', npar, '   nu/omega = ', nu_omega
   write(*,'(/a)') ' species        Z      A       n (m^-3)     f_p (MHz)     f_c (MHz)   omega/Omega_c'
   dens_check = 0.d0
   do k = 1, nsp
      oc = abs(sp(k)%z)*qe*b0/sp(k)%m
      write(*,'(1x,a12,f5.1,f8.4,es14.4,2es14.4,f12.4)') sp(k)%name, sp(k)%z, sp(k)%a, sp(k)%n, &
         sqrt(sp(k)%n*(sp(k)%z*qe)**2/(eps0*sp(k)%m))/(2.d0*pi)/1.d6, oc/(2.d0*pi)/1.d6, omega/oc
      if (k > 1) dens_check = dens_check + sp(k)%z*sp(k)%n
   end do
   if (abs(dens_check - ne) > 1.d-6*ne) &
      write(*,'(a,f8.5)') ' WARNING: quasi-neutrality violated, sum(Z n_i)/n_e = ', dens_check/ne
   write(*,'(a,f7.4,a,i0)') ' minority: omega / (N Omega_c) = ', &
      omega/(nharm*abs(za)*qe*b0/(aa*mp)), '  for nharm = ', nharm

   ! --- point calculation at B = b0 -----------------------------------------
   call stix(sp, all_sp, omega, b0, nu_omega, r, l, p, s, d)
   call stix(sp, no_min, omega, b0, nu_omega, rb, lb, pb, sb, db)
   write(*,'(/a)') ' Stix elements at B0              all species                      bulk only (no minority)'
   call wc(' R ', r, rb); call wc(' L ', l, lb); call wc(' S ', s, sb)
   call wc(' D ', d, db); call wc(' P ', p, pb)

   call nperp2_roots(s, d, p, r, l, npar, x)
   call nperp2_roots(sb, db, pb, rb, lb, npar, xb2)
   write(*,'(/a)') ' Roots n_perp^2 (k_perp along x), all species:'
   call report_root('fast wave  ', s, d, p, r, l, x(1))
   call report_root('slow branch', s, d, p, r, l, x(2))
   write(*,'(/a)') ' Roots n_perp^2, bulk only (minority left out of the dielectric):'
   call report_root('fast wave  ', sb, db, pb, rb, lb, xb2(1))
   call report_root('slow branch', sb, db, pb, rb, lb, xb2(2))
   write(*,'(/a,2es12.4)') ' FW estimate (R-n||^2)(L-n||^2)/(S-n||^2), all species: ', &
      (r - npar**2)*(l - npar**2)/(s - npar**2)

   ! --- comparison with the FP2D RF input -----------------------------------
   write(*,'(/a)') ' Compared with the FP2D RF input:'
   write(*,'(a,f9.3,a)') '   namelist kperp            = ', real(kperp), ' 1/m'
   write(*,'(a,f9.3,a,f9.3,a)') '   fast-wave Re k_perp       = ', real(omega/clight*root_branch(x(1))), &
      ' 1/m (all species),', real(omega/clight*root_branch(xb2(1))), ' 1/m (bulk only)'
   if (abs(eplus) > 0.d0) then
      write(*,'(a,f9.3)') '   namelist |emin/eplus|     = ', abs(emin/eplus)
   end if
   call polar_only(s, d, p, r, l, x(1), msg)
   write(*,'(a)') '   fast-wave |E-/E+|         = ' // trim(msg) // ' (all species)'
   call polar_only(sb, db, pb, rb, lb, xb2(1), msg)
   write(*,'(a)') '                               ' // trim(msg) // ' (bulk only)'

   ! --- cyclotron layers ----------------------------------------------------
   write(*,'(/a,f6.3,a)') ' Ion cyclotron layers, B = B0 R0/R with R0 = ', r0, ' m:'
   do k = 2, nsp
      do nh = 1, nharm_max
         oc = abs(sp(k)%z)*qe*b0/sp(k)%m
         write(*,'(3x,a12,a,i1,a,f7.3,a,f7.4,a)') sp(k)%name, ' N=', nh, '  R = ', &
            r0*nh*oc/omega, ' m   (B = ', omega/(nh*abs(sp(k)%z)*qe/sp(k)%m), ' T)'
      end do
   end do

   ! --- radial scan ---------------------------------------------------------
   if (nscan > 1) then
      if (len_trim(casename) > 0) then
         fout = 'cold_dispersion-' // trim(casename) // '.txt'
      else
         fout = 'cold_dispersion.txt'
      end if
      open(newunit=iu, file=fout, status='replace', action='write')
      write(iu,'(a)') '# cold-plasma dispersion, radial scan B = b0*r0/R, uniform density; case ' // trim(casename)
      write(iu,'(a,f8.4,a,f8.4,a,f9.4,a,f8.3,a,es10.3)') '# b0 = ', b0, ' T  r0 = ', r0, ' m  f = ', &
         frek/1d6, ' MHz  kpar = ', kpar, ' 1/m  nu/omega = ', nu_omega
      write(iu,'(a,f9.3,a,f9.3,a,i0,a,f8.4)') '# namelist kperp = ', real(kperp), '  |emin/eplus| = ', &
         merge(abs(emin/eplus), 0.d0, abs(eplus) > 0.d0), '  nharm = ', nharm, '  minority A/Z = ', aa/za
      write(iu,'(a)') '# columns: 1 R(m) 2 B(T) 3 omega/Omega_min 4-5 Re,Im nperp2 fast 6-7 Re,Im nperp2 slow' // &
         ' 8-9 Re,Im kperp fast(1/m) 10 |E-/E+| fast 11 Re nperp2 fast (bulk only) 12 |E-/E+| fast (bulk only)' // &
         ' 13 Re(S-npar2) 14 Re(L-npar2) 15 Re(R-npar2) 16 Re P'
      do i = 1, nscan
         rr = rmin + (rmax - rmin)*(i - 1)/real(nscan - 1, dp)
         b  = b0*r0/rr
         call stix(sp, all_sp, omega, b, nu_omega, r, l, p, s, d)
         call nperp2_roots(s, d, p, r, l, npar, x)
         call stix(sp, no_min, omega, b, nu_omega, rb, lb, pb, sb, db)
         call nperp2_roots(sb, db, pb, rb, lb, npar, xb2)
         write(iu,'(f9.5,f9.5,f10.5,14es14.5e3)') rr, b, omega/(abs(za)*qe*b/(aa*mp)), &
            x(1), x(2), omega/clight*root_branch(x(1)), abs((l - x(1) - npar**2)/(r - x(1) - npar**2)), &
            real(xb2(1)), abs((lb - xb2(1) - npar**2)/(rb - xb2(1) - npar**2)), &
            real(s - npar**2), real(l - npar**2), real(r - npar**2), real(p)
      end do
      close(iu)
      call find_layers(fout)
      write(*,'(/a,i0,a)') ' Radial scan (', nscan, ' points) written to ' // trim(fout)
   end if

contains

   subroutine wc(lab, z1, z2)
      character(*), intent(in) :: lab
      complex(dp), intent(in) :: z1, z2
      write(*,'(1x,a,2(3x,"(",es13.5,",",es12.4,")"))') lab, z1, z2
   end subroutine wc

   subroutine report_root(lab, s_, d_, p_, r_, l_, x_)
      character(*), intent(in) :: lab
      complex(dp), intent(in) :: s_, d_, p_, r_, l_, x_
      complex(dp) :: emep, eyex, ezex, kp
      real(dp) :: resid
      character(len=12) :: kind_
      call polarisation(s_, d_, p_, r_, l_, npar, x_, emep, eyex, ezex, resid)
      kp = omega/clight*root_branch(x_)
      kind_ = merge('propagating', 'evanescent ', real(x_) > 0.d0)
      write(*,'(3x,a,"  n_perp^2 = (",es12.4,",",es11.3,")  ",a)') lab, x_, trim(kind_)
      if (abs(real(kp)) >= abs(aimag(kp))) then
         write(*,'(5x,"k_perp = (",es11.4,",",es11.3,") 1/m   lambda_perp = ",es10.3," m")') &
            kp, 2.d0*pi/abs(real(kp))
      else
         write(*,'(5x,"k_perp = (",es11.4,",",es11.3,") 1/m   decay length 1/Im k = ",es10.3," m")') &
            kp, 1.d0/abs(aimag(kp))
      end if
      write(*,'(5x,"E-/E+ = (",es11.3,",",es11.3,")  |E-/E+| = ",es10.3)') emep, abs(emep)
      write(*,'(5x,"Ey/Ex = (",es11.3,",",es11.3,")   Ez/Ex = (",es11.3,",",es11.3,")")') eyex, ezex
      write(*,'(5x,"check |det|/scale = ",es9.2)') resid
   end subroutine report_root

   subroutine polar_only(s_, d_, p_, r_, l_, x_, txt)
      complex(dp), intent(in) :: s_, d_, p_, r_, l_, x_
      character(*), intent(out) :: txt
      complex(dp) :: emep, eyex, ezex
      real(dp) :: resid
      call polarisation(s_, d_, p_, r_, l_, npar, x_, emep, eyex, ezex, resid)
      write(txt,'(f9.3)') abs(emep)
   end subroutine polar_only

   ! Zero crossings of S - n_par^2 (ion-ion hybrid / FW resonance), L - n_par^2
   ! and R - n_par^2 (FW cutoffs), P (should be none in ICRF), and sign
   ! changes of Re n_perp^2 of the fast wave, found from the scan table.
   subroutine find_layers(fn)
      character(*), intent(in) :: fn
      integer :: u, ios, n, j, c, jm, jp
      character(len=4) :: kind_
      real(dp), allocatable :: t(:,:)
      real(dp) :: row(16)
      character(len=512) :: line
      character(len=40), parameter :: lab(5) = [character(40) :: &
         'S = n_par^2  (ion-ion hybrid resonance)', 'L = n_par^2  (L cutoff)', &
         'R = n_par^2  (R cutoff)', 'P = 0', 'Re n_perp^2 (fast) changes sign']
      integer, parameter :: col(5) = [13, 14, 15, 16, 4]
      open(newunit=u, file=fn, status='old', action='read')
      n = 0
      do
         read(u,'(a)', iostat=ios) line
         if (ios /= 0) exit
         if (line(1:1) /= '#') n = n + 1
      end do
      allocate(t(16, n))
      rewind(u)
      j = 0
      do
         read(u,'(a)', iostat=ios) line
         if (ios /= 0) exit
         if (line(1:1) == '#') cycle
         read(line, *) row
         j = j + 1; t(:, j) = row
      end do
      close(u)
      write(*,'(/a)') ' Sign changes in the scan (linear interpolation; zero = layer, pole = cold resonance):'
      do c = 1, 5
         do j = 1, n - 1
            if (t(col(c), j)*t(col(c), j+1) < 0.d0) then
               ! a zero: |value| smaller at the crossing than one step outside;
               ! a pole (cold cyclotron resonance): |value| grows towards it
               jm = max(j - 1, 1); jp = min(j + 2, n)
               if (min(abs(t(col(c),j)), abs(t(col(c),j+1))) > &
                   min(abs(t(col(c),jm)), abs(t(col(c),jp)))) then
                  kind_ = 'pole'
               else
                  kind_ = 'zero'
               end if
               write(*,'(3x,a40,1x,a4,"  R = ",f7.3," m  (B = ",f7.4," T)")') lab(c), kind_, &
                  t(1,j) - t(col(c),j)*(t(1,j+1) - t(1,j))/(t(col(c),j+1) - t(col(c),j)), &
                  t(2,j) - t(col(c),j)*(t(2,j+1) - t(2,j))/(t(col(c),j+1) - t(col(c),j))
            end if
         end do
      end do
   end subroutine find_layers

end program cold_dispersion
