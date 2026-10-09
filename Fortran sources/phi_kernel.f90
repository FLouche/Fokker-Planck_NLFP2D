!***********************************************************************
!*  mod_phi_kernel                                                     *
!*                                                                     *
!*  Structured replacement for the dense phi-distance matrix sum_phi.  *
!*                                                                     *
!*  BACKGROUND                                                         *
!*  ----------                                                         *
!*  grid.f90 used to expand the kernel into a dense nbig x nbig matrix  *
!*                                                                     *
!*      sum_phi(ix1, ix2) = kern(|j-jp|, i, ip),  ix = (i-1)*npar + j   *
!*                                                                     *
!*  and nlterm_test.f90 used it in exactly one way: a single DGEMV per  *
!*  time step,  psi = prefactor * sum_phi * g.                          *
!*                                                                     *
!*  With ix blocked by i and running over j inside, that expression is  *
!*  a symmetric block-Toeplitz matrix with symmetric Toeplitz blocks:   *
!*  block (i,ip) is Toeplitz in (j,jp) with first column kern(:,i,ip),  *
!*  and kern(:,i,ip) = kern(:,ip,i) makes the block structure symmetric *
!*  as well.  The dense matrix therefore holds only npar*nperp^2        *
!*  distinct numbers in nperp^2*npar^2 slots -- a factor npar of pure   *
!*  redundancy (11.9 GiB against 61 MiB at 200x200), and the DGEMV had  *
!*  to stream all 11.9 GiB from RAM once per time step.                 *
!*                                                                     *
!*  METHOD                                                             *
!*  ------                                                             *
!*  Each Toeplitz block is embedded in a circulant of length m and      *
!*  applied by FFT.  Writing c = kern(:,i,ip) (length npar), the        *
!*  circulant first column is                                           *
!*                                                                     *
!*      emb(d)   = c(d)      d = 0 .. npar-1                            *
!*      emb(m-d) = c(d)      d = 1 .. npar-1                            *
!*      emb      = 0         elsewhere                                  *
!*                                                                     *
!*  which requires m >= 2*npar-2.  emb is symmetric under d -> m-d, so  *
!*  its transform is REAL; only m/2+1 frequencies are independent.      *
!*  khat is therefore stored as a REAL array (nperp,nperp,0:nf-1),      *
!*  82 MiB at 200x200 against 11.9 GiB -- a factor ~150.                *
!*                                                                     *
!*  Per matvec:                                                        *
!*    - forward FFT of each of the nperp columns of g   (zero-padded)   *
!*    - one small real GEMM per frequency (khat is real, ghat complex)  *
!*    - inverse FFT of each of the nperp result columns                 *
!*                                                                     *
!*  Cost falls from O(nperp^2 npar^2) memory-bound work to              *
!*  O(nperp^2 npar) compute plus O(nperp npar log npar) transforms.     *
!*                                                                     *
!*  m is rounded up to a power of two so a compact self-contained       *
!*  radix-2 FFT can be used; MKL's DFTI would need mkl_dfti.f90 added   *
!*  to the project from the oneAPI install tree, which would tie the    *
!*  build to a particular installation path.                            *
!*                                                                     *
!*  13/08/2026                                                          *
!***********************************************************************

module mod_phi_kernel

  use shared_grid          ! nperp, npar, phi_kern, phi_khat, phi_m, phi_nf

  implicit none

  integer, parameter, private :: dp = kind(1.0d0)

contains

!-----------------------------------------------------------------------
! Smallest power of two >= n (the radix-2 FFT length).
!-----------------------------------------------------------------------
  integer function next_pow2(n)
    integer, intent(in) :: n
    next_pow2 = 1
    do while (next_pow2 < n)
      next_pow2 = 2*next_pow2
    end do
  end function next_pow2

!-----------------------------------------------------------------------
! In-place iterative radix-2 Cooley-Tukey FFT of length n (a power of 2).
! isign = -1 forward, +1 inverse (unnormalised; the caller divides by n).
!-----------------------------------------------------------------------
  subroutine fft_radix2(a, n, isign)
    integer, intent(in)       :: n, isign
    complex(dp), intent(inout) :: a(0:n-1)

    integer     :: i, j, k, m2, step
    real(dp)    :: theta
    complex(dp) :: w, wm, t, u

    ! bit-reversal permutation
    j = 0
    do i = 1, n-1
      k = n/2
      do while (k <= j)
        j = j - k
        k = k/2
      end do
      j = j + k
      if (i < j) then
        t    = a(i)
        a(i) = a(j)
        a(j) = t
      end if
    end do

    step = 1
    do while (step < n)
      m2    = 2*step
      theta = dble(isign) * 3.141592653589793238462643d0 / dble(step)
      wm    = cmplx(cos(theta), sin(theta), dp)
      do k = 0, step-1
        w = cmplx(cos(dble(k)*theta), sin(dble(k)*theta), dp)
        do i = k, n-1, m2
          t        = w * a(i+step)
          u        = a(i)
          a(i)     = u + t
          a(i+step) = u - t
        end do
      end do
      step = m2
    end do

  end subroutine fft_radix2

!-----------------------------------------------------------------------
! Build phi_khat from phi_kern.  Called once, after the kernel is either
! computed or loaded from cache.
!-----------------------------------------------------------------------
  subroutine phi_kernel_transform

    integer     :: i, ip, d, f
    complex(dp), allocatable :: buf(:)

    phi_m  = next_pow2(max(2*npar - 2, 2))
    phi_nf = phi_m/2 + 1

    if (allocated(phi_khat)) deallocate(phi_khat)
    allocate(phi_khat(nperp, nperp, 0:phi_nf-1))
    allocate(buf(0:phi_m-1))

    do i = 1, nperp
      do ip = 1, i                      ! kern(:,i,ip) = kern(:,ip,i)
        buf = (0.0_dp, 0.0_dp)
        do d = 0, npar-1
          buf(d) = cmplx(phi_kern(d, i, ip), 0.0_dp, dp)
        end do
        do d = 1, npar-1                ! wrapped upper triangle
          buf(phi_m-d) = cmplx(phi_kern(d, i, ip), 0.0_dp, dp)
        end do

        call fft_radix2(buf, phi_m, -1)

        ! emb is symmetric under d -> m-d, so the transform is real; the
        ! imaginary part is round-off and is discarded deliberately.
        do f = 0, phi_nf-1
          phi_khat(i,  ip, f) = dble(buf(f))
          phi_khat(ip, i,  f) = dble(buf(f))
        end do
      end do
    end do

    deallocate(buf)

  end subroutine phi_kernel_transform

!-----------------------------------------------------------------------
! psi = sum_phi * g, without ever forming sum_phi.
!
! g and psi are indexed as ix = (i-1)*npar + j, i.e. laid out (j, i).
!-----------------------------------------------------------------------
  subroutine phi_kernel_matvec(g, psi)

    double precision, intent(in)  :: g(nperp*npar)
    double precision, intent(out) :: psi(nperp*npar)

    complex(dp), allocatable :: buf(:), ghat(:,:), phat(:,:)
    double precision, allocatable :: br(:,:), cr(:,:)
    integer :: i, ip, j, f, ix

    allocate(ghat(nperp, 0:phi_nf-1), phat(nperp, 0:phi_nf-1))

    ! Three independent sweeps (over v_perp columns, frequencies, v_perp
    ! rows), each parallel; the implicit barrier after each !$OMP DO orders
    ! them.  buf, br and cr are per-thread work arrays.  dgemm called inside
    ! the parallel region runs single-threaded in MKL, as it should here.
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(buf, br, cr, i, ip, j, f, ix)
    allocate(buf(0:phi_m-1))
    allocate(br(nperp, 2), cr(nperp, 2))

    ! ---- forward transform of each v_perp column of g -----------------
    !$OMP DO SCHEDULE(STATIC)
    do ip = 1, nperp
      buf = (0.0_dp, 0.0_dp)
      do j = 1, npar
        buf(j-1) = cmplx(g((ip-1)*npar + j), 0.0_dp, dp)
      end do
      call fft_radix2(buf, phi_m, -1)
      do f = 0, phi_nf-1
        ghat(ip, f) = buf(f)
      end do
    end do
    !$OMP END DO

    ! ---- one small real GEMM per frequency ----------------------------
    ! khat is real, ghat complex: apply it to the real and imaginary
    ! parts together as a single nperp x nperp by nperp x 2 product.
    !$OMP DO SCHEDULE(STATIC)
    do f = 0, phi_nf-1
      do ip = 1, nperp
        br(ip, 1) = dble(ghat(ip, f))
        br(ip, 2) = aimag(ghat(ip, f))
      end do
      call dgemm('N', 'N', nperp, 2, nperp, 1.0d0, phi_khat(1,1,f), nperp, &
                 br, nperp, 0.0d0, cr, nperp)
      do i = 1, nperp
        phat(i, f) = cmplx(cr(i,1), cr(i,2), dp)
      end do
    end do
    !$OMP END DO

    ! ---- inverse transform, restoring the conjugate-even spectrum -----
    !$OMP DO SCHEDULE(STATIC)
    do i = 1, nperp
      do f = 0, phi_nf-1
        buf(f) = phat(i, f)
      end do
      do f = phi_nf, phi_m-1
        buf(f) = conjg(phat(i, phi_m-f))
      end do
      call fft_radix2(buf, phi_m, +1)
      do j = 1, npar
        ix = (i-1)*npar + j
        psi(ix) = dble(buf(j-1)) / dble(phi_m)
      end do
    end do
    !$OMP END DO

    deallocate(buf, br, cr)
    !$OMP END PARALLEL

    deallocate(ghat, phat)

  end subroutine phi_kernel_matvec

end module mod_phi_kernel
