MODULE mod_fd_stencil_2d_3
  !-----------------------------------------------------------------
  ! 3-point (2nd-order) finite-difference stencil for the 2D FP
  ! operator in cylindrical velocity space (v⊥, v∥).
  !
  ! Companion to mod_fd_stencil_2d (7-point Fornberg stencil).
  ! Used by timefp_upd (CSR/PARDISO-optimised 5-point time solver).
  !
  ! v⊥ : non-uniform grid, 3-point 2nd-order formula
  ! v∥ : uniform grid, centred 3-point formula
  !
  ! FD weights:
  !
  !   v⊥ (h_m = vp(i)-vp(i-1), h_p = vp(i+1)-vp(i)):
  !     df/dvp:   alpha_l = -h_p/(h_m*(h_m+h_p))
  !               alpha_n = (h_p-h_m)/(h_m*h_p)
  !               alpha_r =  h_m/(h_p*(h_m+h_p))
  !     d2f/dvp2: beta_l  =  2/(h_m*(h_m+h_p))
  !               beta_n  = -2/(h_m*h_p)
  !               beta_r  =  2/(h_p*(h_m+h_p))
  !
  !   v∥ (uniform dvpar):
  !     df/dva:   wpa1_l = -1/(2*dvpar),  wpa1_r = +1/(2*dvpar)
  !     d2f/dva2: wpa2_l = +1/dvpar^2, wpa2_n = -2/dvpar^2, wpa2_r = +1/dvpar^2
  !-----------------------------------------------------------------
  USE shared_beam, ONLY: taum

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: fd_stencil_2d_3

  INTEGER, PARAMETER :: dp = KIND(1.0D0)

CONTAINS

SUBROUTINE fd_stencil_2d_3(i, j, nperp, npar, vperp, dvpar, &
                             A_ij, B_ij, C_ij, D_ij, E_ij, F_ij, &
                             col_idx, coeff, n_entries, rhs)
  !-----------------------------------------------------------------
  ! Returns the sparse row for grid point (i,j) under the operator:
  !
  !   A*f + B*df/dvp + C*df/dva + D*d2f/dvp2
  !       + E*d2f/dvp_dva + F*d2f/dva2 - taum*f = 0
  !
  ! Boundary conditions:
  !   i=nperp, j=1, j=npar : f = 0       (Dirichlet)
  !   i=1                   : df/dvp = 0  (Neumann, 2nd-order one-sided)
  !
  ! At most 9 non-zero entries (3 x 3 stencil).
  ! Entries are returned in strictly ascending column order: the
  ! natural (di,dj) loop order guarantees this for the 3-point
  ! stencil when npar >= 4 (always true in practice).
  !-----------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,  INTENT(IN)  :: i, j, nperp, npar
  REAL(dp), INTENT(IN)  :: vperp(nperp), dvpar
  REAL(dp), INTENT(IN)  :: A_ij, B_ij, C_ij, D_ij, E_ij, F_ij
  INTEGER,  INTENT(OUT) :: col_idx(9)
  REAL(dp), INTENT(OUT) :: coeff(9)
  INTEGER,  INTENT(OUT) :: n_entries
  REAL(dp), INTENT(OUT) :: rhs

  REAL(dp) :: h_m, h_p
  REAL(dp) :: alpha_l, alpha_n, alpha_r
  REAL(dp) :: beta_l,  beta_n,  beta_r
  REAL(dp) :: wpa1_l, wpa1_r
  REAL(dp) :: wpa2_l, wpa2_n, wpa2_r
  ! acc(di,dj): di=1->i-1, di=2->i, di=3->i+1
  !             dj=1->j-1, dj=2->j, dj=3->j+1
  REAL(dp) :: acc(3,3)
  INTEGER  :: di, dj, ip, jp

  rhs       = 0.0_dp
  n_entries = 0
  col_idx   = 0
  coeff     = 0.0_dp

  !================================================================
  ! CASE 1: Dirichlet boundary rows  ->  1 * f(i,j) = 0
  !================================================================
  IF (i == nperp .OR. j == 1 .OR. j == npar) THEN
    n_entries  = 1
    col_idx(1) = (i-1)*npar + j
    coeff(1)   = 1.0_dp
    RETURN
  END IF

  !================================================================
  ! CASE 2: Neumann BC at vperp=0 (i=1)  ->  df/dvp = 0
  !   2nd-order one-sided (normalised coefficients):
  !   -1.5*f(1,j) + 2*f(2,j) - 0.5*f(3,j) = 0
  !================================================================
  IF (i == 1) THEN
    n_entries  = 3
    col_idx(1) = j            ! (1-1)*npar + j
    coeff(1)   = -1.5_dp
    col_idx(2) = npar + j     ! (2-1)*npar + j
    coeff(2)   =  2.0_dp
    col_idx(3) = 2*npar + j   ! (3-1)*npar + j
    coeff(3)   = -0.5_dp
    RETURN
  END IF

  !================================================================
  ! CASE 3: Interior points — full 3-point PDE stencil
  !================================================================

  !--- v⊥ weights (non-uniform grid) -----------------------------
  h_m = vperp(i)   - vperp(i-1)
  h_p = vperp(i+1) - vperp(i)

  alpha_l = -h_p          / (h_m*(h_m+h_p))
  alpha_n =  (h_p - h_m)  / (h_m*h_p)
  alpha_r =  h_m          / (h_p*(h_m+h_p))

  beta_l  =  2.0_dp / (h_m*(h_m+h_p))
  beta_n  = -2.0_dp / (h_m*h_p)
  beta_r  =  2.0_dp / (h_p*(h_m+h_p))

  !--- v∥ weights (uniform grid) ---------------------------------
  wpa1_l =  -1.0_dp / (2.0_dp*dvpar)
  wpa1_r =  +1.0_dp / (2.0_dp*dvpar)

  wpa2_l =  1.0_dp / (dvpar*dvpar)
  wpa2_n = -2.0_dp / (dvpar*dvpar)
  wpa2_r =  1.0_dp / (dvpar*dvpar)

  !--- Accumulate PDE terms into 3x3 local array -----------------
  acc = 0.0_dp

  ! A: centre only
  acc(2,2) = acc(2,2) + A_ij

  ! B: df/dvp  (v⊥ first derivative, same v∥)
  acc(1,2) = acc(1,2) + B_ij * alpha_l
  acc(2,2) = acc(2,2) + B_ij * alpha_n
  acc(3,2) = acc(3,2) + B_ij * alpha_r

  ! C: df/dva  (v∥ first derivative, centred; no centre weight)
  acc(2,1) = acc(2,1) + C_ij * wpa1_l
  acc(2,3) = acc(2,3) + C_ij * wpa1_r

  ! D: d2f/dvp2  (v⊥ second derivative, same v∥)
  acc(1,2) = acc(1,2) + D_ij * beta_l
  acc(2,2) = acc(2,2) + D_ij * beta_n
  acc(3,2) = acc(3,2) + D_ij * beta_r

  ! E: d2f/dvp_dva  (outer product of 1st-deriv weights)
  ! v∥ centred weight at dj=0 is zero, so only dj=1,3 contribute.
  ! alpha_n is non-zero for non-uniform grids, hence (i,j±1) entries.
  acc(1,1) = acc(1,1) + E_ij * alpha_l * wpa1_l
  acc(1,3) = acc(1,3) + E_ij * alpha_l * wpa1_r
  acc(2,1) = acc(2,1) + E_ij * alpha_n * wpa1_l
  acc(2,3) = acc(2,3) + E_ij * alpha_n * wpa1_r
  acc(3,1) = acc(3,1) + E_ij * alpha_r * wpa1_l
  acc(3,3) = acc(3,3) + E_ij * alpha_r * wpa1_r

  ! F: d2f/dva2  (v∥ second derivative, same v⊥)
  acc(2,1) = acc(2,1) + F_ij * wpa2_l
  acc(2,2) = acc(2,2) + F_ij * wpa2_n
  acc(2,3) = acc(2,3) + F_ij * wpa2_r

  !--- Subtract taum from centre diagonal -------------------------
  acc(2,2) = acc(2,2) - taum

  !--- Zero contributions that land on Dirichlet boundaries --------
  ! i=nperp -> Dirichlet; i=1 is Neumann (f unknown, keep)
  DO di = 1, 3
    ip = i + (di-2)
    IF (ip == nperp) acc(di,:) = 0.0_dp
  END DO
  ! j=1 and j=npar -> Dirichlet
  DO dj = 1, 3
    jp = j + (dj-2)
    IF (jp == 1 .OR. jp == npar) acc(:,dj) = 0.0_dp
  END DO

  !--- Flatten to sparse output -----------------------------------
  ! Loop order (di=1..3, dj=1..3) traverses columns in ascending
  ! order for any npar >= 4, so no sort is needed.
  DO di = 1, 3
    ip = i + (di-2)
    DO dj = 1, 3
      jp = j + (dj-2)
      IF (acc(di,dj) == 0.0_dp) CYCLE
      n_entries          = n_entries + 1
      col_idx(n_entries) = (ip-1)*npar + jp
      coeff(n_entries)   = acc(di,dj)
    END DO
  END DO

END SUBROUTINE fd_stencil_2d_3

END MODULE mod_fd_stencil_2d_3
