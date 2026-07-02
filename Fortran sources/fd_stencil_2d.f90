MODULE mod_fd_stencil_2d

  USE derivatives_2d, ONLY: dp, fornberg_weights
  USE shared_beam,   ONLY: taum
  USE shared_grid,   ONLY: i_upwind

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: fd_stencil_2d

CONTAINS

SUBROUTINE fd_stencil_2d(i, j, nperp, npar, vperp, dvpar, &
                          A_ij, B_ij, C_ij, D_ij, E_ij, F_ij, &
                          col_idx, coeff, n_entries, rhs)
  !-----------------------------------------------------------------
  ! Assembles the 7-point Fornberg finite-difference stencil for:
  !
  !   A*f + B*df/dvp + C*df/dva + D*d2f/dvp2
  !       + E*d2f/dvp_dva + F*d2f/dva2 - taum*f = 0
  !
  ! v⊥ direction : non-uniform grid, Fornberg 7-point stencil
  ! v∥ direction : uniform spacing dvpar, Fornberg 7-point stencil
  !
  ! taum         : time/damping term from shared_timer, subtracted
  !                from the centre diagonal of all interior rows.
  !
  ! Boundary conditions:
  !   i=nperp, j=1, j=npar : f = 0       (Dirichlet)
  !   i=1                   : df/dvp = 0  (Neumann, unit-scaled)
  !
  ! OUTPUT:
  !   col_idx(k), coeff(k) : sparse row entries (at most 7*7 = 49)
  !   n_entries             : number of non-zero entries
  !   rhs                   : RHS contribution (0, homogeneous BCs)
  !-----------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,  INTENT(IN)  :: i, j, nperp, npar
  REAL(dp), INTENT(IN)  :: vperp(nperp), dvpar
  REAL(dp), INTENT(IN)  :: A_ij, B_ij, C_ij, D_ij, E_ij, F_ij
  INTEGER,  INTENT(OUT) :: col_idx(49)
  REAL(dp), INTENT(OUT) :: coeff(49)
  INTEGER,  INTENT(OUT) :: n_entries
  REAL(dp), INTENT(OUT) :: rhs

  INTEGER  :: Li, Lj, k, l, ip, jp, mi, mj
  REAL(dp) :: xloc_i(7), xloc_j(7)
  REAL(dp) :: wi(7,0:2), wj(7,0:2)
  REAL(dp) :: acc(7,7)
  REAL(dp) :: wmax
  REAL(dp) :: wb(7), dvp, dvm, dv_loc, Pe_perp

  ! Initialise outputs
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
    rhs        = 0.0_dp
    RETURN
  END IF

  !================================================================
  ! CASE 2: Neumann BC at vperp=0 (i=1)  ->  df/dvp = 0
  !
  ! One-sided 7-point Fornberg stencil at node i=1.
  ! Weights normalised by their maximum so this row has the same
  ! order-of-magnitude scaling as the build_ss convention
  ! (-1.5, 2.0, -0.5).
  !================================================================
  IF (i == 1) THEN
    DO k = 1, 7
      xloc_i(k) = vperp(k) - vperp(1)
    END DO
    CALL fornberg_weights(0.0_dp, xloc_i, 7, 1, wi)

    wmax = MAXVAL(ABS(wi(:,1)))
    IF (wmax > 0.0_dp) wi(:,1) = wi(:,1) / wmax

    n_entries = 0
    DO k = 1, 7
      IF (wi(k,1) == 0.0_dp) CYCLE
      n_entries          = n_entries + 1
      col_idx(n_entries) = (k-1)*npar + j
      coeff(n_entries)   = wi(k,1)
    END DO
    rhs = 0.0_dp
    RETURN
  END IF

  !================================================================
  ! CASE 3: Interior points — full PDE stencil
  !================================================================

  !--- v⊥ Fornberg weights (non-uniform) ---------------------------
  Li = MAX(1, i-3)
  Li = MIN(Li, nperp-6)
  DO k = 1, 7
    xloc_i(k) = vperp(Li+k-1) - vperp(i)
  END DO
  CALL fornberg_weights(0.0_dp, xloc_i, 7, 2, wi)

  !--- v∥ Fornberg weights (uniform) --------------------------------
  Lj = MAX(1, j-3)
  Lj = MIN(Lj, npar-6)
  DO k = 1, 7
    xloc_j(k) = REAL(Lj+k-1-j, dp) * dvpar
  END DO
  CALL fornberg_weights(0.0_dp, xloc_j, 7, 2, wj)

  !--- Local index of (i,j) within stencils ------------------------
  mi = i - Li + 1
  mj = j - Lj + 1

  !--- v⊥ drag (B) weights: central by default; Peclet-hybrid upwind
  !    when i_upwind=-1 and the local cell-Peclet |B|dv/D > 2.  The
  !    upwind side is chosen so the semi-discrete eigenvalue has Re<=0
  !    (B>0 => advection speed -B<0 => forward difference, and vice
  !    versa).  This restores dissipativity of the advection-dominated
  !    high-v⊥ boundary layer without changing the sparsity pattern
  !    (the D 2nd-derivative stencil already fills all same-j columns).
  wb(:) = wi(:,1)
  IF (i_upwind == -1) THEN
    dvp     = vperp(i+1) - vperp(i)
    dvm     = vperp(i)   - vperp(i-1)
    dv_loc  = MIN(dvp, dvm)
    Pe_perp = 0.0_dp
    IF (D_ij /= 0.0_dp) Pe_perp = ABS(B_ij) * dv_loc / ABS(D_ij)
    IF (Pe_perp > 2.0_dp) THEN
      wb(:) = 0.0_dp
      IF (B_ij >= 0.0_dp) THEN        ! forward (upwind) difference
        wb(mi)   = -1.0_dp / dvp
        wb(mi+1) = +1.0_dp / dvp
      ELSE                            ! backward (upwind) difference
        wb(mi)   = +1.0_dp / dvm
        wb(mi-1) = -1.0_dp / dvm
      END IF
    END IF
  END IF

  !--- Accumulate PDE terms ----------------------------------------
  acc = 0.0_dp

  DO k = 1, 7
    DO l = 1, 7

      ! A * f  (centre only)
      IF (k == mi .AND. l == mj) &
        acc(k,l) = acc(k,l) + A_ij

      ! B * df/dvp  (v⊥ 1st deriv, fixed j; wb = central or upwind)
      IF (l == mj) &
        acc(k,l) = acc(k,l) + B_ij * wb(k)

      ! C * df/dvpa  (v∥ 1st deriv, fixed i)
      IF (k == mi) &
        acc(k,l) = acc(k,l) + C_ij * wj(l,1)

      ! D * d2f/dvp2  (v⊥ 2nd deriv, fixed j)
      IF (l == mj) &
        acc(k,l) = acc(k,l) + D_ij * wi(k,2)

      ! E * d2f/dvp_dvpa  (outer product of 1st derivs)
      acc(k,l) = acc(k,l) + E_ij * wi(k,1) * wj(l,1)

      ! F * d2f/dvpa2  (v∥ 2nd deriv, fixed i)
      IF (k == mi) &
        acc(k,l) = acc(k,l) + F_ij * wj(l,2)

    END DO
  END DO

  !--- Subtract taum from centre diagonal --------------------------
  acc(mi,mj) = acc(mi,mj) - taum

  !--- Zero Dirichlet boundary node contributions ------------------
  DO k = 1, 7
    ip = Li + k - 1
    IF (ip == nperp) acc(k,:) = 0.0_dp
  END DO
  DO l = 1, 7
    jp = Lj + l - 1
    IF (jp == 1 .OR. jp == npar) acc(:,l) = 0.0_dp
  END DO

  !--- Flatten to sparse output ------------------------------------
  DO k = 1, 7
    ip = Li + k - 1
    IF (ip < 1 .OR. ip > nperp) CYCLE
    DO l = 1, 7
      jp = Lj + l - 1
      IF (jp < 1 .OR. jp > npar) CYCLE
      IF (acc(k,l) == 0.0_dp) THEN
        ! On-axis entries (k==mi or l==mj) only appear when A/B/C/D/F are
        ! non-zero; safe to drop.  Boundary-zeroed entries are also dropped.
        ! Off-axis entries (k/=mi, l/=mj) are contributed to solely by E_ij:
        ! when E_ij=0 their acc is 0, but the pattern was built with sentinel
        ! E_ij=1 so they ARE in ia_L/ja_L.  Dropping them here shifts ptr for
        ! every subsequent row and corrupts aa_L entirely.  Keep them (value 0).
        IF (k == mi .OR. l == mj) CYCLE
        IF (ip == nperp .OR. jp == 1 .OR. jp == npar) CYCLE
      END IF
      n_entries          = n_entries + 1
      col_idx(n_entries) = (ip-1)*npar + jp
      coeff(n_entries)   = acc(k,l)
    END DO
  END DO

END SUBROUTINE fd_stencil_2d

END MODULE mod_fd_stencil_2d
