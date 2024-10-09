subroutine idwt(W_in, V_in, N, P, h_coef, g_coef, X_out)
  implicit none
  ! Declare parameters and variables
  integer :: N, P, time_step, idx_h, idx_g, index, pos, new_m, new_n
  double precision :: W_in(N/2), V_in(N/2), h_coef(P), g_coef(P), X_out(N)

  ! Initialize new_m and new_n
  new_m = -1
  new_n = 0

  ! Main loop for inverse discrete wavelet transform
  do time_step = 1, N
    new_m = new_m + 2
    new_n = new_n + 2
    pos = time_step
    idx_h = 2
    idx_g = 1
    X_out(new_m) = h_coef(idx_h) * W_in(pos) + g_coef(idx_h) * V_in(pos)
    X_out(new_n) = h_coef(idx_g) * W_in(pos) + g_coef(idx_g) * V_in(pos)
    
    if (P > 2) then
      do index = 1, P/2 - 1
        pos = pos + 1
        if (pos > N) pos = 1
        idx_h = idx_h + 2
        idx_g = idx_g + 2
        X_out(new_m) = X_out(new_m) + h_coef(idx_h) * W_in(pos) + g_coef(idx_h) * V_in(pos)
        X_out(new_n) = X_out(new_n) + h_coef(idx_g) * W_in(pos) + g_coef(idx_g) * V_in(pos)
      end do
    end if
  end do
end subroutine idwt
