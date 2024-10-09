subroutine dwt(V_in, N, P, h_coef, g_coef, W_out, V_out)
  implicit none
  ! Declare parameters and variables
  integer :: N, P, time_step, filter_index, index
  double precision :: V_in(N), h_coef(P), g_coef(P), W_out(N/2), V_out(N/2)

  ! Main loop for discrete wavelet transform
  do time_step = 1, N/2
    index = 2 * time_step
    W_out(time_step) = h_coef(1) * V_in(index)
    V_out(time_step) = g_coef(1) * V_in(index)
    
    do filter_index = 2, P
      index = index - 1
      if (index < 1) index = N
      W_out(time_step) = W_out(time_step) + h_coef(filter_index) * V_in(index)
      V_out(time_step) = V_out(time_step) + g_coef(filter_index) * V_in(index)
    end do
  end do
end subroutine dwt
