subroutine modwt(V_in, N, scale_level, P, h_tilde, g_tilde, W_out, V_out)
  implicit none
  ! Declare parameters and variables
  integer :: N, scale_level, P, time_step, filter_index, index
  double precision :: V_in(N), h_tilde(P), g_tilde(P), W_out(N), V_out(N)

  ! Main loop for maximal overlap discrete wavelet transform
  do time_step = 1, N
    index = time_step
    W_out(time_step) = h_tilde(1) * V_in(index)
    V_out(time_step) = g_tilde(1) * V_in(index)
    
    do filter_index = 2, P
      index = index - int(2**(scale_level - 1))
      if (index < 1) index = index + N
      W_out(time_step) = W_out(time_step) + h_tilde(filter_index) * V_in(index)
      V_out(time_step) = V_out(time_step) + g_tilde(filter_index) * V_in(index)
    end do
  end do
end subroutine modwt
