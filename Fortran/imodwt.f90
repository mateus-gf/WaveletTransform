subroutine imodwt(W_in, V_in, N, scale_level, P, h_tilde, g_tilde, V_out)
  implicit none
  ! Declare parameters and variables
  integer :: N, scale_level, P, time_step, filter_index, index
  double precision :: W_in(N), V_in(N), h_tilde(P), g_tilde(P), V_out(N)

  ! Main loop for inverse maximal overlap discrete wavelet transform
  do time_step = 1, N
    index = time_step
    V_out(time_step) = h_tilde(1) * W_in(index) + g_tilde(1) * V_in(index)
    
    do filter_index = 2, P
      index = index + int(2**(scale_level - 1))
      if (index > N) index = index - N
      V_out(time_step) = V_out(time_step) + h_tilde(filter_index) * W_in(index) + g_tilde(filter_index) * V_in(index)
    end do
  end do
end subroutine imodwt
