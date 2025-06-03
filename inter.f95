module inter
    public :: interpolate

contains

! function to intepolate the ham files between two nuclear time steps
subroutine interpolate(ham, lindblad, nelect, n_states, dtnucl, n_files, dtelec, ham_intrpl, lindblad_intrpl)

    implicit none 

    ! general input for the function that stays constant throughout
    integer, intent(in) :: nelect, n_states, dtnucl, n_files
    real, intent(in) :: dtelec

    ! ham is the hamiltonian read from the input files
    real, intent(in) :: ham(n_files, n_states, n_states)
    real, intent(in) :: lindblad(n_files, n_states, n_states)

    ! Note : interpolation only upto (n_files - 1) is possible
    ! Example, for the nuclear steps 1 to 2, the electronic steps are
    ! stored from 0 to 100 (or 1000 or nelect), instead of 101 to 200
    ! (if nelect = 100)
    real, intent(out) :: ham_intrpl((n_files-1)*nelect, n_states, n_states)
    real, intent(out) :: lindblad_intrpl((n_files-1)*nelect, n_states, n_states)

    ! telec is the time elapsed since the last nuclear step
    real :: telec, diff_re, diff_im

    ! dummy variables relevant only to this function
    integer :: electronic_step, nuclear_step, state1, state2, istep

    ! for beyond the files containing different nuclear steps, using
    ! forward difference method
    do nuclear_step = 1, n_files-1
        telec = 0
        do electronic_step = 1, nelect
            istep = (nuclear_step-1)*nelect + electronic_step
            do state1 = 1, n_states
              do state2 = 1, n_states

                diff_re = ham(nuclear_step + dtnucl, state1, state2) - ham(nuclear_step , state1, state2)
                ham_intrpl(istep, state1, state2) = ham(nuclear_step, state1, state2) + (diff_re/real(dtnucl))*telec
                                                    
                diff_im = lindblad(nuclear_step + dtnucl, state1, state2) - lindblad(nuclear_step , state1, state2)
                lindblad_intrpl(istep, state1, state2) = lindblad(nuclear_step, state1, state2) + (diff_im/real(dtnucl))*telec                                 

              end do
            end do
          telec = telec + dtelec
        end do
    end do

end subroutine

end module inter
