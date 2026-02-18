module readHnL
    public :: read_ham

contains

! function to read the hamiltonain files
subroutine read_ham(Ham_re_prefix, Ham_re_suffix, Ham_im_prefix, Ham_im_suffix, path, n_files_available, &
                n_files, states_in_Ham, n_states, state_low, state_high, ham, lindblad)

    implicit none

    ! general user input defined
    integer, intent(in) :: n_files_available, n_files, states_in_Ham, n_states, state_low, state_high
    character (len = 100), intent(in) :: Ham_re_prefix, Ham_im_prefix, Ham_re_suffix, Ham_im_suffix, path

    ! The diagonal elements of Ham_re file are the real part of the actual hamiltonian
    ! The non-diagonal elements of the Ham_im file are the imaginary part of the actual hamiltonian
    real, dimension(n_files, n_states, n_states) :: ham(n_files, n_states, n_states)
    real, dimension(n_files, n_states, n_states) :: lindblad(n_files, n_states, n_states)

    ! since all the Ham files are considered to be in rydberg unit, they are converted to eV unit
    real :: rydberg_to_eV = 13.6056980659

    ! Defining the dummy variables that are only relevant to this function
    integer :: current_file, state, point, state1, state2
    real, dimension(states_in_Ham) :: dummy1(states_in_Ham), dummy2(states_in_Ham)
    character (len = 120) :: file_re, file_im, between

    do current_file = 1, n_files_available

        ! converting the 'current_file' type to character
        write(between, '(I5)') current_file

        ! trim(adjustl()) will remove the extra space of the characters
        file_re = trim(adjustl(path))//trim(adjustl(Ham_re_prefix))//trim(adjustl(between))//trim(adjustl(Ham_re_suffix))
        file_im = trim(adjustl(path))//trim(adjustl(Ham_im_prefix))//trim(adjustl(between))//trim(adjustl(Ham_im_suffix))

        ! unit = 6 is not accepted by fortran
        ! opening both the files, one by one
        open(current_file+7, file = file_re, action = 'read')
        open(current_file+100000, file = file_im, action = 'read')

        ! iterating over all the rows of both the files
        do state = 1, states_in_Ham

            ! reading all the rows in both the file one by one; row stored in dummy1 and dummy2
            read(current_file+7, *) dummy1
            read(current_file+100000, *) dummy2

            ! only stopping at the row corresponding to the required state
            if (state .EQ. state_low) then

                ! accept only the columns in the row corresponding to required state
                ham(current_file, 1, 1) = dummy1(state_low) * rydberg_to_eV
                ham(current_file, 1, 2) = 0.0000

                lindblad(current_file, 1, 1) = 0.0000
                lindblad(current_file, 1, 2) = dummy2(state_high) * rydberg_to_eV

            else if (state .EQ. state_high) then

                ham(current_file, 2, 2) = dummy1(state_high) * rydberg_to_eV
                ham(current_file, 2, 1) = 0.0000

                lindblad(current_file, 2, 2) = 0.0000
                lindblad(current_file, 2, 1) = dummy2(state_low) * rydberg_to_eV

            end if
        end do

        ! closing both the files without modifying them
        close(current_file+7, status = 'keep')
        close(current_file+100000, status = 'keep')

    end do

    ! when the files have to be repeated
    if (n_files .GT. n_files_available) then
        
        ! n_files_available are already read, storing data beyond that
        do current_file = (n_files_available+1), n_files
           
             point = mod(current_file, n_files_available)

             ! because ham0 does not exist, we assume
             if (point .EQ. 0) then
                 point = n_files_available
             end if

             do state1 = 1, n_states
                 do state2 = 1, n_states
            
                     ham(current_file, state1, state2) = ham(point, state1, state2)
                     lindblad(current_file, state1, state2) = lindblad(point, state1, state2)
                
                end do
            end do
        
        end do
   
    end if
end subroutine

end module readHnL