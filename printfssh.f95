module printfssh
    public func

contains

subroutine print_se(nsteps, nelect, sepop, filename, n_states)

    implicit none
    integer, intent(in) :: nsteps, nelect, n_states
    integer :: istep
    real, dimension(nsteps*nelect, n_states) :: sepop(nsteps*nelect, n_states)
    real, dimension(nsteps*nelect) :: time(nsteps*nelect)
    character(len=*), intent(in) :: filename

    open(1, file = filename, status = 'new')
        do istep = 1, nsteps*nelect
            if (istep .EQ. 1) then
                time(istep) = 0.0
            else
                time(istep) = real(istep)/real(nelect)
            end if
            write(1,*) time(istep), sepop(istep, 1), sepop(istep, 2)
        end do
    close(1)

end subroutine print_se

subroutine print_sh(nsteps, dtnucl, shpop, filename, n_states)

    implicit none
    integer, intent(in) :: nsteps, dtnucl, n_states
    integer :: istep
    real, dimension(nsteps, n_states) :: shpop(nsteps, n_states)
    real, dimension(nsteps) :: time(nsteps)
    character(len=*), intent(in) :: filename

    open(1213, file = filename, status = 'new')
        do istep = 1, nsteps
            if (istep .EQ. 1) then
                time(istep) = 0.0
            else
                time(istep) = real(istep)/real(dtnucl)
            end if

            write(1213,*) time(istep), shpop(istep, 1), shpop(istep, 2)

        end do

    close(1213)

end subroutine print_sh
end module printfssh