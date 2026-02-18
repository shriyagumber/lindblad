program main

    use readHnL       ! calling the module where all the required functions are define
    use inter         ! use inter2 for central difference method 
    use solve_me2
    use printfssh
    use solve_tdse
    use fssh

    implicit none

    ! the required information to stay constant throughout
    integer, parameter :: n_files_available = 2           ! number of Ham files available, Ham0 will not be considered
    integer, parameter :: n_files = 15005                    ! total number of files needed for dynamics
                                                            ! should always be 1 more than nsteps for dynamics
                                                            ! should be high enough considering all the ICs
    integer, parameter :: states_in_Ham = 22                ! total number of states in Ham file
    integer, parameter :: n_states = 2                      ! total number of states I want
    integer, parameter :: state_low = 16                    ! for a two state system, state number 1
    integer, parameter :: state_high = 17                   ! for a two state system, state number 1
    integer, parameter :: dtnucl = 1                        ! electronic time step
    integer, parameter :: nelect = 100                      ! number of electronic steps per nuclear step
    integer, parameter :: nsteps = 500                      ! perform dynamics for nsteps total nuclear steps
    real, parameter :: temperature = 300                    ! temperature, required for boltzmann
    integer, parameter :: realizations = 100                  ! number of stochastic realizations per initial condition
    integer, parameter :: boltz = 1                         ! 1 if you need boltzmann factor in fssh/dish
    integer, parameter :: istate = 2                        ! where the electron is initially
    real, parameter :: dtelec = real(dtnucl)/real(nelect)   ! electronic time step
    real, parameter :: gamma = 1.0                          ! for lindblad equation

    integer :: i = 1, file_no = 1, lenic, current_ic, re, total_steps, istep, j, k, m, n
    integer, dimension(:), allocatable :: ic                  ! array for initial conditions

    ! initializing the variables
    character (len = 100) :: Ham_re_prefix, Ham_im_prefix, Ham_re_suffix, Ham_im_suffix, path

    ! initializing the ham and interpolated ham 3-dimensional arrays
    real, dimension(n_files, n_states, n_states) :: ham(n_files, n_states, n_states)
    real, dimension(n_files, n_states, n_states) :: lindblad(n_files, n_states, n_states)

    ! interpolate the hamiltonian
    real, dimension((n_files-1)*nelect, n_states, n_states) :: ham_intrpl((n_files-1)*nelect, n_states, n_states)
    real, dimension((n_files-1)*nelect, n_states, n_states) :: lindblad_intrpl((n_files-1)*nelect, n_states, n_states)

    ! dimensions equal to sstens*nelect for each IC
    real, dimension(nsteps*nelect, n_states, n_states) :: hamic(nsteps*nelect, n_states, n_states)
    real, dimension(nsteps*nelect, n_states, n_states) :: lindbladic(nsteps*nelect, n_states, n_states)

    ! dimensions equal to ham=nsteps*nelect
    real, dimension(nsteps*nelect, n_states, n_states) :: den_lindblad(nsteps*nelect, n_states, n_states)
    real, dimension(nsteps*nelect, n_states, n_states) :: den_tdse(nsteps*nelect, n_states, n_states)
    real, dimension(nsteps*nelect, n_states) :: sepop_lindblad(nsteps*nelect, n_states)
    real, dimension(nsteps*nelect, n_states) :: sepop_tdse(nsteps*nelect, n_states)
    real, dimension(nsteps*nelect, n_states) :: seavg_lindblad(nsteps*nelect, n_states)
    real, dimension(nsteps*nelect, n_states) :: seavg_tdse(nsteps*nelect, n_states) 
    
    real, dimension(nsteps, n_states) :: shpopu(nsteps, n_states)
    real, dimension(nsteps, n_states) :: shavg(nsteps, n_states)
    integer, dimension(nsteps, n_states) :: sh(nsteps, n_states)

    total_steps = nsteps * nelect

    ! The Ham files are assumed to be in rydberg units
    path = "./avg-ham/"                                         ! path where the ham files are stored
    
    ! suffix and prefix of the names of the ham files
    Ham_re_prefix = "Ham"
    Ham_re_suffix = "_re"
    Ham_im_prefix = "Ham"
    Ham_im_suffix = "_im"

    ! stating the initial conditions
    allocate (ic(0))
    do while (file_no .le. n_files_available)
        ic = [ic, file_no]
        file_no  = file_no+1
    end do

    lenic = size(ic)

    sepop_lindblad(:,:) = 0.0
    sepop_tdse(:,:) = 0.0
    shpopu(:,:) = 0.0


    ! The below part will remain constant, DO NOT CHANGE ANYTHING
    ! Reading the ham files
    call read_ham(Ham_re_prefix, Ham_re_suffix, Ham_im_prefix, Ham_im_suffix, path, n_files_available, &
                    n_files, states_in_Ham, n_states, state_low, state_high, ham, lindblad)

    ! interpolating the ham files in a 3d array
    call interpolate(ham, lindblad, nelect, n_states, dtnucl, n_files, dtelec, ham_intrpl, lindblad_intrpl)

    ! to average over all the initial conditions
    do i = 1, lenic

        current_ic = ic(i)

        ! The hamiltonian corresponding to the initial condition is selected from the part of TD ham       
        ! since the ham between 2 to 1 is stored from 1 to 1000 and btwn 3 to 2
        ! between 1001 to 2000, the (current_ic - 1) is done, + 1 is done
        ! because Ham0 is not included. To be consistent with the size, -2 is
        ! done at the end
        hamic(:,:,:) = ham_intrpl(((current_ic-1)*nelect+1):(((current_ic-1)*nelect)+(nsteps*nelect)-2), :, :)
        lindbladic(:,:,:) = lindblad_intrpl(((current_ic-1)*nelect+1):(((current_ic-1)*nelect)+(nsteps*nelect)-2), :, :)

        do re = 1, realizations
            ! Evolving TDSE for SE population
            ! den = tdse(hamic, n_states, nelect, dtelec, n_files, nsteps, istate)
            ! sepopu = sepopu + den

            ! Evolving Lindblad ME
            call me(hamic, lindbladic, boltz, temperature, n_states, nelect, dtelec, n_files, nsteps, istate, gamma, den_lindblad)
            call tdse(hamic, lindbladic, n_states, nelect, dtelec, n_files, nsteps, istate, den_tdse)  
            call fsshpop(den_tdse, hamic, n_states, nelect, nsteps, temperature, boltz, istate, sh)

            ! extract and accumulate populations (the diagonal of den)
            do istep = 1, total_steps
                do j = 1, n_states
                    sepop_lindblad(istep,j) = sepop_lindblad(istep,j) + den_lindblad(istep,j,j)
                    sepop_tdse(istep,j) = sepop_tdse(istep,j) + den_tdse(istep,j,j)
                end do
            end do

            do m = 1, nsteps
                do n = 1, n_states
                    shpopu(m,n) = shpopu(m,n) + sh(m,n)
                end do
            end do
        end do
    end do

    seavg_lindblad = sepop_lindblad/(lenic*realizations)
    seavg_tdse = sepop_tdse/(lenic*realizations)
    shavg = shpopu/(lenic*realizations)

    call print_se(nsteps, nelect, seavg_lindblad, "lindblad_pop.out", n_states)
    call print_se(nsteps, nelect, seavg_tdse, "tdse_pop.out", n_states)
    call print_sh(nsteps, dtnucl, shavg, "sh_pop.out", n_states)

end program main