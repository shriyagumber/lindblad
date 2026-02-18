module solve_tdse
    public :: tdse

contains

! Hamic: ham wrt initial condition
! n_states: 2 always
! nelect: number of electronic steps per 1 nuclear step
! dtelect: electronic time-step, 
! n_files: number od original ham
! nsteps: dynamics length
! istate: initial state for electron
! To solve the time dependent Lindblad master equation and return the density matrix
subroutine tdse(hamic, lindbladic, n_states, nelect, dtelec, n_files, nsteps, istate, den)

    implicit none

    ! only the part of the hamiltonian corresponding to the initial condition is stored
    real, intent(in) :: hamic(nsteps*nelect, n_states, n_states)
    real, intent(in) :: lindbladic(nsteps*nelect, n_states, n_states)

    ! the general input for this function
    integer, intent(in) :: n_states, nelect, n_files, istate, nsteps
    real, intent(in) :: dtelec

    ! defining the density matrix
    real, intent(out) :: den(nsteps*nelect, n_states, n_states)

    ! dummy variables that are relavant only to this function
    real :: factor_re, factor_im
    integer :: state, istep

    ! — locals —
    integer :: total_steps, i, j
    real, parameter :: hbar = 0.6582119510   ! eV·fs

    complex, dimension(n_states,n_states) :: rho(n_states,n_states)
    complex, dimension(n_states,n_states) :: rho_next(n_states,n_states)

    complex, allocatable :: drho(:,:), H(:,:), L(:,:)

    ! total electronic steps
    total_steps = nsteps * nelect

    ! 1) initialize ρ(0) = |istate><istate|
    rho = (0.00, 0.00)
    rho(istate, istate) = (1.000, 0.000)
    den(1,:,:) = real(rho)
    
    ! 2) time‐propagation loop (forward Euler)
    do istep = 1, total_steps-1

        ! allocate complex matrices
        allocate(drho(n_states,n_states))
        allocate(H(n_states,n_states))
        allocate(L(n_states,n_states))

        ! 2a) build H and L at this step
        do i = 1, n_states
            do j = 1, n_states
                H(i,j) = CMPLX( hamic(istep, i, j), 0.00)
                L(i,j) = CMPLX( 0.00, lindbladic(istep, i, j))
            end do
        end do

        ! 2b) compute drho = −(i/ħ)[H,ρ] + γ (LρL† − ½{L†L,ρ})
        ! call compute_drho(rho, H, L, hbar, n_states, drho)

        ! 2c) Euler update
        ! rho_next = rho + dtelec * drho

        ! New RK4
        call rk4_tdse(rho, H, L, hbar, dtelec, n_states, rho_next)

        ! 2d) renormalise Tr(ρ) → 1
        call renormalise(n_states, rho_next)

        ! 2e) save for output and advance
        den(istep+1,:,:) = real(rho_next)
        rho = rho_next
        
        ! deallocate scratch
        deallocate(drho, H, L)

    end do

end subroutine tdse

! basicaaly, you pass it rho, and get rho_next
subroutine rk4_tdse(rho, H, L, hbar, dt, n_states, rho_next)

    implicit none
    integer, intent(in) :: n_states
    real, intent(in) :: hbar, dt
    complex, intent(in) :: rho(n_states,n_states), H(n_states,n_states), L(n_states,n_states)
    complex, intent(out) :: rho_next(n_states,n_states)

    complex, allocatable :: k1(:,:), k2(:,:), k3(:,:), k4(:,:), tmp(:,:)
    allocate(k1(n_states,n_states), k2(n_states,n_states), k3(n_states,n_states), k4(n_states,n_states), tmp(n_states,n_states))

    ! adding the diagonal energies and off-diagonal complex NACs
    call compute_commutator(H+L, rho, hbar, k1)
    tmp = rho + 0.50 * dt * k1
    call compute_commutator(H+L, tmp, hbar, k2)
    tmp = rho + 0.50 * dt * k2
    call compute_commutator(H+L, tmp, hbar, k3)
    tmp = rho + dt * k3
    call compute_commutator(H+L, tmp, hbar, k4)

    rho_next = rho + dt/6.00 * (k1 + 2*k2 + 2*k3 + k4)
    deallocate(k1, k2, k3, k4, tmp)

end subroutine

subroutine compute_commutator(H, rho, hbar, k_out)
    implicit none
    integer :: n_states
    complex, intent(in) :: H(:,:), rho(:,:)
    complex, intent(out) :: k_out(:,:)
    real, intent(in) :: hbar
    complex :: iC
    n_states = size(H, 1)
    iC = (0.0d0, 1.0d0)
    k_out = -iC / hbar * (matmul(H, rho) - matmul(rho, H))
end subroutine

! Using Euler's method
subroutine compute_drho(rho, H, L, hbar, n_states, drho)
    implicit none
    integer :: n, i, n_states
    complex, intent(in)  :: rho(n_states,n_states), H(n_states,n_states), L(n_states,n_states)
    real, intent(in)  :: hbar
    complex, intent(out) :: drho(n_states,n_states)

    ! some local variables
    complex, allocatable :: comm(:,:)
    complex :: iC

    ! so everything is of the size n_statesxn_states
    allocate(comm(n_states,n_states))
    iC = (0.000, 1.000)

    ! term1: commutator [H,ρ]
    comm = matmul(H+L, rho) - matmul(rho, H+L)
    print *, H
    print *, L

    drho = -(iC/hbar)*comm

    deallocate(comm)

end subroutine

subroutine renormalise(n_states, rho)
    implicit none

    integer, intent(in) :: n_states
    complex, intent(inout) :: rho(n_states,n_states)
    integer :: i
    real :: tr
    tr = 0.00
    do i = 1, n_states
      tr = tr + real(rho(i,i))
    end do
    rho = rho / tr
end subroutine

end module solve_tdse