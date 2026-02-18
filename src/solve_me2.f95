module solve_me2
    public :: me

contains

! Hamic: ham wrt initial condition
! n_states: 2 always
! nelect: number of electronic steps per 1 nuclear step
! dtelect: electronic time-step, 
! n_files: number od original ham
! nsteps: dynamics length
! istate: initial state for electron
! To solve the time dependent Lindblad master equation and return the density matrix
subroutine me(hamic, lindbladic, boltz, temperature, n_states, nelect, dtelec, n_files, nsteps, istate, gamma, den)

    implicit none

    ! only the part of the hamiltonian corresponding to the initial condition is stored
    real, intent(in) :: hamic(nsteps*nelect, n_states, n_states)
    real, intent(in) :: lindbladic(nsteps*nelect, n_states, n_states)

    ! the general input for this function
    integer, intent(in) :: n_states, nelect, n_files, istate, nsteps, boltz
    real, intent(in) :: dtelec
    real, intent(in) :: gamma  ! decay rate

    ! defining the density matrix
    real, intent(out) :: den(nsteps*nelect, n_states, n_states)

    ! dummy variables that are relavant only to this function
    real :: factor_re, factor_im, boltz_factor, temperature
    integer :: state, istep

    ! — locals —
    integer :: total_steps, i, j
    real, parameter :: hbar = 0.6582119510   ! eV·fs
    real, parameter :: kb = 0.0000862   ! eV/K

    complex, dimension(n_states,n_states) :: rho(n_states,n_states)
    complex, dimension(n_states,n_states) :: rho_next(n_states,n_states)

    complex, allocatable :: drho(:,:), H(:,:), Leg(:,:), nac(:,:), Lge(:,:)
    real, allocatable :: fermi_rate

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
        allocate(nac(n_states,n_states))
        allocate(Leg(n_states,n_states))
        allocate(Lge(n_states,n_states))
        allocate(fermi_rate)

        ! 2a) build H and L at this step
        do i = 1, n_states
            do j = 1, n_states
                H(i,j) = CMPLX( hamic (istep, i, j), 0.00 )

                nac(i,j) = CMPLX( 0.00, lindbladic(istep, i, j))
                ! nac(i,j) = CMPLX( 0.00, 0.00)
                
                ! ! Boltzmann factor
                ! if (boltz .EQ. 1) then
                !     if (real(H(j,j)) > real(H(i,i))) then
                !         boltz_factor = exp(-kb/temperature) * (H(j,j) - H(i,i))
                !         L(i,j) = L(i,j)*sqrt(boltz_factor/hbar)
                !     else
                !         L(i,j) = L(i,j)*sqrt(1/hbar)
                !     end if
                ! end if
            end do
        end do

        ! rate = abs(lindbladic(istep, 1, 2)/hbar)

        fermi_rate = 3.14159*lindbladic(istep, 1, 2)*lindbladic(istep, 1, 2)/((H(2,2)-H(1,1))*hbar)

        Leg(1,1)=0.0
        Leg(1,2)=0.0
        Leg(2,1)=0.0
        Leg(2,2)=0.0

        Lge(1,1)=0.0
        Lge(1,2)=sqrt(fermi_rate)
        Lge(2,1)=0.0
        Lge(2,2)=0.0

        ! manually using the boltzmann factor scaling
        !boltz_factor = exp((-kb/temperature) * (H(2,2) - H(1,1)))

        !Leg(2,1) = sqrt(1/hbar);
        !Lge(1,2) = sqrt(nac(1,2)*boltz_factor/hbar);
        !Leg(2,1) = sqrt(nac(2,1)/hbar);
        
        ! 2b) compute drho = −(i/ħ)[H,ρ] + γ (LρL† − ½{L†L,ρ})
        ! call compute_drho(rho, H, L, gamma, hbar, n_states, drho)
        call rk4_lindblad(rho, H, Lge, Leg, hbar, dtelec, n_states, rho_next)

        ! 2c) Euler update
        ! rho_next = rho + dtelec * drho

        ! 2d) renormalise Tr(ρ) → 1
        call renormalise(n_states, rho_next)

        ! 2e) save for output and advance
        den(istep+1,:,:) = real(rho_next)
        rho = rho_next
        
        ! deallocate scratch
        deallocate(drho, H, Leg, Lge, fermi_rate, nac)

    end do

end subroutine me

subroutine rk4_lindblad(rho, H, Lge, Leg, hbar, dt, n, rho_out)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: dt, hbar
    complex, intent(in) :: rho(n,n), H(n,n), Lge(n,n), Leg(n,n)
    complex, intent(out) :: rho_out(n,n)

    complex, allocatable :: k1(:,:), k2(:,:), k3(:,:), k4(:,:), tmp(:,:)

    allocate(k1(n,n), k2(n,n), k3(n,n), k4(n,n), tmp(n,n))

    call compute_drho_lindblad(rho, H, Lge, Leg, hbar, n, k1)
    tmp = rho + 0.50 * dt * k1
    call compute_drho_lindblad(tmp, H, Lge, Leg, hbar, n, k2)
    tmp = rho + 0.50 * dt * k2
    call compute_drho_lindblad(tmp, H, Lge, Leg, hbar, n, k3)
    tmp = rho + dt * k3
    call compute_drho_lindblad(tmp, H, Lge, Leg, hbar, n, k4)

    rho_out = rho + dt/6.00 * (k1 + 2*k2 + 2*k3 + k4)

    deallocate(k1, k2, k3, k4, tmp)

end subroutine

subroutine compute_drho_lindblad(rho, H, Lge, Leg, hbar, n, drho)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: hbar
    complex, intent(in) :: rho(n,n), H(n,n), Lge(n,n), Leg(n,n)
    complex, intent(out) :: drho(n,n)

    complex :: iC
    complex, dimension(n,n) :: comm, LrhoLdag, LdagL, anticomm, Ldag

    iC = (0.0d0, 1.0d0)
    comm = matmul(H, rho) - matmul(rho, H)

    ! Dissipator for Leg
    Ldag = transpose(conjg(Leg))
    LrhoLdag = matmul(Leg, matmul(rho, Ldag))
    LdagL = matmul(Ldag, Leg)
    anticomm = matmul(LdagL, rho) + matmul(rho, LdagL)

    drho = -(iC/hbar)*comm + LrhoLdag - 0.50 * anticomm

    ! Dissipator for Lge
    Ldag = transpose(conjg(Lge))
    LrhoLdag = matmul(Lge, matmul(rho, Ldag))
    LdagL = matmul(Ldag, Lge)
    anticomm = matmul(LdagL, rho) + matmul(rho, LdagL)

    drho = drho + LrhoLdag - 0.50 * anticomm

end subroutine

! using Euler's method
subroutine compute_drho(rho, H, L, gamma, hbar, n_states, drho)
    implicit none
    integer :: n, i, n_states
    complex, intent(in)  :: rho(n_states,n_states), H(n_states,n_states), L(n_states,n_states)
    real, intent(in)  :: gamma, hbar
    complex, intent(out) :: drho(n_states,n_states)

    ! some local variables
    complex, allocatable :: comm(:,:), dissip(:,:), LdL(:,:), tmp1(:,:), tmp2(:,:)
    complex :: iC

    ! so everything is of the size n_statesxn_states
    allocate(comm(n_states,n_states), dissip(n_states,n_states), LdL(n_states,n_states), tmp1(n_states,n_states), tmp2(n_states,n_states))
    iC = (0.000, 1.000)

    ! term1: commutator [H,ρ]
    comm = matmul(H, rho) - matmul(rho, H)

    ! L†L
    LdL = matmul((transpose(conjg(L))), L)

    ! dissipator = LρL† − ½(L†L·ρ + ρ·L†L)
    tmp1   = matmul(L, rho)
    tmp2   = matmul(tmp1, (transpose(L)))  ! LρL†
    dissip = tmp2 - 0.5*(matmul(LdL, rho) + matmul(rho, LdL))

    ! full right-hand side
    drho = -(iC/hbar)*comm + gamma * dissip

    deallocate(comm, dissip, LdL, tmp1, tmp2)
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

end module solve_me2