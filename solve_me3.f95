module solve_me3
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
subroutine me(hamic, nacic, boltz, temperature, n_states, nelect, dtelec, n_files, nsteps, istate, gamma, den)

    implicit none

    ! only the part of the hamiltonian corresponding to the initial condition is stored
    real, intent(in) :: hamic(nsteps*nelect, n_states, n_states)
    real, intent(in) :: nacic(nsteps*nelect, n_states, n_states)

    ! the general input for this function
    integer, intent(in) :: n_states, nelect, n_files, istate, nsteps, boltz
    real, intent(in) :: dtelec
    real, intent(in) :: gamma  ! decay rate

    ! defining the density matrix
    real, intent(out) :: den(nsteps*nelect, n_states, n_states)

    ! dummy variables that are relavant only to this function
    real :: factor_re, factor_im, boltz_factor, temperature, rate_ge, rate_eg
    integer :: state, istep

    ! — locals —
    integer :: total_steps, i, j
    real, parameter :: hbar = 0.6582119510   ! eV·fs
    real, parameter :: kb = 0.00008617333   ! eV/K

    complex, dimension(n_states,n_states) :: rho(n_states,n_states)
    complex, dimension(n_states,n_states) :: rho_next(n_states,n_states)
    real, dimension(n_states,n_states) :: L_ge(n_states,n_states)
    real, dimension(n_states,n_states) :: L_eg(n_states,n_states)

    complex, allocatable :: drho(:,:), H(:,:), nac(:,:), Ham(:,:)

    ! total electronic steps
    total_steps = nsteps * nelect

    ! 1) initialize ρ(0) = |istate><istate|
    rho = (0.00, 0.00)
    rho(istate, istate) = (1.000, 0.000)
    den(1,:,:) = real(rho)
    
    ! 2) time‐propagation loop (forward Euler)
    do istep = 1, total_steps-1

        ! defining a constant L_ge for check
        L_ge(1,1)=0.0
        L_ge(1,2)=1.0
        L_ge(2,1)=0.0
        L_ge(2,2)=0.0

        ! defining a constant L_eg for check
        L_eg(1,1)=0.0
        L_eg(1,2)=0.0
        L_eg(2,1)=1.0
        L_eg(2,2)=0.0

        rate_ge = 0.0
        rate_eg = 0.0

        ! allocate complex matrices
        allocate(drho(n_states,n_states))
        allocate(H(n_states,n_states))
        allocate(nac(n_states,n_states))


        ! 2a) build H and L at this step
        do i = 1, n_states
            do j = 1, n_states
                H(i,j) = CMPLX( hamic (istep, i, j), 0.00 )
                nac(i,j) = CMPLX( 0.00, nacic(istep, i, j))

                ! ! Boltzmann factor
                ! if (boltz .EQ. 1) then
                !     if (real(H(j,j)) > real(H(i,i))) then
                !         boltz_factor = exp(-kb/temperature) * (H(j,j) - H(i,i))
                !         L_ge(i,j) = L_ge(i,j)*sqrt(boltz_factor)
                !         L_eg(i,j) = L_eg(i,j)*sqrt(boltz_factor)
                !     end if
                ! end if
            end do
        end do

        rate_eg = (2*3.14/hbar)*nacic(istep, 2, 1)*nacic(istep, 2, 1)*1

        if (real(H(2,2)) > real(H(1,1))) then
            rate_ge = (2*3.14/hbar)*nacic(istep, 1, 2)*nacic(istep, 1, 2)*(exp(-kb/temperature) * (H(2,2) - H(1,1)))
        end if 

        ! 2b) compute drho = −(i/ħ)[H,ρ] + γ (LρL† − ½{L†L,ρ})
        ! call compute_drho(rho, H, L, gamma, hbar, n_states, drho)
        call rk4_lindblad(rho, H, L_ge, rate_ge, L_eg, rate_eg, hbar, dtelec, n_states, rho_next)

        ! 2c) Euler update
        ! rho_next = rho + dtelec * drho

        ! 2d) renormalise Tr(ρ) → 1
        call renormalise(n_states, rho_next)

        ! 2e) save for output and advance
        den(istep+1,:,:) = real(rho_next)
        rho = rho_next
        
        ! deallocate scratch
        deallocate(drho, H, nac)

    end do

end subroutine me

subroutine rk4_lindblad(rho, H, L_ge, rate_ge, L_eg, rate_eg, hbar, dt, n, rho_out)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: dt, hbar, rate_ge, rate_eg
    complex, intent(in) :: rho(n,n), H(n,n)
    real, intent(in) :: L_ge(n,n), L_eg(n,n)
    complex, intent(out) :: rho_out(n,n)

    ! g to e
    complex, allocatable :: k1(:,:), k2(:,:), k3(:,:), k4(:,:), tmp(:,:)

    allocate(k1(n,n), k2(n,n), k3(n,n), k4(n,n), tmp(n,n))

    call compute_drho_lindblad(rho, H, L_ge, rate_ge, L_eg, rate_eg, hbar, n, k1)
    tmp = rho + 0.50 * dt * k1
    call compute_drho_lindblad(tmp, H, L_ge, rate_ge, L_eg, rate_eg, hbar, n, k2)
    tmp = rho + 0.50 * dt * k2
    call compute_drho_lindblad(tmp, H, L_ge, rate_ge, L_eg, rate_eg, hbar, n, k3)
    tmp = rho + dt * k3
    call compute_drho_lindblad(tmp, H, L_ge, rate_ge, L_eg, rate_eg, hbar, n, k4)

    rho_out = rho + dt/6.00 * (k1 + 2*k2 + 2*k3 + k4)
    
    deallocate(k1, k2, k3, k4, tmp)

end subroutine

subroutine compute_drho_lindblad(rho, H, L_ge, rate_ge, L_eg, rate_eg, hbar, n, drho)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: hbar, rate_ge, rate_eg
    complex, intent(in) :: rho(n,n), H(n,n)
    real, intent(in) :: L_ge(n,n), L_eg(n,n)
    complex, intent(out) :: drho(n,n)

    complex :: iC
    complex, dimension(n,n) :: comm, LrhoLdag_eg, LdagL_eg, anticomm_eg, Ldag_eg, LrhoLdag_ge, LdagL_ge, anticomm_ge, Ldag_ge

    iC = (0.0d0, 1.0d0)
    comm = matmul(H, rho) - matmul(rho, H)

    ! Dissipator_eg
    Ldag_eg = transpose((L_eg))
    LrhoLdag_eg = matmul(L_eg, matmul(rho, Ldag_eg))
    LdagL_eg = matmul(Ldag_eg, L_eg)
    anticomm_eg = matmul(LdagL_eg, rho) + matmul(rho, LdagL_eg)

    ! Dissipator_ge
    Ldag_ge = transpose((L_ge))
    LrhoLdag_ge = matmul(L_ge, matmul(rho, Ldag_ge))
    LdagL_ge = matmul(Ldag_ge, L_ge)
    anticomm_ge = matmul(LdagL_ge, rho) + matmul(rho, LdagL_ge)

    drho = -(iC/hbar)*comm + rate_eg*(LrhoLdag_eg - 0.50 * anticomm_eg) + rate_ge*(LrhoLdag_ge - 0.50 * anticomm_ge)

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
    LdL = matmul((transpose(L)), L)

    ! dissipator = LρL† − ½(L†L·ρ + ρ·L†L)
    tmp1   = matmul(L, rho)
    tmp2   = matmul(tmp1, transpose(L))  ! LρL†
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

end module solve_me3