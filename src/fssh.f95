module fssh
    public func

contains

function calc_boltz(E_new, E_old, temp, boltzmann) result(boltz_factor)

    implicit none

    real, intent(in) :: E_new, E_old, temp
    integer, intent(in) :: boltzmann 
    real :: dE, args, boltz_factor, kb

    dE = E_new - E_old

    if (boltzmann .EQ. 0) then
        boltz_factor = 1.0

    else
        kb = 8.617333262 * (10.**(-5)) !units- eV K-1
        if (dE .GT. 0.0000000000000) then
            args = dE/(kb*temp)

            if (args .GT. 50.0) then
                boltz_factor = 0.0
            else
                boltz_factor = (exp(-args))
            end if

        else
            boltz_factor = 1.0

        end if
    end if

return
end function calc_boltz

subroutine fsshpop(den, hamic, n_states, nelect, nsteps, temp, boltzmann, istate, shpop)

    implicit none

    ! general parameters from the main function
    integer, intent(in) :: n_states, nelect, nsteps, boltzmann, istate
    real, intent(in) :: temp
 
    ! hamiltonian matrix required for boltzmann calculation
    real, dimension(nsteps*nelect, n_states, n_states) :: hamic(nsteps*nelect, n_states, n_states)    

    ! the density matrix is required for probability calculation
    real, intent(in) :: den(nsteps*nelect, n_states, n_states)

    ! the array to store probability
    real, dimension(nsteps-1, n_states) :: prob(nsteps-1, n_states)

    ! the array for sh population
    integer, intent(out) :: shpop(nsteps, n_states)

    ! general variables
    integer :: istep, ist, state, fromstate, tostate
    real :: boltz_factor, rn, E_new, E_old

    shpop(:,:) = 0
    shpop(1, istate) = 1

    do ist = 1, (nsteps-1)
        istep = ist*nelect

        do state = 1, n_states

            ! using the straight formula given in SI of pyxaid paper
            prob(ist, state) = (den(istep, state, state) - den((istep+nelect), state, state))/den(istep, state, state)

            ! negative probabilities are equated to 0
            if (prob(ist, state) .LT. 0) then
                prob(ist, state) = 0
            end if

            ! boltzmann from the current state to the other state
            fromstate = state
           
            if (state .EQ. 2) then
                tostate = 1
            else if (state .EQ. 1) then
                tostate = 2          
            end if

            ! the probability is scaled by boltzmann factor
            E_new = real(hamic(istep, tostate, tostate))
            E_old = real(hamic(istep, fromstate, fromstate))

            ! 
            boltz_factor = calc_boltz(E_new, E_old, temp, boltzmann)
            prob(ist, state) = prob(ist, state)*boltz_factor

            ! the probability is compared to a random number  
            if (shpop(ist, state) .EQ. 1) then
                call random_number(rn)
                if (rn .LT. prob(ist, state)) then
                    shpop(ist+1, :) = 1
                    shpop(ist+1, state) = 0
                else
                    shpop(ist+1, :) = shpop(ist, :)
                end if
            end if
        end do
    end do

end subroutine fsshpop

end module fssh
