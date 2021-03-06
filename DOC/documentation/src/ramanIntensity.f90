module para
  integer, parameter :: dp = selected_real_kind(15, 307)

  real(kind = dp), parameter :: ev = 1.6021766d-19
    !! Conversion factor from eV to J
  real(kind = dp), parameter :: hbar = 1.0545718d-34
    !! \(\bar{h}\)
  real(kind = dp), parameter :: kB = 1.38064852d-23
    !! Boltzmann constant
  real(kind = dp), parameter :: mev = 1.6021766d-22
    !! Conversion factor from meV to J
  real(kind = dp), parameter :: mevtocm = 8.0655438354 
    !! Conversion factor from meV to cm\(^{-1}\)
  real(kind = dp), parameter :: pi  =  3.1415926535897932_dp
    !! \(\pi\)
  real(kind = dp), parameter :: tpi  =  2.0*3.1415926535897932_dp
    !! \(2\pi\)
  complex(kind = dp), parameter :: I = cmplx(0.0d0, 1.0d0, dp)  
    !! \(i\)
end module

program ramanIntensity
  !! The program calculates the normalized Raman intensity
  !! at given energy shifts for given laser energies.
  !!
  !! Input variables:
  !!
  !!  * `nIntSteps` -- number of integration steps
  !!  * `limit` -- limit to truncate the integration at
  !!  * `gamma_p` -- the lifetime \(\gamma\) of the electronic 
  !!     state \(|n\rangle\)
  !!  * `alpha` -- smearing parameter
  !!  * `elevel` -- energy \(E_a\) of the intermediate state
  !!  * `elaser(elaser_num)` -- laser energies
  !!  * `eshift(eshift_num)` -- energy shifts to calculate 
  !!    intensities at
  !!  * `Sj` -- \(S_j = \dfrac{\omega_j^2}{2\hbar}\delta q_j^2\)
  !!  * `omega_j` -- phonon frequencies \(\omega_j\) in the 
  !!     ground state \(|0\rangle\)
  !!  * `omega_{nj}` -- phonon frequencies \(\omega_{nj}\) in
  !!    the electronic states \(|n\rangle\)
  !!
  !! Output variables:
  !!
  !!  * `elaser` -- same laser energies from input
  !!  * `eshift` (meV) -- same energy shifts from input
  !!  * `eshift` (cm\(^{-1}\)) -- converted energy shifts
  !!  * Normalized Raman intensity for each laser energy and 
  !!    energy shift
  !!
  !! <h2>Walkthrough</h2>
  !!

use para
use mpi

implicit none

integer :: elaser_num
  !! Number of laser energies used
integer :: eshift_num
  !! Number of energy shifts to calculate
integer :: ierror
  !! Error for MPI
integer :: id
  !! MPI process id
integer :: ilaserE
  !! Loop index over laser energies
integer :: imode
  !! Loop index over phonon modes
integer :: indexI
  !! Calculated index for constructing
  !! integration intervals and linear
  !! extrapolation of \(e^{i\omega_j t}\)
integer :: interval_t
  !! Interval for time integration
integer :: ishiftE
  !! Loop index over energy shifts
integer :: iX, iY, iT
  !! Loop indices over \(x\), \(y\), and \(t\)
integer :: j, k
  !! Misc loop indicies
integer :: nExpSteps
  !! Hardcoded to be 100000; used to generate `ex1`
integer :: nIntSteps
  !! Number of integration steps to take without
  !! considering the limit cutoff or smearing
integer :: nmode
  !! Number of phonon modes
integer :: nprocs
  !! Number of MPI processes

real(kind = dp) :: alpha
  !! Smearing factor
real(kind = dp) :: beta
  !! \(\beta = 1/k_{B}T\)
real(kind = dp) :: count1
  !! Count used to set interval limits
  !! for each process
real(kind = dp) :: elevel
  !! Energy \(E_a\) of the intermediate state
real(kind = dp) :: gamma_p
  !! The lifetime of the electronic state 
  !! \(|n\rangle\). Input as \(\gamma\) and
  !! changed to \(\gamma/\hbar\)
real(kind = dp) :: limit
  !! Limit to truncate integration at since
  !! the integrand is exponentially decaying
real(kind = dp) :: loglimit
  !! \(\log(\text{limit})\)
real(kind = dp) :: omega_a
  !! \(E_a/\hbar\)
real(kind = dp) :: intStep
  !! Step size for the integration
real(kind = dp) :: expStep
  !! Step to go from 0 to \(2\pi\) in `nExpSteps`
real(kind = dp) :: scalingFactor
  !! Factor to scale down inputs to ensure that
  !! integration scale is small enough to give 
  !! reasonable results
real(kind = dp) :: t
  !! \(t\)
real(kind = dp) :: temperature
  !! The temperature
real(kind = dp) :: tmp_r
  !! Temporary variable used in linear 
  !! interpolation of \(t\) exponential
real(kind = dp) :: x
  !! \(x\)
real(kind = dp) :: y
  !! \(y\)

complex(kind = dp) :: zfactor
  !! The pre-factor in the product in equation 42
  !! and something else

character(len = 256) :: SjOutputFile
  !! Ouput file name from the \(S_j\) calculation

integer,allocatable :: interval(:)
  !! Defines the integration interval for each process

real(kind = dp), allocatable :: count2(:)
  !! Count used to set interval limits
  !! for each process
real(kind = dp), allocatable :: hbarOmegaBeta(:)
  !! \(\hbar\omega_j\beta\)
real(kind = dp), allocatable :: domega(:)
  !! \(\delta\omega_{nj} = \omega_{nj} - \omega_j\)
real(kind = dp), allocatable :: elaser(:)
  !! Laser energies \(E_L\)
real(kind = dp), allocatable :: eshift(:)
  !! Energy shifts to calculate the intensity for
real(kind = dp), allocatable :: omega_j(:)
  !! \(\omega_j\)
real(kind = dp), allocatable :: omega_l(:)
  !! \(E_L/\hbar\)
real(kind = dp), allocatable :: omega_nj(:)
  !! \(\omega_{nj}\)
real(kind = dp), allocatable :: omega_s(:)
  !! \(E_s/\hbar\)
real(kind = dp), allocatable :: Sj(:)
  !! \(S_j = \dfrac{\omega_j^2}{2\hbar}\delta q_j^2\)
  !! from output file of previous \(S_j\) calculation

complex(kind = dp), allocatable :: ex1(:)
  !! \(e^{i\theta}\) where \(\theta\) goes from 0 to \(2\pi\)
  !! used in linear interpolation of \(e^{i\omega_j t}\)
complex(kind = dp), allocatable :: expForFj(:)
  !! An exponential form to more quickly
  !! calculate \(F_j\)
complex(kind = dp), allocatable :: expT(:)
  !! \(e^{i\omega t}\) used in calculating \(F_j\)
complex(kind = dp), allocatable :: expX(:)
  !! \(e^{i\omega x}\) used in calculating \(F_j\)
complex(kind = dp), allocatable :: expY(:)
  !! \(e^{-i\omega y}\) used in calculating \(F_j\)
complex(kind = dp), allocatable :: Fj(:)
  !! Function \(F_j\) from equation 44
complex(kind = dp), allocatable :: FjFractionFactor(:)
  !! Fraction factor in front of cosine terms in \(F_j\)
complex(kind = dp), allocatable :: global_sum(:,:)
  !! Total integral from all processes
complex(kind = dp), allocatable :: s1(:)
  !! Integral over \(t\) for a single process
complex(kind = dp), allocatable :: s2(:,:)
  !! Integral over \(y\) and \(t\) for a single
  !! process
complex(kind = dp), allocatable :: s3(:,:)
  !! Integral over \(x\), \(y\), and \(t\) for
  !! a single process
complex(kind = dp), allocatable :: theta(:)
  !! Serves as the argument for the sines in the fraction
  !! factor `FjFractionFactor` and the exponentials in
  !! `TFractionFactor` 
complex(kind = dp), allocatable :: TFractionFactor(:)
  !! Fraction factor in equation 42
complex(kind = dp), allocatable :: partitionFunction(:)
  !! Partition function \(Z = \dfrac{e^{\frac{1}{2}\beta\hbar\omega_j}}{e^{\beta\hbar\omega_j} - 1}\)

! Define a namelist to read in all of the input variables from the input file
namelist /ramanInput/ temperature, nIntSteps, limit, gamma_p, alpha, elevel, &
                      elaser_num, eshift_num, SjOutputFile


call MPI_INIT(ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, id, ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
  !! * Initialize MPI pool

scalingFactor = 1.0d14
  !! * Define a number to scale all inputs down by. This is currently 
  !!   set to 2 orders of magnitude larger than the phonon frequency
  !!   scale so that the time step is small enough to get reasonable 
  !!   results.
nExpSteps = 100000
  !! * Define the number of exponentials to pre-calculate in order
  !!   to do the linear interpolation of \(e^{i\omega_j t}\)

if(id == 0) then
  !! * If root process
  !!    * Read temperature, `nIntSteps`, `limit`, `gamma_p`, `alpha`, 
  !!      `elevel`, `elaser_num`, `eshift_num`, and `SjOutputFile`
  !!    * Open `SjOutputFile` to read the number of phonon modes
  !!      which is needed to allocate variables
  !!    * Calculate \(\beta = 1/k_{B}T\)

  read(5, ramanInput) 
  read(5,*)
    ! Read the `ramanInput` namelist then skip the next blank line

  open(11, file=trim(SjOutputFile), Action='read', status='old')
  read(11,*)
  read(11,*) nmode

  beta = 1/(kB*temperature)
endif

call MPI_BCast( nmode, 1, MPI_Integer, 0, MPI_COMM_WORLD, ierror)
call MPI_BCast( eshift_num, 1, MPI_Integer, 0, MPI_COMM_WORLD, ierror)
call MPI_BCast( elaser_num, 1, MPI_Integer, 0, MPI_COMM_WORLD, ierror)
  !! * Broadcast array sizes to other modes for allocation

allocate(eshift(eshift_num), elaser(elaser_num), omega_s(eshift_num), omega_l(elaser_num))
allocate(s1(eshift_num), s2(eshift_num,elaser_num), s3(eshift_num,elaser_num), global_sum(eshift_num,elaser_num))
allocate(Sj(nmode), omega_j(nmode), omega_nj(nmode), hbarOmegaBeta(nmode), FjFractionFactor(nmode))
allocate(Fj(nmode), expForFj(nmode), expT(nmode), expX(nmode), expY(nmode))
allocate(domega(nmode), theta(nmode), TFractionFactor(nmode), partitionFunction(nmode), ex1(0:nExpSteps+1), interval(2), count2(2))

expStep = tpi/float(nExpSteps)
  !! * Calculate the step size for the exponential pre-calculation

do j = 0, nExpSteps+1
  !! * Pre-calculate the exponential terms, \(e^{i\theta}\) where \(\theta\)
  !!   goes from 0 to \(2\pi\), used in the linear interpolation of \(e^{i\omega_j t}\)

   ex1(j)=exp(I*j*expStep)

end do

if(id == 0) then
  !! * If root process
  !!    * Scale the phonon frequencies down by 100
  !!    * For each mode, read index (currently not being used), 
  !!      \(S_j\), \(\omega_j\), and \(\omega_{nj}\)
  !!    * Calculate \(\text{hbarOmegaBeta}=\hbar\omega\beta\omega_j\) for
  !!      each \(\omega_j\)
  !!    * Calculate the partition function for each \(\omega_j\)
  !!    * Calculate \(\delta\omega_{nj} = \omega_{nj} - \omega_j\) for all
  !!      \(\omega_j\)/\(\omega_{nj}\) pair
  !!    * Read in the energy shifts and laser energies
  !!      @note 
  !!        Reading the input file assumes the form `&ramanInput ... /`,
  !!        blank line, `ESHIFT`, list of energy shifts to calculate
  !!        intesities for, blank line, `ELASER`, laser energie(s) to use
  !!      @endnote

  omega_j(:) = omega_j(:)/100.0d0
  omega_nj(:) = omega_nj(:)/100.0d0     
  

  do imode = 1, nmode
     
     read(11,*) j, Sj(imode), omega_j(imode), omega_nj(imode)

  end do

     hbarOmegaBeta(:) = (hbar*scalingFactor)*omega_j(:)*beta
     partitionFunction(:) =exp(0.5*hbarOmegaBeta(:)) / ( exp(hbarOmegaBeta(:)) - 1 )
     
     domega(:) = omega_nj(:)  - omega_j(:)
  
     read(5,*)
     read(5,*) eshift(:)
     read(5,*)
     read(5,*)
     read(5,*) elaser(:)
      ! Make sure to skip the blank lines and headers. 

endif

call MPI_Bcast( nIntSteps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( limit, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( gamma_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( elaser, elaser_num, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( elevel, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( eshift, eshift_num, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( Sj, nmode, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( omega_j, nmode, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( omega_nj, nmode, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( hbarOmegaBeta, nmode, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror) 
call MPI_Bcast( partitionFunction, nmode, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror) 
call MPI_Bcast( domega, nmode, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror) 
call MPI_Barrier(MPI_COMM_WORLD,ierror)
  !! * Broadcast input variables to other processes

omega_l(:) = (elaser(:)*ev)/(hbar*scalingFactor)
omega_a = (elevel*ev)/(hbar*scalingFactor)
omega_s(:) = (eshift(:)*mev)/(hbar*scalingFactor)
gamma_p = (gamma_p*mev)/(hbar*scalingFactor)
alpha = (alpha*mev)/(hbar*scalingFactor)
  !! * Convert the input energies to J, then define frequencies by dividing the energies 
  !!   by \(\hbar\) that has been scaled to ensure the results to ensure that integration 
  !!   scale is small enough to get a reasonable result
  !!   @note 
  !!      \(\alpha\) and \(\gamma\) are also divided by \(\hbar\) here so 
  !!      that it doesn't have to be done in the final exponential.
  !!   @endnote

if(id /= 0) then
  deallocate(eshift, elaser)
endif

intStep = tpi/float(nIntSteps)
loglimit = -log(limit)
count1 = 0.0
  ! Set loop variables

!> * Calculate the total number of integration steps needed and
!>   distribute the intervals between the processes
do j = 0, int(loglimit/gamma_p/intStep)
   do k = 0, int(loglimit/gamma_p/intStep-j)
      count1 = count1 + int((loglimit/intStep - gamma_p*j - gamma_p*k)/alpha)
   enddo
end do

count2(1) = float(id)/float(nprocs)*count1
count2(2) = float(id+1)/float(nprocs)*count1
indexI = 1
count1 = 0.0

do j = 0, int(loglimit/gamma_p/intStep)

   do k = 0, int(loglimit/gamma_p/intStep)-j
      count1 = count1 + int((loglimit/intStep - gamma_p*j - gamma_p*k)/alpha)
   enddo

   if(count1 >= count2(indexI)) then

      interval(indexI) = j
      indexI = indexI + 1
   
      if(indexI == 3 ) then
         exit
      endif

   endif
end do

interval(2) = interval(2) - 1


call MPI_Barrier(MPI_COMM_WORLD,ierror)
  !! * Make sure that all processes get here before
  !! moving forward

s3 = 0.0d0
do iX = interval(1), interval(2)
  !! * Begin integration over \(x\)

   if(id == 0) then
      !! * If root process, output your progress

      write(*,*) "The root process is on step ", iX, " of ", interval(2)

   endif

   s2 = 0.0d0

   x = (iX + 0.5) * intStep
    !! * Define \(x\) for the current step

   expX(:) = cos(omega_nj(:)*x) + I*sin(omega_nj(:)*x)
    !! * Calculate \(e^{i\omega_{nj}x}\) or 
    !!   \(e^{i\omega_j x^{\prime}}\) for the current \(x\)
 
   do iY = 0, int(loglimit/gamma_p/intStep) - iX 
      !! * Begin the integration over \(y\)

      s1 = 0.0d0

      y = (iY + 0.5) * intStep
        !! * Define \(y\) for the current step
      
      expY(:) = cos( omega_nj(:)*y ) - I * sin(omega_nj(:)*y)
        !! * Calculate \(e^{-i\omega_{nj}y}\) or 
        !!   \(e^{-i\omega_j y^{\prime}}\) for the current \(y\)
      

      theta(:) = hbarOmegaBeta(:) + I*domega(:) * ( x - y ) 
      TFractionFactor(:) = exp(0.5*theta(:)) / ( exp(theta(:)) - 1 )        
      zfactor = product(TFractionFactor(:)/partitionFunction(:))
        !! * Calculate `zfactor` which is
        !!   \(\dfrac{1}{Z}\dfrac{e^{\frac{1}{2}(i\delta\omega_{nj}(x - y)+\beta\hbar\omega_j)}}{e^{(i\delta\omega_{nj}(x - y)+\beta\hbar\omega_j)} - 1}\)


      theta(:) = domega(:)*( x - y ) - I*hbarOmegaBeta(:)
      FjFractionFactor(:) = sin(theta(:))/( 1 - cos(theta(:))  )
       !! * Calculate the fractional factor in front of the cosines in
       !!   \(F_j\) (equation 44)

      interval_t = int((loglimit/intStep - gamma_p*iX - gamma_p*iY)/alpha)

      do iT = 0, interval_t
         !! * Begin the integration over \(t\)

         t = (iT + 0.5)*intStep
          !! * Define \(t\) for the current step
 
         do imode = 1, nmode
            !! * Use a linear interpolation for \(e^{i\omega_j t}\)
            !!   to get a slightly more accurate value
            !!   @note 
            !!     This process needs to stay in the innermost loop; it cannot
            !!     be moved out of the loop because `tmp_r` represents a random
            !!     number that will change slightly each loop.
            !!   @endnote

            tmp_r = -omega_j(imode)*t/tpi
            tmp_r = (tmp_r - floor(tmp_r))*float(nExpSteps)
            indexI = floor(tmp_r)
            tmp_r = tmp_r - indexI
            expT(imode) = (1.0 - tmp_r)*ex1(indexI) + tmp_r*ex1(indexI+1)

         enddo

         expForFj(:) = -expT(:)*(1 - expX(:))*(1 - expY(:)) - (expX(:) + expY(:))
          !! * Calculate `expForFj`\( = -e^{-i\omega t}(1-e^{i\omega x})(1-e^{-i\omega y})-(e^{i\omega x} + e^{-i\omega y})\).
          !!   This is used as a trick to be able to calculate \(F_j\) quicker as the expontentials include both the
          !!   sines and cosines needed 

         Fj(:) = Aimag(expForFj(:)) + FjFractionFactor(:)*Real(2.0d0 + expForFj(:))           
          !! * Calculate \(F_j = \text{Im}(\)`expForFj`\() + \)`FjFractionFactor`\(\text{Re}(2 + \)`expForFj`\()\)
          !! @todo Add detailed derivation of this in a separate page @endtodo

         s1(:) = s1(:) + exp(I*sum(Fj(:)*Sj(:)) - I*omega_s(:)*t - alpha*abs(t))
          !! * Increment the innermost sum that will end up being
          !!   \(\sum\left(\prod e^{iS_jF_j}\right)e^{-\frac{i}{\hbar}E_s t}e^{-\alpha/\hbar|t|}\)
          !!   which is the portion of the integrand in equation 30 that depends on \(t\)
          !!   @note
          !!      The \(e^{-\alpha/\hbar |t|}\) term adds smearing. Even though the equation 
          !!      only has `alpha`, the division by \(\hbar\) is implicit because we already 
          !!      did the division when scaling the input variables.
          !!   @endnote

      end do

      do ilaserE = 1, elaser_num
        !! * For each laser energy input, increment the middle sum that will end up being
        !!   \(\sum\)`zfactor`\(e^{\frac{-\gamma-iE_L+iE_a}{\hbar}y}\)`s1`
        !!   which is the portion of equation 30 that depends on \(y\) multiplied by the 
        !!   inner \(t\) sum
        !!   @note
        !!      Here the fraction factor in the product in equation 42 has been pulled 
        !!      out of the integration over time, but the \(e^{iS_jF_j(x,y,t)}\)
        !!      term must stay inside the \(t\) integral.
        !!   @endnote
        !!   @note
        !!      We are technically doing multiple integrals at once for each of the 
        !!      different laser energies \(E_L\). Because of this, `s2` and `s3` need
        !!      to be two-dimensional, but the laser energy doesn't enter the integral
        !!      over time, so `s1` can be one-dimensional.
        !!   @endnote


        s2(:,ilaserE) = s2(:,ilaserE) + s1(:)*exp(-(I*omega_l(ilaserE) - I*omega_a + gamma_p)*y) * zfactor

      enddo

   end do

   do ilaserE = 1, elaser_num
    !! * For each laser energy input, increment the outermost sum that will end up being
    !!   \(\sum e^{\frac{-\gamma+iE_L-iE_a}{\hbar}x}\)
    !!   which is the portion of equation 30 that depends on \(x\) multiplied by the 
    !!   inner sums over \(y\) and \(t\)

     s3(:,ilaserE) = s3(:,ilaserE) + s2(:,ilaserE)*exp(-(-I*omega_l(ilaserE) + I*omega_a + gamma_p)*x)
   
   enddo
enddo

deallocate(omega_j, omega_nj, omega_s, omega_l)
deallocate(Sj, s1, s2, hbarOmegaBeta, FjFractionFactor)
deallocate(Fj, expForFj, expX, expY, expT, domega)
deallocate(theta, TFractionFactor, partitionFunction, ex1, interval, count2)


call MPI_Barrier(MPI_COMM_WORLD,ierror) 
call MPI_Reduce( s3, global_sum, eshift_num*elaser_num, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
  !! * Once all sums are complete on each process, combine
  !!   the results into `global_sum`

deallocate(s3)

if(id == 0) then
  !! * If root process, output the calculated
  !!   intensity for each laser energy and energy
  !!   shift

  write(*,*)"calculation finalized"

  do ilaserE = 1, elaser_num
    write(*,*) "Laser energy: ", elaser(ilaserE)
    
    do ishiftE = 1, eshift_num 
      write(*,*) eshift(iShiftE), eshift(iShiftE)*mevtocm, intStep**3*Real(global_sum(iShiftE,ilaserE))*2.0
        ! When outputting the global sum, multiply by \(\Delta x\Delta y\Delta t\) to make it 
        ! equivalent to the integral. Also multiply by 2 for some reason?
        !! @note 
        !!  The final intensity is the `global_sum` multiplied by 2*\(\Delta x\Delta y\Delta t\).
        !!  The step sizes convert the sum into an integral and the 2 accounts for the fact that
        !!  we only integrated over the positive half of the total time integral. We do not multiply
        !!  by \(\dfrac{\kappa^{\prime}}{2\pi} = \dfrac{|\beta_{1n}\beta_{2n}|^2}{\hbar\sum_ie^{-\beta\Theta_i}}\)
        !!  because there are many unknown parameters and we only care about the normalized signal
        !!  intensity.
        !! @endnote
    enddo

  enddo

  deallocate(eshift, elaser)

endif


deallocate(global_sum)

call MPI_FINALIZE(ierror)

end program ramanIntensity

         
