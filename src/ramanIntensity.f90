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
integer :: n1
  !! Number of integration steps to take without
  !! considering the limit cutoff or smearing
integer :: n2
  !! Hardcoded to be 100000; used to generate `ex1`
integer :: nmode
  !! Number of phonon modes
integer :: nprocs
  !! Number of MPI processes
integer :: tmp_i
real(kind = dp) :: beta
  !! \(\beta = 1/k_{B}T\)
real(kind = dp) :: count1
real(kind = dp) :: elevel
  !! \(E_n\)
real(kind = dp) :: gamma_p
real(kind = dp) :: alpha
real(kind = dp) :: loglimit
real(kind = dp) :: limit
real(kind = dp) :: omega
real(kind = dp) :: step
real(kind = dp) :: dstep
  !! Step to go from 0 to \(\2\pi) in `n2` steps
real(kind = dp) :: t
real(kind = dp) :: temperature
real(kind = dp) :: tmp_r
  !! Temporary variable used in linear 
  !! interpolation of \(t\) exponential
real(kind = dp) :: x
real(kind = dp) :: y

complex(kind = dp) :: expForFj
  !! An exponential form to more quickly
  !! calculate \(F_j\)
complex(kind = dp) :: expT
  !! \(e^{i\omega t}\) used in calculating \(F_j\)
complex(kind = dp) :: Fj
  !! Function \(F_j\) from equation 44
complex(kind = dp) :: tmp_exp
complex(kind = dp) :: zfactor

character(len = 256) :: SjOutputFile
  !! Input file name
character(len = 256) :: dummy

integer,allocatable :: interval(:)

real(kind = dp), allocatable :: count2(:)
real(kind = dp), allocatable :: hbarOmegaBeta(:)
real(kind = dp), allocatable :: domega(:)
  !! \(\delta\omega_{nj} = \omega_{nj} - \omega_j\)
real(kind = dp), allocatable :: elaser(:)
  !! Laser energies \(E_L\)
real(kind = dp), allocatable :: eshift(:)
real(kind = dp), allocatable :: omega_s(:)
real(kind = dp), allocatable :: omega_j(:)
  !! \(\omega_j\)
real(kind = dp), allocatable :: omega_l(:)
real(kind = dp), allocatable :: omega_nj(:)
  !! \(\omega_{nj}\)
real(kind = dp), allocatable :: Sj(:)
  !! \(S_j = \dfrac{\omega_j^2}{2\hbar}\delta q_j^2\)
  !! Taken from input file

complex(kind = dp), allocatable :: ex1(:)
  !! \(e^{i\theta}\) where \(\theta\) goes from 0 to \(2\pi\)
complex(kind = dp), allocatable :: expX(:)
  !! \(e^{i\omega x}\) used in calculating \(F_j\)
complex(kind = dp), allocatable :: expY(:)
  !! \(e^{-i\omega y}\) used in calculating \(F_j\)
complex(kind = dp), allocatable :: FjFractionFactor(:)
  !! Fraction factor in front of cosine terms in \(F_j\)
complex(kind = dp), allocatable :: global_sum(:,:)
complex(kind = dp), allocatable :: s1(:)
complex(kind = dp), allocatable :: s2(:,:)
complex(kind = dp), allocatable :: s3(:,:)
complex(kind = dp), allocatable :: theta(:)
  !! Serves as the argument for the sines in the fraction
  !! factor `FjFractionFactor` and the exponentials in
  !! `zfactor1` 
complex(kind = dp), allocatable :: zfactor1(:)
complex(kind = dp), allocatable :: zfactor2(:)

namelist /ramanInput/ temperature, n1, limit, gamma_p, alpha, elevel, &
                      elaser_num, eshift_num, SjOutputFile


call MPI_INIT(ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, id, ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
  !! * Initialize MPI pool

omega=1.0d14
  !! * Define a number to scale all inputs down by. This is currently 
  !!   set to 2 orders of magnitude larger than the phonon frequency
  !!   scale so that the time step is small enough to get reasonable 
  !!   results.
n2 = 100000
  !! * Define the number of exponentials to pre-calculate in order
  !!   to do the linear interpolation of \(e^(i\omega_j t)\)

if(id == 0) then
  !! * If root process
  !!    * Read temperature, `n1`, `limit`, `gamma_p`, `alpha`, 
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
allocate(Sj(nmode), omega_j(nmode), omega_nj(nmode), hbarOmegaBeta(nmode), FjFractionFactor(nmode), expX(nmode), expY(nmode))
allocate(domega(nmode), theta(nmode), zfactor1(nmode), zfactor2(nmode), ex1(0:n2+1), interval(2), count2(2))
  !! * Allocate space for variables on all processes

dstep = tpi/float(n2)
  !! * Calculate the step size for the exponential pre-calculation

do j = 0, n2+1
  !! * Pre-calculate the exponential terms, \(e^{i\theta}\) where \(\theta\)
  !!   goes from 0 to \(2\pi\), used in the linear interpolation of \(e^{i\omega_j t}\)

   ex1(j)=exp(I*j*dstep)

end do

if(id == 0) then
  !! * If root process
  !!    * Scale the phonon frequencies down by 100
  !!    * For each mode, read index (currently not being used), 
  !!      \(S_j\), \(\omega_j\), and \(\omega_{nj}\)
  !!    * Calculate \(\text{hbarOmegaBeta}=\hbar\omega\beta\omega_j\) for
  !!      each \(\omega_j\)
  !!    * Calculate \(\delta\omega_{nj} = \omega_{nj} - \omega_j\) for all
  !!      \(\omega_j\)/\(\omega_{nj}\) pair
  !!    * Read in the energy shifts and laser energies
  !!      @note 
  !!        Reading the input file assumes the form `&ramanInput ... /`,
  !!        blank line, ESHIFT, list of energy shifts to calculate
  !!        intesities for, blank line, ELASER, laser energie(s) to use
  !!      @endnote

  omega_j(:) = omega_j(:)/100.0d0
  omega_nj(:) = omega_nj(:)/100.0d0     
  

  do imode = 1, nmode
     
     read(11,*) j, Sj(imode), omega_j(imode), omega_nj(imode)

  end do

     hbarOmegaBeta(:) = hbar*omega_j(:)*omega*beta
     
     domega(:) = omega_nj(:)  - omega_j(:)
  
     read(5,*)
     read(5,*) eshift(:)
     read(5,*)
     read(5,*)
     read(5,*) elaser(:)
      ! Make sure to skip the blank lines and headers. 

endif

call MPI_Bcast( n1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
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
call MPI_Bcast( domega, nmode, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror) 
call MPI_Barrier(MPI_COMM_WORLD,ierror)
  !! * Broadcast input variables to other processes

omega_l(:) = (elaser(:) - elevel)*ev/(hbar*omega)
  ! \((E_L-E_n)/\hbar\omega\)
  !! @todo Break this off into `omega_l` and `omega_a` for clarity @endtodo
omega_s(:) = eshift(:)*mev/(hbar*omega)
gamma_p = gamma_p*mev/(hbar*omega)
alpha = alpha*mev/(hbar*omega)
  !! * Scale down `omega_l`, `omega_s`, `gamma_p`, and `alpha`
  !!   to ensure that integration scale is small enough to
  !!   get a reasonable result

step = tpi/float(n1)
loglimit = -log(limit)
count1 = 0.0
  ! Set loop variables

!> * Calculate the total number of integration steps needed and
!>   distribute the intervals between the processes
do j = 0, int(loglimit/gamma_p/step)
   do k = 0, int(loglimit/gamma_p/step-j)
      count1 = count1 + int((loglimit/step - gamma_p*j - gamma_p*k)/alpha)
   enddo
end do

count2(1) = float(id)/float(nprocs)*count1
count2(2) = float(id+1)/float(nprocs)*count1
indexI = 1
count1 = 0.0

! Make sure the boundaries between processes are set properly
do j = 0, int(loglimit/gamma_p/step)

   do k = 0, int(loglimit/gamma_p/step)-j
      count1 = count1 + int((loglimit/step - gamma_p*j - gamma_p*k)/alpha)
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

   if(id==0) then
      write(*,*) id, iX
   endif

   s2 = 0.0d0
   x = (iX + 0.5) * step

   expX(:) = cos(omega_nj(:)*x) + I*sin(omega_nj(:)*x)
 
   do iY = 0, int(loglimit/gamma_p/step) - iX 
      s1 = 0.0d0
      y = (iY + 0.5) * step
      
      expY(:) = cos( omega_nj(:)*y ) - I * sin(omega_nj(:)*y)
      
      theta(:) = hbarOmegaBeta(:) + I*domega(:) * ( x - y ) 
      zfactor1(:) = exp(0.5*theta(:)) / ( exp(theta(:)) - 1 )        
       !! Calculate the first fraction in equation 42
      zfactor2(:) =exp(0.5*hbarOmegaBeta(:)) / ( exp(hbarOmegaBeta(:)) - 1 )
       !! @todo Figure out where `zfactor2` comes from @endtodo
      zfactor = product(zfactor1(:)/zfactor2(:))
      !zfactor = product(exp( 0.5*I*domega(:)*(x-y) ) * ( exp( hbarOmegaBeta(:) ) - 1 ) / ( exp( theta(:) ) - 1 ))

      theta(:) = domega(:)*( x - y ) - I*hbarOmegaBeta(:)
      FjFractionFactor(:) = sin(theta(:))/( 1 - cos(theta(:))  )
       !! * Calculate the fractional factor in front of the cosines in
       !!   \(F_j\)

      interval_t = int((loglimit/step - gamma_p*iX - gamma_p*iY)/alpha)

      do iT = 0, interval_t
         t = (iT + 0.5)*step
         tmp_exp = 0.0d0
 
         do imode = 1, nmode

            tmp_r = -omega_j(imode)*t/tpi
            tmp_r = (tmp_r - floor(tmp_r))*float(n2)
            indexI = floor(tmp_r)
            tmp_r = tmp_r - indexI
            expT = (1.0 - tmp_r)*ex1(indexI) + tmp_r*ex1(indexI+1)
              !! Use a linear interpolation for \(e^{i\omega_j t}\)

            expForFj = -expT*(1 - expX(imode))*(1 - expY(imode)) - (expX(imode) + expY(imode))
              !! * Calculate `expForFj`\( = -e^{-i\omega t}(1-e^{i\omega x})(1-e^{-i\omega y})-(e^{i\omega x} + e^{-i\omega y})\).
              !!   This is used as a trick to be able to calculate \(F_j\) quicker as the expontentials include both the
              !!   sines and cosines needed 
            Fj = Aimag(expForFj) + FjFractionFactor(imode)*Real(2.0d0 + expForFj)           
              !! * Calculate \(F_j = \text{Im}(\)`expForFj`\() + \)`FjFractionFactor`\(\text{Re}(2 + \)`expForFj`\()\)
              !! @todo Add detailed derivation of this in a separate page @endtodo
            tmp_exp = tmp_exp + Fj*Sj(imode)
         enddo

         s1(:) = s1(:) + exp(I*tmp_exp - I*omega_s(:)*t - alpha*abs(t))

      end do

      do ilaserE = 1, elaser_num

        s2(:,ilaserE) = s2(:,ilaserE) + s1(:)*exp(-(I*omega_l(ilaserE) + gamma_p)*y) * zfactor

      enddo

   end do

   do ilaserE = 1, elaser_num

     s3(:,ilaserE) = s3(:,ilaserE) + s2(:,ilaserE)*exp(-(-I*omega_l(ilaserE) + gamma_p)*x)
   
   enddo
enddo


call MPI_Barrier(MPI_COMM_WORLD,ierror) 
call MPI_Reduce( s3, global_sum, eshift_num*elaser_num, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierror)



if(id == 0) then
  write(12,*)"calculation finalized"

  do ilaserE = 1, elaser_num
    write(12,*) "Laser energy: ", elaser(ilaserE)
    
    do ishiftE = 1, eshift_num 
      write(12,*) eshift(iShiftE), eshift(iShiftE)*mevtocm, step**3*Real(global_sum(iShiftE,ilaserE))*2.0
    enddo

  enddo
endif
call MPI_FINALIZE(ierror)

end program ramanIntensity

         
