var tipuesearch = {"pages":[{"text":"Raman Intensities Installation Dependencies Install all of the required dependencies (commands given for Unix environment): sudo apt-get update\nsudo apt install gfortran\nsudo apt install libmpich-dev\nsudo apt install libopenmpi-dev\nsudo apt-get install python-pip python-dev build-essential\nsudo pip install ford\nsudo apt-get install graphviz\nsudo apt-get install git Compile Source Clone the repository using git clone https://github.com/laurarnichols/raman.git Change into the raman directory Compile using [ADD DIRECTIONS ONCE COMPILE] Documentation This documentation is generated by FORD . Note This documentation is in progress. Please report any issues you find. How to Run I don't know how to run yet. Input and Output","tags":"home","loc":"index.html","title":" Raman Intensities "},{"text":"Contents Programs ramanIntensity Modules para Source Code ramanIntensity.f90 Source Code module para integer , parameter :: dp = selected_real_kind ( 15 , 307 ) real ( kind = dp ), parameter :: ev = 1.6021766d-19 !! Conversion factor from eV to J real ( kind = dp ), parameter :: hbar = 1.0545718d-34 !! \\bar{h} real ( kind = dp ), parameter :: kB = 1.38064852d-23 !! Boltzmann constant real ( kind = dp ), parameter :: mev = 1.6021766d-22 !! Conversion factor from meV to J real ( kind = dp ), parameter :: mevtocm = 8.0655438354 !! Conversion factor from meV to cm&#94;{-1} real ( kind = dp ), parameter :: pi = 3.1415926535897932_dp !! \\pi real ( kind = dp ), parameter :: tpi = 2.0 * 3.1415926535897932_dp !! 2\\pi complex ( kind = dp ), parameter :: I = cmplx ( 0.0d0 , 1.0d0 , dp ) !! i end module program ramanIntensity use para use mpi implicit none integer :: elaser_num !! Number of laser energies used integer :: eshift_num !! Number of energy shifts to calculate integer :: ierror !! Error for MPI integer :: id !! MPI process id integer :: ilaserE !! Loop index over laser energies integer :: imode !! Loop index over phonon modes integer :: index1 integer :: interval_a integer :: interval_t integer :: j integer :: k integer :: l integer :: max_k integer :: max_l integer :: n1 integer :: n2 !! Hardcoded to be 100000; used to generate `ex1` integer :: nmode !! Number of phonon modes integer :: nprocs !! Number of MPI processes integer :: tmp_i real ( kind = dp ) :: beta !! \\beta = 1/k_{B}T real ( kind = dp ) :: count1 real ( kind = dp ) :: elevel !! E_n real ( kind = dp ) :: gamma_p real ( kind = dp ) :: alpha real ( kind = dp ) :: loglimit real ( kind = dp ) :: limit real ( kind = dp ) :: omega real ( kind = dp ) :: step1 real ( kind = dp ) :: step2 !! Step to go from 0 to \\2\\pi) in `n2` steps real ( kind = dp ) :: t real ( kind = dp ) :: temperature real ( kind = dp ) :: tmp_r1 real ( kind = dp ) :: tmp_r2 real ( kind = dp ) :: x real ( kind = dp ) :: y complex ( kind = dp ) :: expForFj !! An exponential form to more quickly !! calculate F_j complex ( kind = dp ) :: expT !! e&#94;{i\\omega t} used in calculating F_j complex ( kind = dp ) :: Fj !! Function F_j from equation 44 complex ( kind = dp ) :: tmp_exp complex ( kind = dp ) :: zfactor character ( len = 256 ) :: Inputfile !! Input file name character ( len = 256 ) :: dummy integer , allocatable :: interval (:) real ( kind = dp ), allocatable :: count2 (:) real ( kind = dp ), allocatable :: hbarOmegaBeta (:) real ( kind = dp ), allocatable :: domega (:) !! \\delta\\omega_{nj} = \\omega_{nj} - \\omega_j real ( kind = dp ), allocatable :: elaser (:) !! Laser energies E_L real ( kind = dp ), allocatable :: eshift (:) real ( kind = dp ), allocatable :: omega2 (:) real ( kind = dp ), allocatable :: omega_j (:) !! \\omega_j real ( kind = dp ), allocatable :: omega_l (:) real ( kind = dp ), allocatable :: omega_nj (:) !! \\omega_{nj} real ( kind = dp ), allocatable :: Sj (:) !! S_j = \\dfrac{\\omega_j&#94;2}{2\\hbar}\\delta q_j&#94;2 !! Taken from input file complex ( kind = dp ), allocatable :: ex1 (:) !! e&#94;{i\\theta} where \\theta goes from 0 to 2\\pi complex ( kind = dp ), allocatable :: expX (:) !! e&#94;{i\\omega x} used in calculating F_j complex ( kind = dp ), allocatable :: expY (:) !! e&#94;{-i\\omega y} used in calculating F_j complex ( kind = dp ), allocatable :: FjFractionFactor (:) !! Fraction factor in front of cosine terms in F_j complex ( kind = dp ), allocatable :: global_sum (:,:) complex ( kind = dp ), allocatable :: s1 (:) complex ( kind = dp ), allocatable :: s2 (:,:) complex ( kind = dp ), allocatable :: s3 (:,:) complex ( kind = dp ), allocatable :: theta (:) !! Serves as the argument for the sines in the fraction !! factor `FjFractionFactor` and the exponentials in !! `zfactor1` complex ( kind = dp ), allocatable :: zfactor1 (:) complex ( kind = dp ), allocatable :: zfactor2 (:) call MPI_INIT ( ierror ) call MPI_COMM_RANK ( MPI_COMM_WORLD , id , ierror ) call MPI_COMM_SIZE ( MPI_COMM_WORLD , nprocs , ierror ) !! * Initialize MPI pool Inputfile = 'Sj.out' !! * Define input file name !! @todo Make this an input variable rather than hardcode @endtodo omega = 1.0d14 n2 = 100000 !! @todo Figure out the purpose of these variables @endtodo if ( id == 0 ) then !! * If root process !!    * Open `Inputfile` to read the number of phonon modes !!      which is needed to allocate variables !!    * Open `input.txt` file to read temperature, `n1`, !!      `limit`, `gamma_p`, `alpha`, the energy of the laser, !!      `elevel`, and `eshift_num` !!    * Open `output.txt` for output !! @todo Change `input.txt` to be read from input file @endtodo !! @todo Change output to go to command line like QE @endtodo open ( 11 , file = trim ( Inputfile ), Action = 'read' , status = 'old' ) read ( 11 , * ) read ( 11 , * ) nmode open ( 12 , file = 'input.txt' , Action = 'read' , status = 'old' ) read ( 12 , * ) read ( 12 , * ) temperature , n1 , limit , gamma_p , alpha , elevel , elaser_num , eshift_num open ( 13 , file = \"output.txt\" , Action = \"write\" , status = \"replace\" ) endif call MPI_BCast ( nmode , 1 , MPI_Integer , 0 , MPI_COMM_WORLD , ierror ) call MPI_BCast ( eshift_num , 1 , MPI_Integer , 0 , MPI_COMM_WORLD , ierror ) call MPI_BCast ( elaser_num , 1 , MPI_Integer , 0 , MPI_COMM_WORLD , ierror ) call MPI_Bcast ( temperature , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD , ierror ) !! * Broadcast the number of modes, `eshift_num`, and temperature to !!  other processes beta = 1 / ( kB * temperature ) !! * Calculate \\beta = 1/k_{B}T allocate ( eshift ( eshift_num ), elaser ( elaser_num ), omega2 ( eshift_num ), omega_l ( elaser_num )) allocate ( s1 ( eshift_num ), s2 ( eshift_num , elaser_num ), s3 ( eshift_num , elaser_num ), global_sum ( eshift_num , elaser_num )) allocate ( Sj ( nmode ), omega_j ( nmode ), omega_nj ( nmode ), hbarOmegaBeta ( nmode ), FjFractionFactor ( nmode ), expX ( nmode ), expY ( nmode )) allocate ( domega ( nmode ), theta ( nmode ), zfactor1 ( nmode ), zfactor2 ( nmode ), ex1 ( 0 : n2 + 1 ), interval ( 2 ), count2 ( 2 )) !! * Allocate space for variables on all processes step2 = tpi / float ( n2 ) !! * Calculate `step2`= 2\\pi/n2 where `n2=100000` do j = 0 , n2 + 1 !! * Calculate e&#94;{i\\theta} where \\theta goes from 0 to 2\\pi !! @todo Figure out what this is used for @endtodo ex1 ( j ) = exp ( I * j * step2 ) end do if ( id == 0 ) then !! * If root process !!    * For each mode !!       * Read index (currently not being used), S_j, !!         \\omega_j, and \\omega_{nj} !!       * Divide the initial and final phonon frequencies by 100 !!       * Calculate \\text{hbarOmegaBeta}=\\hbar\\omega\\beta\\omega_j !!       * Calculate `FjFractionFactor`={\\sin(-i\\hbar\\omega\\beta)}{1-\\cos(-i\\hbar\\omega\\beta)} !!       * Calculate \\delta\\omega_{nj} = \\omega_{nj} - \\omega_j !!       * Write out the id and phonon frequencies !!    * Read in the energy shifts omega_j (:) = omega_j (:) / 10 0.0d0 omega_nj (:) = omega_nj (:) / 10 0.0d0 do imode = 1 , nmode read ( 11 , * ) j , Sj ( imode ), omega_j ( imode ), omega_nj ( imode ) !write(*,*) id, Sj(imode), omega_j(imode), omega_nj(imode) end do hbarOmegaBeta (:) = hbar * omega_j (:) * omega * beta domega (:) = omega_nj (:) - omega_j (:) read ( 12 , * ) eshift (:) read ( 12 , * ) elaser (:) endif call MPI_Bcast ( n1 , 1 , MPI_INTEGER , 0 , MPI_COMM_WORLD , ierror ) call MPI_Bcast ( limit , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD , ierror ) call MPI_Bcast ( gamma_p , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD , ierror ) call MPI_Bcast ( alpha , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD , ierror ) call MPI_Bcast ( elaser , elaser_num , MPI_DOUBLE , 0 , MPI_COMM_WORLD , ierror ) call MPI_Bcast ( elevel , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD , ierror ) call MPI_Bcast ( eshift , eshift_num , MPI_DOUBLE , 0 , MPI_COMM_WORLD , ierror ) call MPI_Bcast ( Sj , nmode , MPI_DOUBLE , 0 , MPI_COMM_WORLD , ierror ) call MPI_Bcast ( omega_j , nmode , MPI_DOUBLE , 0 , MPI_COMM_WORLD , ierror ) call MPI_Bcast ( omega_nj , nmode , MPI_DOUBLE , 0 , MPI_COMM_WORLD , ierror ) call MPI_Bcast ( hbarOmegaBeta , nmode , MPI_DOUBLE , 0 , MPI_COMM_WORLD , ierror ) call MPI_Bcast ( domega , nmode , MPI_DOUBLE , 0 , MPI_COMM_WORLD , ierror ) call MPI_Barrier ( MPI_COMM_WORLD , ierror ) !> Maybe unit conversions? omega_l (:) = ( elaser (:) - elevel ) * ev / hbar / omega !! (E_L-E_n)/\\hbar\\omega omega2 (:) = eshift (:) * mev / hbar / omega gamma_p = gamma_p * mev / hbar / omega alpha = alpha * mev / hbar / omega !> Set loop variables? step1 = tpi / float ( n1 ) loglimit =- log ( limit ) count1 = 0.0 write ( * , * ) \"id\" , eshift (:) !write(*,*),step1 !> Do some sort of sum? What are `gamma_p` and `alpha`? They come from the input file !> I think this is figuring out the intervals for each node to integrate over do j = 0 , int ( loglimit / gamma_p / step1 ) do k = 0 , int ( loglimit / gamma_p / step1 - j ) count1 = count1 + int (( loglimit / step1 - gamma_p * j - gamma_p * k ) / alpha ) enddo end do count2 ( 1 ) = float ( id ) / float ( nprocs ) * count1 count2 ( 2 ) = float ( id + 1 ) / float ( nprocs ) * count1 !write(*,*)\"count1\",count1 index1 = 1 count1 = 0.0 write ( * , * ) \"iteration2\" do j = 0 , int ( loglimit / gamma_p / step1 ) do k = 0 , int ( loglimit / gamma_p / step1 ) - j count1 = count1 + int (( loglimit / step1 - gamma_p * j - gamma_p * k ) / alpha ) enddo if ( count1 >= count2 ( index1 )) then interval ( index1 ) = j index1 = index1 + 1 if ( index1 == 3 ) then exit endif endif end do interval ( 2 ) = interval ( 2 ) - 1 call MPI_Barrier ( MPI_COMM_WORLD , ierror ) !! Make sure that all processes get here before !! moving forward s3 = 0.0d0 do j = interval ( 1 ), interval ( 2 ) !! Integrate over x if ( id == 0 ) then write ( 13 , * ) id , j endif s2 = 0.0d0 x = ( j + 0.5 ) * step1 expX (:) = cos ( omega_nj (:) * x ) + I * sin ( omega_nj (:) * x ) !! @todo Figure out exponentials for F_j use \\omega_{nj} @endtodo do k = 0 , int ( loglimit / gamma_p / step1 ) - j s1 = 0.0d0 y = ( k + 0.5 ) * step1 expY (:) = cos ( omega_nj (:) * y ) - I * sin ( omega_nj (:) * y ) theta (:) = domega (:) * ( x - y ) - I * hbarOmegaBeta (:) FjFractionFactor (:) = sin ( theta (:)) / ( 1 - cos ( theta (:)) ) !! * Calculate the fractional factor in front of the cosines in !!   F_j theta (:) = hbarOmegaBeta (:) + I * domega (:) * ( x - y ) zfactor1 (:) = exp ( 0.5 * theta (:)) / ( exp ( theta (:)) - 1 ) !! Calculate the first fraction in equation 42 zfactor2 (:) = exp ( 0.5 * hbarOmegaBeta (:)) / ( exp ( hbarOmegaBeta (:)) - 1 ) !! @todo Figure out where `zfactor2` comes from @endtodo zfactor = product ( zfactor1 (:) / zfactor2 (:)) interval_t = int (( loglimit / step1 - gamma_p * j - gamma_p * k ) / alpha ) do l = 0 , interval_t t = ( l + 0.5 ) * step1 tmp_exp = 0.0d0 do imode = 1 , nmode tmp_r1 =- omega_j ( imode ) * t / tpi tmp_r1 = ( tmp_r1 - floor ( tmp_r1 )) * float ( n2 ) tmp_i = floor ( tmp_r1 ) tmp_r2 = tmp_r1 - tmp_i expT = ( 1.0 - tmp_r2 ) * ex1 ( tmp_i ) + tmp_r2 * ex1 ( tmp_i + 1 ) !! Use a linear interpolation for e&#94;{i\\omega t} expForFj =- expT * ( 1 - expX ( imode )) * ( 1 - expY ( imode )) - ( expX ( imode ) + expY ( imode )) !! * Calculate `expForFj` = -e&#94;{-i\\omega t}(1-e&#94;{i\\omega x})(1-e&#94;{-i\\omega y})-(e&#94;{i\\omega x} + e&#94;{-i\\omega y}). !!   This is used as a trick to be able to calculate F_j quicker as the expontentials include both the !!   sines and cosines needed Fj = Aimag ( expForFj ) + FjFractionFactor ( imode ) * Real ( 2.0d0 + expForFj ) !! * Calculate F_j = \\text{Im}(`expForFj`) + `FjFractionFactor`\\text{Re}(2 + `expForFj`) !! @todo Add detailed derivation of this in a separate page @endtodo tmp_exp = tmp_exp + Fj * Sj ( imode ) enddo s1 (:) = s1 (:) + exp ( I * tmp_exp - I * omega2 (:) * t - alpha * abs ( t )) end do do ilaserE = 1 , elaser_num s2 (:, ilaserE ) = s2 (:, ilaserE ) + s1 (:) * exp ( - ( I * omega_l ( ilaserE ) + gamma_p ) * y ) * zfactor enddo end do do ilaserE = 1 , elaser_num s3 (:, ilaserE ) = s3 (:, ilaserE ) + s2 (:, ilaserE ) * exp ( - ( - I * omega_l ( ilaserE ) + gamma_p ) * x ) enddo enddo call MPI_Barrier ( MPI_COMM_WORLD , ierror ) call MPI_Reduce ( s3 , global_sum , eshift_num * elaser_num , MPI_DOUBLE_COMPLEX , MPI_SUM , 0 , MPI_COMM_WORLD , ierror ) if ( id == 0 ) then write ( 13 , * ) \"calculation finalized\" do ilaserE = 1 , elaser_num write ( 13 , * ) \"Laser energy: \" , elaser ( ilaserE ) do j = 1 , eshift_num write ( 13 , * ) eshift ( j ), eshift ( j ) * mevtocm , step1 ** 3 * Real ( global_sum ( j , ilaserE )) * 2.0 enddo enddo endif call MPI_FINALIZE ( ierror ) end program ramanIntensity","tags":"","loc":"sourcefile/ramanintensity.f90.html","title":"ramanIntensity.f90 – Raman Intensities"},{"text":"Contents Variables dp ev hbar kB mev mevtocm pi tpi I Variables Type Visibility Attributes Name Initial integer, public, parameter :: dp = selected_real_kind(15, 307) real(kind=dp), public, parameter :: ev = 1.6021766d-19 Conversion factor from eV to J real(kind=dp), public, parameter :: hbar = 1.0545718d-34 \\bar{h} real(kind=dp), public, parameter :: kB = 1.38064852d-23 Boltzmann constant real(kind=dp), public, parameter :: mev = 1.6021766d-22 Conversion factor from meV to J real(kind=dp), public, parameter :: mevtocm = 8.0655438354 Conversion factor from meV to cm &#94;{-1} real(kind=dp), public, parameter :: pi = 3.1415926535897932_dp \\pi real(kind=dp), public, parameter :: tpi = 2.0*3.1415926535897932_dp 2\\pi complex(kind=dp), public, parameter :: I = cmplx(0.0d0, 1.0d0, dp) i","tags":"","loc":"module/para.html","title":"para – Raman Intensities"},{"text":"Uses para mpi Initialize MPI pool Define input file name Todo Make this an input variable rather than hardcode Todo Figure out the purpose of these variables If root process Open Inputfile to read the number of phonon modes\n  which is needed to allocate variables Open input.txt file to read temperature, n1 , limit , gamma_p , alpha , the energy of the laser, elevel , and eshift_num Open output.txt for output","tags":"","loc":"program/ramanintensity.html","title":"ramanIntensity – Raman Intensities"},{"text":"Welcome to the README for the Raman Intensities project! This tab contains more information\non the theory behind the code, the input and output files used by the code, and how to install/run\nthe code as well as a Todo list for cleaning up and developing the code.","tags":"","loc":"page//index.html","title":"README – Raman Intensities"},{"text":"../src/ramanIntensity.f90:137: Make this an input variable rather than hardcode ../src/ramanIntensity.f90:140: Figure out the purpose of these variables ../src/ramanIntensity.f90:150: Change input.txt to be read from input file ../src/ramanIntensity.f90:151: Change output to go to command line like QE ../src/ramanIntensity.f90:188: Figure out what this is used for ../src/ramanIntensity.f90:303: Figure out exponentials for F_j use \\omega_{nj} ../src/ramanIntensity.f90:320: Figure out where zfactor2 comes from ../src/ramanIntensity.f90:345: Add detailed derivation of this in a separate page","tags":"","loc":"page/./todo.html","title":"Todo – Raman Intensities"}]}