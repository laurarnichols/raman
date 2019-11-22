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

integer :: eshift_num
integer :: ierror
  !! Error for MPI
integer :: id
  !! MPI process id
integer :: imode
integer :: index1
integer :: interval_a
integer :: interval_b
integer :: j
integer :: k
integer :: l
integer :: max_k
integer :: max_l
integer :: n1
integer :: n2
  !! Hardcoded to be 100000; used to generate `ex1`
integer :: nmode
  !! Number of phonon modes
integer :: nprocs
  !! Number of MPI processes
integer :: tmp_i
integer :: tmp_j

real(kind = dp) :: beta
  !! \(\beta = 1/k_{B}T\)
real(kind = dp) :: count1
real(kind = dp) :: domega
real(kind = dp) :: elaser
  !! \(E_L\)
real(kind = dp) :: elevel
  !! \(E_n\)
real(kind = dp) :: gamma1
real(kind = dp) :: gamma2
real(kind = dp) :: loglimit
real(kind = dp) :: limit
real(kind = dp) :: omega
real(kind = dp) :: omega1
real(kind = dp) :: step1
real(kind = dp) :: step2
  !! Step to go from 0 to \(\2\pi) in `n2` steps
real(kind = dp) :: t
real(kind = dp) :: temperature
real(kind = dp) :: tmp_r1
real(kind = dp) :: tmp_r2
real(kind = dp) :: x
real(kind = dp) :: y

complex(kind = dp) :: expForFj
  !! An exponential form to more quickly
  !! calculate \(F_j\)
complex(kind = dp) :: expT
  !! \(e^{i\omega t}\) used in calculating \(F_j\)
complex(kind = dp) :: Fj
  !! Function \(F_j\) from equation 44
complex(kind = dp) :: theta
  !! Serves as the argument for the fraction
  !! factor `FjFractionFactor`
complex(kind = dp) :: tmp
complex(kind = dp) :: tmp1
complex(kind = dp) :: tmp_exp
complex(kind = dp) :: zfactor
complex(kind = dp) :: zfactor1
complex(kind = dp) :: zfactor2 

character(len = 256) :: Inputfile
  !! Input file name
character(len = 256) :: dummy

integer,allocatable :: interval(:)

real(kind = dp),allocatable :: count2(:)
real(kind = dp), allocatable :: eshift(:)
real(kind = dp),allocatable :: factor(:)
real(kind = dp), allocatable :: omega2(:)
real(kind = dp), allocatable :: omega_j(:)
  !! \(\omega_j\)
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
complex(kind = dp), allocatable :: global_sum(:)
complex(kind = dp), allocatable :: s1(:)
complex(kind = dp), allocatable :: s2(:)
complex(kind = dp), allocatable :: s3(:)



call MPI_INIT(ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, id, ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
  !! * Initialize MPI pool

Inputfile = 'Sj.out'
  !! * Define input file name
  !! @todo Make this an input variable rather than hardcode @endtodo
omega=1.0d14
n2=100000

if(id == 0) then
  !! * If root process
  !!    * Open `Inputfile` to read the number of phonon modes
  !!      which is needed to allocate variables
  !!    * Open `input.txt` file to read temperature, `n1`, 
  !!      `limit`, `gamma1`, `gamma2`, the energy of the laser,
  !!      `elevel`, and `eshift_num`
  !!    * Open `output.txt` for output
  !! @todo Change `input.txt` to be read from input file @endtodo
  !! @todo Change output to go to command line like QE @endtodo

  open(11, file=trim(Inputfile), Action='read', status='old')
  read(11,*)
  read(11,*)nmode

  open(12 , file='input.txt', Action='read', status='old')
  read(12,*) dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy
    !! @todo Change this to `read(12,*)` @endtodo
  read(12,*) temperature, n1, limit, gamma1, gamma2, elaser, elevel, eshift_num

  open(13,file="output.txt",Action="write",status="replace")
endif


call MPI_BCast( nmode, 1, MPI_Integer, 0, MPI_COMM_WORLD, ierror)
call MPI_BCast( eshift_num, 1, MPI_Integer, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( temperature, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
  !! * Broadcast the number of modes, `eshift_num`, and temperature to
  !!  other processes


beta=1/(kB*temperature)
  !! * Calculate \(\beta = 1/k_{B}T\)


allocate(eshift(eshift_num), omega2(eshift_num), s1(eshift_num), s2(eshift_num), s3(eshift_num), global_sum(eshift_num))
allocate(Sj(nmode), omega_j(nmode), omega_nj(nmode), factor(nmode), FjFractionFactor(nmode), expX(nmode), expY(nmode))
allocate(ex1(0:n2+1), interval(2), count2(2))
  !! * Allocate space for variables on all processes
  !! @todo Figure out why `interval` and `count2` are allocatable @endtodo

step2=tpi/float(n2)
  !! * Calculate `step2`\(= 2\pi/n2\) where `n2=100000`

do j=0,n2+1
  !! * Calculate \(e^{i\theta}\) where \(\theta\) goes from 0 to \(2\pi\)

   ex1(j)=exp(I*j*step2)
end do

if(id == 0) then
  !! * If root process
  !!    * For each mode
  !!       * Read index (currently not being used), \(S_j\),
  !!         \(\omega_j\), and \(\omega_{nj}\)
  !!       * Divide the initial and final phonon frequencies by 100
  !!       * Calculate \(\text{factor}=\hbar\omega\beta\text{phonon}(2)\)
  !!       * Calculate `FjFractionFactor`\(={\sin(-i\text{factor})}{1-\cos(-i\text{factor})}\)
  !!       * Write out the id and phonon frequencies

  do imode=1,nmode
     
     read(11,*) j, Sj(imode), omega_j(imode), omega_nj(imode)

     omega_j(imode)=omega_j(imode)/100.0d0
     omega_nj(imode)=omega_nj(imode)/100.0d0     

     factor(imode)=hbar*omega_j(imode)*omega*beta
      
     write(*,*)id, Sj(imode), omega_j(imode), omega_nj(imode)

  end do
  
     read(12,*) eshift(:)

endif

call MPI_Bcast( n1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( limit, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( gamma1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( gamma2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( elaser, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( elevel, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( eshift, eshift_num, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( Sj, nmode, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( omega_j, nmode, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( omega_nj, nmode, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( factor, nmode, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror) 
call MPI_Barrier(MPI_COMM_WORLD,ierror)

!> Maybe unit conversions?
omega1=(elaser-elevel)*ev/hbar/omega
  !! \((E_L-E_n)/\hbar\omega\)
omega2(:)=eshift(:)*mev/hbar/omega
gamma1=gamma1*mev/hbar/omega
gamma2=gamma2*mev/hbar/omega

!> Set loop variables?
step1=tpi/float(n1)
loglimit=-log(limit)
count1=0.0

write(*,*)"id", eshift(:)
!write(*,*),step1

!> Do some sort of sum? What are gamma1 and gamma2? They come from the input file
!> I think this is figuring out the intervals for each node to integrate over
do j=0,int(loglimit/gamma1/step1)
   do k=0,int(loglimit/gamma1/step1-j)
      count1=count1+int((loglimit/step1-gamma1*j-gamma1*k)/gamma2)
   enddo
end do

count2(1) = float(id)/float(nprocs)*count1
count2(2) = float(id+1)/float(nprocs)*count1
!write(*,*)"count1",count1
index1=1
count1=0.0
write(*,*)"iteration2"
do j=0,int(loglimit/gamma1/step1)
   do k=0,int(loglimit/gamma1/step1)-j
      count1=count1+int((loglimit/step1-gamma1*j-gamma1*k)/gamma2)
   enddo

   if(count1>=count2(index1)) then
      interval(index1)=j
      index1=index1+1
      if(index1==3 ) then
         exit
      endif
   endif
end do

interval(2)=interval(2)-1


!write(*,*)"interval1", id, interval(1)
!write(*,*)"interval2", id,  interval(2)
call MPI_Barrier(MPI_COMM_WORLD,ierror)
  !! Make sure that all processes get here before
  !! moving forward
!max_k=floor(loglimit/gamma1/step1)
!max_l=floor(loglimit/gamma2/step1)
!interval_a = floor(loglimit/gamma1/step1*float(id)/float(nprocs))
!interval_b = floor(loglimit/gamma1/step1*float(id+1)/float(nprocs))-1

s3=0.0d0
do j=interval(1),interval(2)
  !! Integrate over \(x\)

   if(id==0) then
      write(13,*)id,j
   endif

   s2 = 0.0d0
   x = (j+0.5) * step1

   do imode=1,nmode
 
      expX(imode)=cos(omega_nj(imode)*x)+I*sin(omega_nj(imode)*x)
 
   end do
 
   do k=0,int(loglimit/gamma1/step1)-j 
      !write(*,*)id,k,int(loglimit/gamma1/step1)-j
      s1 = 0.0d0
      y = (k+0.5) * step1


      zfactor = 1.0
      do imode=1,nmode
          !! @todo Take this out of the loop @endtodo
         domega = omega_nj(imode)  - omega_j(imode)
         expY(imode) = cos( omega_nj(imode)*y ) - I * sin(omega_nj(imode)*y)
         
         tmp = factor(imode)  
         tmp1 = tmp + I*domega * ( x - y ) 
         
         zfactor1 = exp(0.5*tmp1) / ( exp(tmp1) - 1 )        
          !! Calculate the first fraction in equation 42
         zfactor2 =exp(0.5*tmp) / ( exp(tmp) - 1 )
          !! @todo Figure out where `zfactor2` comes from
         zfactor = zfactor * zfactor1 / zfactor2
          !! @todo Since multiple them, figure out why have `exp(0.5*tmp)` as it cancels @endtodo


         theta = domega*( x - y ) - I*theta
         FjFractionFactor(imode) = sin(theta)/( 1 - cos(theta)  )
          !! * Calculate the fractional factor in front of the cosines in
          !!   \(F_j\)
      enddo
  
      
              
      

      interval_b=int((loglimit/step1-gamma1*j-gamma1*k)/gamma2)
!      interval_a=-interval_b

      do l= 0, interval_b!-max_l,max_l
         t=(l+0.5)*step1
         tmp_exp=0.0d0
 
         do imode=1,nmode

            tmp_r1=-omega_j(imode)*t/tpi
            tmp_r1=(tmp_r1-floor(tmp_r1))*float(n2)
            tmp_i=floor(tmp_r1)
            tmp_r2=tmp_r1-tmp_i
            expT=(1.0-tmp_r2)*ex1(tmp_i)+tmp_r2*ex1(tmp_i+1)
              !! Use a linear interpolation for \(e^{i\omega t}\)

            expForFj=-expT*(1-expX(imode))*(1-expY(imode))-(expX(imode)+expY(imode))
              !! * Calculate `expForFj`\( = -e^{-i\omega t}(1-e^{i\omega x})(1-e^{-i\omega y})-(e^{i\omega x} + e^{-i\omega y})\).
              !!   This is used as a trick to be able to calculate \(F_j\) quicker as the expontentials include both the
              !!   sines and cosines needed 
            Fj=Aimag(expForFj)+FjFractionFactor(imode)*Real(2.0d0+expForFj)           
              !! * Calculate \(F_j = \text{Im}(\)`expForFj`\() + \)`FjFractionFactor`\(\text{Re}(2 + \)`expForFj`\()\)
              !! @todo Add detailed derivation of this in a separate page @endtodo
            tmp_exp=tmp_exp+Fj*Sj(imode)
         enddo

         s1(:)=s1(:)+exp(I*tmp_exp-I*omega2(:)*t-gamma2*abs(t))
      end do


      s2(:)=s2(:)+s1(:)*exp(-(I*omega1+gamma1)*y) * zfactor
   end do


   s3(:)=s3(:)+s2(:)*exp(-(-I*omega1+gamma1)*x)
enddo


call MPI_Barrier(MPI_COMM_WORLD,ierror) 
call MPI_Reduce( s3, global_sum, eshift_num, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierror)



if(id == 0) then
  write(13,*)"calculation finalized"
  do j = 1, eshift_num 
     write(13,*) eshift(j), eshift(j)*mevtocm, step1**3*Real(global_sum(j))*2.0
  end do
endif
call MPI_FINALIZE(ierror)

end program ramanIntensity

         
