module para
  integer, parameter :: dp = selected_real_kind(15, 307)
  real(kind = dp), parameter :: pi  =  3.1415926535897932_dp
    !! \(\pi\)
  real(kind = dp), parameter :: tpi  =  2.0*3.1415926535897932_dp
    !! \(2\pi\)
  real(kind = dp), parameter :: hbar = 1.0545718d-34
    !! \(\bar{h}\)
  real(kind = dp), parameter :: ev = 1.6021766d-19
    !! Conversion factor from eV to J
  real(kind = dp), parameter :: mev = 1.6021766d-22
    !! Conversion factor from meV to J
  real(kind = dp), parameter :: mevtocm = 8.0655438354 
    !! Conversion factor from meV to cm\(^{-1}\)
  complex(kind = dp), parameter :: I = cmplx(0.0d0, 1.0d0, dp)  
    !! \(i\)
end module

program lsf
use para
implicit none
include "mpif.h"

integer :: eshift_num
integer :: ierror
integer :: id
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
integer :: nmode
integer :: nprocs
integer :: tmp_i
integer :: tmp_j

real(kind = dp) :: beta
real(kind = dp) :: count1
real(kind = dp) :: domega
real(kind = dp) :: elaser
real(kind = dp) :: elevel
real(kind = dp) :: gamma1
real(kind = dp) :: gamma2
real(kind = dp) :: loglimit
real(kind = dp) :: limit
real(kind = dp) :: omega
real(kind = dp) :: omega1
real(kind = dp) :: omega_tmp
real(kind = dp) :: step1
real(kind = dp) :: step2
real(kind = dp) :: t
real(kind = dp) :: temperature
real(kind = dp) :: tmp_r1
real(kind = dp) :: tmp_r2
real(kind = dp) :: x
real(kind = dp) :: y

complex(kind = dp) :: T2! T1,T3
complex(kind = dp) :: tmp
complex(kind = dp) :: tmp1
complex(kind = dp) :: tmp_exp
complex(kind = dp) :: zfactor
complex(kind = dp) :: zfactor1
complex(kind = dp) :: zfactor2 

character(len = 256) :: Inputfile
character(len = 256) :: dummy

integer,allocatable :: interval(:)

real(kind = dp),allocatable :: count2(:)
real(kind = dp), allocatable :: eshift(:)
real(kind = dp),allocatable :: factor(:)
real(kind = dp), allocatable :: omega2(:)
real(kind = dp),allocatable :: phonon(:,:)

complex(kind = dp), allocatable :: ex1(:)
complex(kind = dp), allocatable :: factor1(:)
complex(kind = dp), allocatable :: global_sum(:)
complex(kind = dp), allocatable :: s1(:)
complex(kind = dp), allocatable :: s2(:)
complex(kind = dp), allocatable :: s3(:)
complex(kind = dp), allocatable :: T1(:)
complex(kind = dp), allocatable :: T3(:)



call MPI_INIT(ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, id, ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)


Inputfile = 'Sj.out'
omega=1.0d14
n2=100000

if(id == 0) then
  open(11, file=trim(Inputfile), Action='read', status='old')
  open(12 , file='input.txt', Action='read', status='old')
  read(11,*)
  read(11,*)nmode
  read(12,*) dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy
  read(12,*) temperature, n1, limit, gamma1, gamma2, elaser, elevel, eshift_num
  open(13,file="output.txt",Action="write",status="replace")
endif


call MPI_BCast( nmode, 1, MPI_Integer, 0, MPI_COMM_WORLD, ierror)
call MPI_BCast( eshift_num, 1, MPI_Integer, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( temperature, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)



beta=1/(temperature*1.38064852d-23)


allocate(eshift(1:eshift_num))
allocate(omega2(1:eshift_num))
allocate(s1(1:eshift_num))
allocate(s2(1:eshift_num))
allocate(s3(1:eshift_num))
allocate(global_sum(1:eshift_num))
allocate(phonon(nmode,3))
allocate(factor(1:nmode))
allocate(factor1(1:nmode))
allocate(T1(1:nmode))
allocate(T3(1:nmode))
allocate(ex1(0:n2+1))
allocate(interval(1:2))
allocate(count2(1:2))

step2=tpi/float(n2)

do j=0,n2+1
   ex1(j)=exp(I*j*step2)
end do

if(id == 0) then
  do imode=1,nmode
     
     read(11,*)j,phonon(imode,1),phonon(imode,2), phonon(imode, 3)

     phonon(imode,2)=phonon(imode,2)/100.0d0
     phonon(imode,3)=phonon(imode,3)/100.0d0     
     factor(imode)=hbar*phonon(imode,2)*omega*beta
     factor1(imode)=sin(-I*factor(imode))/(1-cos(-I*factor(imode)))
      
     write(*,*)id, phonon(imode,1),phonon(imode,2),phonon(imode,3)

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
call MPI_Bcast( phonon, 3*nmode, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
call MPI_Bcast( factor, nmode, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror) 
call MPI_Bcast( factor1, nmode, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)
call MPI_Barrier(MPI_COMM_WORLD,ierror)
omega1=(elaser-elevel)*ev/hbar/omega
omega2(:)=eshift(:)*mev/hbar/omega
gamma1=gamma1*mev/hbar/omega
gamma2=gamma2*mev/hbar/omega
step1=tpi/float(n1)
loglimit=-log(limit)
count1=0.0
write(*,*)"id", eshift(:)
!write(*,*),step1
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
!max_k=floor(loglimit/gamma1/step1)
!max_l=floor(loglimit/gamma2/step1)
!interval_a = floor(loglimit/gamma1/step1*float(id)/float(nprocs))
!interval_b = floor(loglimit/gamma1/step1*float(id+1)/float(nprocs))-1
s3=0.0d0
do j=interval(1),interval(2)
   if(id==0) then
      write(13,*)id,j
   endif

   s2 = 0.0d0
   x = (j+0.5) * step1

   do imode=1,nmode
 
      omega_tmp=phonon(imode,3)
      T1(imode)=cos(omega_tmp*x)+I*sin(omega_tmp*x)
 
   end do
 
   do k=0,int(loglimit/gamma1/step1)-j 
      !write(*,*)id,k,int(loglimit/gamma1/step1)-j
      s1 = 0.0d0
      y = (k+0.5) * step1


      zfactor = 1.0
      do imode=1,nmode
         omega_tmp = phonon(imode,3)
         domega = omega_tmp  - phonon(imode,2)
         T3(imode) = cos( omega_tmp*y ) - I * sin(omega_tmp*y)
         
         tmp = factor(imode)  
         tmp1 = tmp + I*domega * ( x - y ) 
         
         zfactor1 = exp(0.5*tmp1) / ( exp(tmp1) - 1 )        
         zfactor2 =exp(0.5*tmp) / ( exp(tmp) - 1 )
         zfactor = zfactor * zfactor1 / zfactor2


         tmp = domega*( x - y ) - I*tmp
         factor1(imode) = sin(tmp)/( 1 - cos(tmp)  )
      enddo
  
      
              
      

      interval_b=int((loglimit/step1-gamma1*j-gamma1*k)/gamma2)
!      interval_a=-interval_b

      do l= 0, interval_b!-max_l,max_l
         t=(l+0.5)*step1
         tmp_exp=0.0d0
 
         do imode=1,nmode
            omega_tmp = phonon(imode,2)
            tmp_r1=-omega_tmp*t/tpi
            tmp_r1=(tmp_r1-floor(tmp_r1))*float(n2)
            tmp_i=floor(tmp_r1)
            tmp_r2=tmp_r1-tmp_i
            T2=(1.0-tmp_r2)*ex1(tmp_i)+tmp_r2*ex1(tmp_i+1)
            tmp=-T2*(1-T1(imode))*(1-T3(imode))-(T1(imode)+T3(imode))
            tmp1=Aimag(tmp)+factor1(imode)*Real(2.0d0+tmp)           
            tmp_exp=tmp_exp+tmp1*phonon(imode,1)
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
     write(13,*),eshift(j), eshift(j)*mevtocm, step1**3*Real(global_sum(j))*2.0
  end do
endif
call MPI_FINALIZE(ierror)

end program

