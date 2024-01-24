!
!    mlff_select: selection of an effective set of 
!      local reference configurations (basis sets) for 
!      a VASP machine-learning force field from a given
!      ML_AB file. The procedure is similar to the VASP 
!      command ML_MODE = select, but is faster and has more
!      options to be adjusted.
!    Part of VASP4CLINT
!     Julien Steffen, 2024 (julien.steffen@fau.de)
!
program mlff_select
implicit none
integer::i,j,k,l   ! loop indices
integer::inc,inc2,inc3  ! incremented index
integer::readstat
integer::nbasis  ! desired number of local refs. in new ML_AB
integer::conf_num  ! number of configurations (full structures)
integer::act_conf   ! the current configuration number
integer::natoms_max  ! maximum number of atoms per configuration
integer::natoms_sum  ! sum of all atoms within all configurations
integer::nelems  ! current number of elements
integer::nelems_all  ! total number of different elements
integer::nelems_glob  ! number of elements in current ML_AB
integer::ind_act ! current element index
integer::max_environ  ! maximum number of atoms in environment
integer::max_around  ! maximum actual number of atoms around
integer::ngrid  ! number of grid points for RDF calculation
integer,allocatable::natoms(:)  ! atom numbers of configurations
real(kind=8),allocatable::xyz(:,:,:) ! the geometries
real(kind=8),allocatable::xyz_dir(:,:,:) ! the geometries (direct coords)
real(kind=8),allocatable::cells(:,:,:) ! the unit cells
real(kind=8),allocatable::angles(:)  ! list of angles for current environment
integer,allocatable::inds(:,:)  ! the element indices (core charges)
integer::ind_list(50)   ! list of different element indices
character(len=2),allocatable::el_list(:)  ! current list of elements
character(len=2)::el_list_glob(50)  ! global list of elements 
character(len=2)::el_act   ! the current element symbol for printing
character(len=120)::arg  ! command line argument
real(kind=8)::train_div   ! desired trainset diversity (0: none, 1: full)
real(kind=8),allocatable::energies(:)  ! the energies 
real(kind=8),allocatable::grads(:,:,:)  ! the gradient vectors
real(kind=8),allocatable::environ(:,:,:)  ! the local environment of each atom
real(kind=8),allocatable::dist_list(:)  ! list of distances in environment
real(kind=8),allocatable::rdf_all(:,:)  ! radial distribution functions for atoms
real(kind=8),allocatable::adf_all(:,:)  ! angular distribution functions for atoms
real(kind=8)::cell_inv(3,3)  ! inverted coordinate unit cell
real(kind=8)::dist_act   ! scalar distance between two atoms
real(kind=8)::cutoff    ! the distance cutoff during the ML-FF learning
real(kind=8)::dist_vec(3)  ! distance vector between two atoms
real(kind=8)::dx   ! distance between two gridpoints (RDF)
real(kind=8)::da   ! distance between two gridpoints (ADF)
real(kind=8)::alpha_r   ! exponent for RDF Gaussian functions
real(kind=8)::alpha_a   ! exponent for ADF Gaussian functions
real(kind=8)::x_act   ! current x-value for RDF buildup
real(kind=8)::pi  ! the Pi
integer,allocatable::el_nums(:)  ! current numbers of elements
integer,allocatable::confnum_all(:)  ! numbers of configurations 
integer,allocatable::nat_all(:) ! number of atom in configurations
integer,allocatable::ind_all(:)  ! element indices of all atoms
integer,allocatable::num_around(:,:) ! number of atoms around the atom
integer,allocatable::ind_env(:,:)  ! core charges of atoms in environments
integer,allocatable::neigh_global(:,:)  ! list with numbers of neighbors (global)
logical::eval_stat(10)  ! the progress for evaluation loops

real(kind=8),allocatable::gradnorm_all(:)  ! gradient forms for all atoms

character(len=130)::a130


write(*,*) 
write(*,*) "PROGRAM mlff_select: selection of an effective set of"
write(*,*) " local reference configurations (basis sets) for a VASP"
write(*,*) " machine-learning force field from a given ML-AB file."
write(*,*) "The procedure is similar to the VASP command ML_MODE = select,"
write(*,*) " but is faster (since purely based on geometrical comparisons)"
write(*,*) " and has more options to be adjusted (global approach)"
write(*,*) "Usage: mlff_select -nbasis=[number], where [number] = desired"
write(*,*) " number of local reference configurations (basis functions) "
write(*,*) " per element in the newly written ML_AB file."
write(*,*)


!
!     Read in command line arguments
!
!     The desired number of reference environments (basis functions)
!
nbasis=0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:8))  .eq. "-nbasis=") then
      read(arg(9:),*,iostat=readstat) nbasis
      if (readstat .ne. 0) then
         write(*,*) "The format of the -nbasis=[number] command seems to be corrupted!"
         write(*,*)
         stop
      end if
   end if
end do
if (nbasis .lt. 1) then
   write(*,*) "Please give the number of local reference configurations after selection"
   write(*,*) " with the keyword -nbasis=[number]! Recommended are 2000-6000"
   write(*,*)
   stop
end if       
!
!     The training set diversity (1: max, 0 min)
!
train_div=-1.d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:11))  .eq. "-train_div=") then
      read(arg(12:),*,iostat=readstat) nbasis
      if (readstat .ne. 0) then
         write(*,*) "The format of the -train_div=[value] command seems to be corrupted!"
         write(*,*)
         stop
      end if
   end if
end do
if (nbasis .lt. 0.0 .or. train_div .gt. 1.0) then
   write(*,*) "Please give the desired trainset diversity with the keyword -train_div=[value]"
   write(*,*) " where the maximum diversity is 1.0 and the minimum is 0.0"
   write(*,*)
   stop
end if





write(*,'(a,i7)') " Number of desired local reference configuraitons per element:",nbasis

!
!     Define the Pi
!
pi=4d0*atan(1.0d0)
!
!     Default value for the cutoff: 5 Ang
!
cutoff=5.d0
!
!     Default values for global elements array
!
el_list_glob="XX"
!
!     Default values for global core charge array
!
ind_list=0
nelems_all=0
!
!     Maximum number of atoms in environment (assumed)
!
max_environ=100
!
!     Number of grid points of numerical radial distribution functions
!
ngrid=500
!
!     Length between two x values on RDF grid (Ang.)
!
dx=cutoff/real(ngrid)
!
!     Length between two angle values on ADF grid (degrees)
!
da=180.d0/real(ngrid)
!
!     Exponent for RDF Gaussian
!
alpha_r=20.0
!
!     Exponent for ADF Gaussian
!
alpha_a=0.5d0
!
!     Open the ML_AB file and check if it's there
!
open(unit=56,file="ML_AB",iostat=readstat)
if (readstat .ne. 0) then
   write(*,*) "The file ML_AB could not been found!"
   stop
end if
write(*,*) "Read in the ML_AB file..."
!
!     Read in the ML_AB file
!
!     First, read in the header and determine number of atoms and 
!      configurations
read(56,*)
read(56,*)
read(56,*)
read(56,*)
read(56,*) conf_num
read(56,*)
read(56,*)
read(56,*)
read(56,*) nelems_glob
read(56,*)
read(56,*)
read(56,*)
read(56,*,iostat=readstat) el_list_glob
!
!     Determine different element indices
!
outer: do i=1,50
   if (el_list_glob(i) .eq. "XX") exit
   call elem(el_list_glob(i),ind_act)
   inner: do j=1,nelems_all
      if (ind_act .eq. ind_list(j)) cycle outer
   end do inner 
   ind_list(i) = ind_act
   nelems_all=nelems_all+1
end do outer

read(56,*)
read(56,*)
read(56,*)
read(56,*) natoms_max
!
!     Allocate global arrays 
!
!     The number of atoms in all reference structures
allocate(natoms(conf_num)) 
!     The geometries (cartesian coordinates)
allocate(xyz(3,natoms_max,conf_num))
!     The geometries (direct coordinates
allocate(xyz_dir(3,natoms_max,conf_num))
!     The unit cell shapes
allocate(cells(3,3,conf_num))
!     The element indices (core charges), more effective then names
allocate(inds(natoms_max,conf_num))
!     The reference energies
allocate(energies(conf_num))
!     The reference gradients
allocate(grads(3,natoms_max,conf_num))
!     
act_conf=0
do 
   read(56,'(a)',iostat=readstat) a130
   if (readstat .ne. 0) then
      exit
   end if
!
!     If a new configuration has been found, read it in
!
   if (index(a130,"Configuration") .ne. 0) then
      act_conf=act_conf+1
      do i=1,7
         read(56,*)
      end do 
!
!     Read the list and numbers of elements to fulle the inds array
!
      read(56,*) nelems
      if (allocated(el_list)) deallocate(el_list)
      if (allocated(el_nums)) deallocate(el_nums)
      allocate(el_list(nelems))
      allocate(el_nums(nelems))
      do i=1,3
         read(56,*)
      end do
!
!     The current number of atoms
!
      read(56,*) natoms(act_conf)
      read(56,*)
      read(56,*)
      read(56,*)
      do i=1,nelems
         read(56,*) el_list(i),el_nums(i) 
      end do
!
!     Fill element index array
!
      inc=0
      do i=1,nelems
         do j=1,el_nums(i)
            inc=inc+1
            call elem(el_list(i),ind_act)
            inds(inc,act_conf)=ind_act
         end do
      end do      
      do i=1,7
         read(56,*) 
      end do
!
!     Read in the current unit cell vectors
!
      do i=1,3
         read(56,*) cells(:,i,act_conf)
      end do
      read(56,*)
      read(56,*)
      read(56,*)
!
!     Read in the current geometry
!
      do i=1,natoms(act_conf)
         read(56,*) xyz(:,i,act_conf)
      end do
!
!     Invert the unit cell matrix and calculate the geometry 
!        in direct coordinates
!
      call matinv3(cells(:,:,act_conf),cell_inv)
      do i=1,natoms(act_conf)
         xyz_dir(:,i,act_conf)=xyz(1,i,act_conf)*cell_inv(1,:)+ &
              & xyz(2,i,act_conf)*cell_inv(2,:)+xyz(3,i,act_conf)*cell_inv(3,:)
      end do 
      read(56,*)
      read(56,*)
      read(56,*)
!
!     Read in the current energy
!
      read(56,*) energies(act_conf)
      read(56,*)
      read(56,*)
      read(56,*)
!
!     Read in the current gradient
!
      do i=1,natoms(act_conf)
         read(56,*) grads(:,i,act_conf)
      end do

   end if

end do
close(56)
write(*,*) " ... done!"
!
!     Setup classifier arrays for the atoms
!     All atoms (no matter to which configuation they belong)
!     are stored into one large array for each classifier!
!
!     First, allocation of arrays:
!
!     Total number of atoms within any of the configurations
natoms_sum=sum(natoms)
!     The configuration number of the element
allocate(confnum_all(natoms_sum))
!     The atom number in the respective configuration
allocate(nat_all(natoms_sum))
!     The element index (core charge) of each atom
allocate(ind_all(natoms_sum))     
!     The gradient norm for each atom
allocate(gradnorm_all(natoms_sum))
!     The number of atoms in the surrounding (sorted by core charges)
allocate(num_around(nelems,natoms_sum))
!     The coordinates of atoms in the surrounding (for xyz printout)
allocate(environ(3,max_environ,natoms_sum))
!     The core charges/elements of atoms in surrounding
allocate(ind_env(max_environ,natoms_sum))
!     The list of distances to atoms within surrounding (local)
allocate(dist_list(max_environ))
!     The radial distribution functions around all atoms
allocate(rdf_all(ngrid,natoms_sum))
!     The angular distribution functions around all atoms
allocate(adf_all(ngrid,natoms_sum))
!     The local list of angles in the environment
allocate(angles(max_environ*max_environ))
!     The global histogram with the number of neighbors
allocate(neigh_global(nelems,max_environ))

num_around=0
rdf_all=0.d0
adf_all=0.d0
write(*,*) "Total number of atoms (possible basis functions):",natoms_sum

write(*,*)
write(*,*) "Calculate classifiers for all atoms in the ML_AB file..."
eval_stat=.false.
inc=0
do i=1,conf_num
   do j=1,10
      if (real(i)/real(conf_num) .gt. real(j)*0.1d0) then
         if (.not. eval_stat(j)) then
            write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
            eval_stat(j) = .true.
         end if
      end if
   end do

   do j=1,natoms(i)
      inc=inc+1
      confnum_all(inc)=conf_num
      nat_all(inc)=j
      ind_all(inc)=inds(j,i)
!
!     Calculate the gradient norm
!
      gradnorm_all(inc)=sqrt(grads(1,j,i)*grads(1,j,i)+grads(2,j,i)*&
                     & grads(2,j,i)+grads(3,j,i)*grads(3,j,i))
!
!     Calculate the number of atoms in the direct surrounding
!
!     The atom itself is located in the origin for environment
!     If an environment is smaller than the largest one, set all missing
!     atoms to the origin as well, same element as central atom
!
      environ(:,:,inc)=0.d0
      ind_env(:,inc)=ind_all(inc)
!
!     All distances below cutoff
!
      inc2=1
      do k=1,natoms(i)
         if (k .ne. j) then
            dist_vec=xyz_dir(:,j,i)-xyz_dir(:,k,i)
!
!     Apply image vectors for periodicity
!
!     a-component 
            do while (abs(dist_vec(1)) .gt. 0.5d0) 
               dist_vec(1)=dist_vec(1)-sign(1.d0,dist_vec(1))
            end do
!     b-component                 
            do while (abs(dist_vec(2)) .gt. 0.5d0)
               dist_vec(2)=dist_vec(2)-sign(1.d0,dist_vec(2))
            end do   
!     c-component                 
            do while (abs(dist_vec(3)) .gt. 0.5d0)
               dist_vec(3)=dist_vec(3)-sign(1.d0,dist_vec(3))
            end do
!   
!     Finally, convert distance vector back to cartesian coordinates
!
            dist_vec(:)=dist_vec(1)*cells(1,:,act_conf)+dist_vec(2)*cells(2,:,act_conf)+ &
                       & dist_vec(3)*cells(3,:,act_conf)
!
!     Calculate the distance
!
            dist_act=sqrt(dot_product(dist_vec,dist_vec))            
!
!     Determine if atom is within the cutoff. If yes, add it to the list of 
!      surrounding atoms
!      Fill the respective atom (its coordinates and index) into the local 
!      environment cluster structure
!
            if (dist_act .lt. cutoff) then
               do l=1,nelems_all
                  if (inds(k,i) .eq. ind_list(l)) then                          
                     num_around(l,inc)=num_around(l,inc)+1
                     inc2=inc2+1
                     environ(:,inc2,inc)=-dist_vec(:)
                     ind_env(inc2,inc)=inds(k,i)
                     dist_list(inc2-1)=dist_act
                  end if
               end do
            end if            
         end if
      end do
!
!     Generate radial distribution function for current environment
!      Precompute it for each atom!
!
      do k=1,ngrid
         x_act=k*dx
         do l=2,inc2-1
            rdf_all(k,inc)=rdf_all(k,inc)+exp(-alpha_r*(x_act-dist_list(l))**2)
         end do 
      end do  
!
!     Normalize radial distribution function to number of atoms in environment
!
      rdf_all(:,inc)=rdf_all(:,inc)/(inc2-1)
!
!     Generate angular distribution function for current environment
!     
      angles=0.d0
      inc3=1
      do k=2,inc2-1
         do l=k+1,inc2-1
            angles(inc3)=acos(dot_product(environ(:,k,inc),environ(:,l,inc))/&
                            &sqrt(dist_list(k-1)**2*dist_list(l-1)**2))*180d0/pi 
            inc3=inc3+1        
         end do
      end do   
      do k=l+1,ngrid
         x_act=k*da        
         do l=1,inc3-1
            adf_all(k,inc)=adf_all(k,inc)+exp(-alpha_a*(x_act-angles(l))**2)
         end do
      end do
!
!     Normalize angular distribution function to number of angles in enviroment
!
      adf_all(:,inc)=adf_all(:,inc)/(inc3-1)

   end do
end do

!
!     Write local environments to trajectory file
!
max_around=maxval(sum(num_around,dim=1)) 
open(unit=56,file="environments.xyz",status="replace")
do i=1,1000
   write(56,*) max_around
   write(56,'(a,i8)') " environment No. ",i
   do j=1,max_around
      call atomname(ind_env(j,i),el_act) 
      write(56,*) el_act,environ(:,j,i)
   end do  

end do
close(56)

write(*,*)
write(*,*) "File 'environments.xyz' with local environments written."
write(*,*)

!
!     Test write out for radial distribution functions
!
open(unit=57,file="rdf_test.dat",status="replace")
do i=1,ngrid
   write(57,*) i*dx,rdf_all(i,1:5)
end do
close(57)

!
!     Test write out angular distribution functions
!
open(unit=58,file="adf_test.dat",status="replace")
do i=1,ngrid
   write(58,*) i*da,adf_all(i,1:5)
end do
close(58)



!
!    Now compare the atoms based on their environments!
!    To avoid huge scaling with atom number, do the simplest classifications 
!    first
!
!    A: element
!



!
!    B: Total number of neighbors (no matter which element) 
!
neigh_global=0
do i=1,natoms_sum
   do j=1,nelems
      if (ind_all(i) .eq. ind_list(j)) then
         neigh_global(j,sum(num_around(:,i)))=neigh_global(j, &
                        & sum(num_around(:,i)))+1
      end if
   end do           
end do

open(unit=59,file="neighbor_nums.dat",status="replace")
write(59,*) "#Here, the number of neighbors around atoms in ML_AB are"
write(59,*) "#summed and packed into bins, where the number of atoms with "
write(59,*) "#with an environment of certain size (total number) is listed,"
write(59,*) "#sorted by element."
do i=1,max_environ
   write(59,*) i,real(neigh_global(:,i))

end do

close(59)

!
!    C: Neighborhood diversity (multiply numbers of elements within the current 
!      total neighbor class
!




end program mlff_select

!
!     subroutine elem: read in character with element symbol and
!       give out the number
!
!     part of QMDFF (taken from Caracal)
!
subroutine elem(key1, nat)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CHARACTER(len=*)::KEY1
CHARACTER(len=2)::ELEMNT(107),E

DATA ELEMNT/'h ','he', &
 & 'li','be','b ','c ','n ','o ','f ','ne', &
 & 'na','mg','al','si','p ','s ','cl','ar', &
 & 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu', &
 & 'zn','ga','ge','as','se','br','kr', &
 & 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag', &
 & 'cd','in','sn','sb','te','i ','xe', &
 & 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy', &
 & 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt', &
 & 'au','hg','tl','pb','bi','po','at','rn', &
 & 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','xx', &
 & 'fm','md','cb','xx','xx','xx','xx','xx'/

nat=0
e='  '
do i=1,len(key1)
   if (key1(i:i).ne.' ') L=i
end do
k=1
DO J=1,L
   if (k.gt.2) exit
   N=ICHAR(key1(J:J))
   if (n.ge.ichar('A') .and. n.le.ichar('Z') ) then
      e(k:k)=char(n+ICHAR('a')-ICHAR('A'))
      k=k+1
   end if
   if (n.ge.ichar('a') .and. n.le.ichar('z') ) then
      e(k:k)=key1(j:j)
      k=k+1
   end if
end do

DO I=1,107 
   if (e.eq.elemnt(i)) then
      NAT=I
      RETURN
   END IF
END DO

return
end subroutine elem

!
!     subroutine atomname: is an atomic charge is given, extract the
!        corresponding element label
!
!     part of EVB
!
subroutine atomname(nat,key1)
implicit none
CHARACTER(len=2)::KEY1  ! the name of the element
CHARACTER(len=2)::ELEMNT(107)  ! array with all names
integer::nat  ! atomic charge

DATA ELEMNT/'h ','he', &
 & 'li','be','b ','c ','n ','o ','f ','ne', &
 & 'na','mg','al','si','p ','s ','cl','ar', &
 & 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu', &
 & 'zn','ga','ge','as','se','br','kr', &
 & 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag', &
 & 'cd','in','sn','sb','te','i ','xe', &
 & 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy', &
 & 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt', &
 & 'au','hg','tl','pb','bi','po','at','rn', &
 & 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','xx', &
 & 'fm','md','cb','xx','xx','xx','xx','xx'/

!
!     Return the i'th entry of the elements array of i is the given
!     atomic charge of the atom
!
key1=elemnt(nat)


return
end subroutine atomname



!
!    subroutine matinv3: inversion of 3x3 matrices
!    modified from https://fortranwiki.org/fortran/show/Matrix+inversion
!
subroutine matinv3(A,B)
!! Performs a direct calculation of the inverse of a 3×3 matrix.
real(kind=8), intent(in) :: A(3,3)   !! Matrix
real(kind=8), intent(out)   :: B(3,3)   !! Inverse matrix
real(kind=8)             :: detinv

! Calculate the inverse determinant of the matrix
detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
          - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
          + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

! Calculate the inverse of the matrix
B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))

return
end subroutine
