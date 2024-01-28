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
integer::i,j,k,l,m   ! loop indices
integer::inc,inc2,inc3,inc4  ! incremented index
integer::inc_remain  ! remainder index for basis functions
integer,allocatable::incs(:)  ! moredimensional increment
integer::readstat
integer::nbasis  ! desired number of local refs. in new ML_AB
integer::neigh_min_bas  ! minimum number of confs. per neighbor bin
integer::nbas_grad  ! number of basis functions for max.  gradient comps
integer::neigh_classes ! number of neighborhood subdivision classes
integer::grad_pre    ! number of basis functions for large gradnorms
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
integer::spots_avail  ! number of basis functions for final k-means
integer::remain_class  ! index of the remainder class for basis functions
integer,allocatable::natoms(:)  ! atom numbers of configurations
real(kind=8),allocatable::xyz(:,:,:) ! the geometries
real(kind=8),allocatable::xyz_dir(:,:,:) ! the geometries (direct coords)
real(kind=8),allocatable::cells(:,:,:) ! the unit cells
real(kind=8),allocatable::angles(:)  ! list of angles for current environment
real(kind=8)::grad_min,grad_max  ! minimum and maximum gradient norm (global)
real(kind=8)::grad_step  ! distance between two gradient bins
real(kind=8)::normfac   ! normalizaion factor for RDFs and ADFs
real(kind=8)::rdum   ! dummy real number
integer,allocatable::inds(:,:)  ! the element indices (core charges)
integer::ind_list(50)   ! list of different element indices
character(len=2),allocatable::el_list(:)  ! current list of elements
character(len=2)::el_list_glob(50)  ! global list of elements 
character(len=2)::el_act   ! the current element symbol for printing
character(len=120)::arg  ! command line argument
character(len=20)::bas_scale   ! linear or sqare root scaling of basis classes
real(kind=8)::train_div   ! desired trainset diversity (0: none, 1: full)
real(kind=8)::grad_frac   ! desired fraction of basis functions for large gradnorm
real(kind=8),allocatable::energies(:)  ! the energies 
real(kind=8),allocatable::grads(:,:,:)  ! the gradient vectors
real(kind=8),allocatable::environ(:,:,:)  ! the local environment of each atom
real(kind=8),allocatable::dist_list(:)  ! list of distances in environment
real(kind=8),allocatable::rdf_all(:,:)  ! radial distribution functions for atoms
real(kind=8),allocatable::adf_all(:,:)  ! angular distribution functions for atoms
real(kind=8),allocatable::atsel_grad_val(:,:)  ! gradnorms of largest gradnorms 
real(kind=8),allocatable::rdf_tmp(:,:)  ! temporary RDF array
real(kind=8),allocatable::mat_overlap(:,:)  ! overlap matrix for RDF combinations
real(kind=8)::cell_inv(3,3)  ! inverted coordinate unit cell
real(kind=8)::dist_act   ! scalar distance between two atoms
real(kind=8)::cutoff    ! the distance cutoff during the ML-FF learning
real(kind=8)::dist_vec(3)  ! distance vector between two atoms
real(kind=8)::dx   ! distance between two gridpoints (RDF)
real(kind=8)::da   ! distance between two gridpoints (ADF)
real(kind=8)::alpha_r   ! exponent for RDF Gaussian functions
real(kind=8)::alpha_a   ! exponent for ADF Gaussian functions
real(kind=8)::x_act   ! current x-value for RDF buildup
real(kind=8)::rdf_overlap  ! the current RDF overlap integral
real(kind=8)::pi  ! the Pi
integer,allocatable::el_nums(:)  ! current numbers of elements
integer,allocatable::confnum_all(:)  ! numbers of configurations 
integer,allocatable::nat_all(:) ! number of atom in configurations
integer,allocatable::ind_all(:)  ! element indices of all atoms
integer,allocatable::num_around(:,:) ! number of atoms around the atom
integer,allocatable::ind_env(:,:)  ! core charges of atoms in environments
integer,allocatable::neigh_global(:,:)  ! list with numbers of neighbors (global)
integer,allocatable::neigh_bas(:,:)  ! number of basis functions per environment
integer,allocatable::grad_histo(:)   ! histogram with gradient norm ranges
integer,allocatable::final_choice(:,:)  ! the final selection of basis functions
integer,allocatable::atsel_grad(:,:)  ! the atom indices of the largest gradnorms
integer,allocatable::neighnum_global(:)  ! sum of all neighbors for all atoms
integer,allocatable::neigh_div1(:,:),neigh_div2(:,:)  ! the neighborhood diversities
integer,allocatable::basis_sizes(:)  ! sizes of all basis classes (per element)
integer,allocatable::basis_classes(:,:)  ! list of atoms in each basis class
integer,allocatable::k_number(:)   ! list of k-means spots for basis classes
integer,allocatable::smallest(:)  ! atom list with smallest mutual overlaps
integer,allocatable::hierarchy(:,:)  ! the time-dependent cluster indices
integer::mat_coord(2)  ! the current matrix indices
logical::eval_stat(10)  ! the progress for evaluation loops
logical,allocatable::atom_used(:)  ! boolean mask for blocking of treated atoms
!   for Lapack call (matrix diagonalization)
character(len=1)::JOBZ,UPLO
integer::Nn,LDA,LDU,LDV,LWORK,INFO
real(kind=8),dimension(:),allocatable::WORK,W
real(kind=8),dimension(:,:),allocatable::A_mat


real(kind=8),allocatable::gradnorm_all(:)  ! gradient forms for all atoms

character(len=130)::a130


write(*,*) 
write(*,*) "PROGRAM mlff_select: selection of an effective set of"
write(*,*) " local reference configurations (basis sets) for a VASP"
write(*,*) " machine-learning force field from a given ML-AB file."
write(*,*) "The procedure is similar to the VASP command ML_MODE = select,"
write(*,*) " but is faster (since purely based on geometrical comparisons)"
write(*,*) " and has more options to be adjusted (global approach)"
write(*,*) "The following settings must be given by command line arguments:"
write(*,*) " -nbasis=[number] : the [number] is the desired number of local"
write(*,*) "    reference configurations (basis functions) per element in the"
write(*,*) "    newly written ML_AB file."
write(*,*) " -grad_frac=[value] : percentage of basis functions to be "
write(*,*) "    allocated for the largest gradient components in the "
write(*,*) "    references. The [value]*nbasis atoms with the largest "
write(*,*) "    will be taken all."
write(*,*) " -train_div=[value] : the training set diversity (maximum: 1.0,"
write(*,*) "    minimum: 0.0). The larger the value, the higher percentage"
write(*,*) "    will be chosen based on different number of neighbor atoms,"
write(*,*) "    such that the tails of the neighbor distribution will be "
write(*,*) "    will be weighted larger (larger diversity)"
write(*,*) " -neigh_classes=[number] : Number of diversity-based neighborhood"
write(*,*) "    classes (per number of atoms in the environment, the number"
write(*,*) "    subdivisions, sorted by the diversity (respect to elements)"
write(*,*) "    of the neighborhoods"
write(*,*) " -bas_scale=('linear' or 'root') : How the number of basis functions"
write(*,*) "    allocated to a certain neighborhood class shall scale with"
write(*,*) "    the number of atoms contained into it."
write(*,*) "    possibe are: 'linear' (linear scaling) or 'root' (square root"
write(*,*) "    scaling). The linear scaling will give more spots to classes "
write(*,*) "    with many atoms, whereas the root scaling will favor those"
write(*,*) "    with lower number of atoms."
write(*,*) "Usage: mlff_select [list of command line arguments]"
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
!     The fraction of basis functions based on the largest gradient norms
!
grad_frac=-1.d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:11))  .eq. "-grad_frac=") then
      read(arg(12:),*,iostat=readstat) grad_frac
      if (readstat .ne. 0) then
         write(*,*) "The format of the -grad_frac=[value] command seems to be corrupted!"
         write(*,*)
         stop
      end if
   end if
end do
if (grad_frac .lt. 0.0 .or. grad_frac .gt. 1.0) then
   write(*,*) "Please give the desired fraction of basis functions for the largest "
   write(*,*) " gradient norm atoms with the keyword -grad_frac=[value]"
   write(*,*) " where the maximum fraction is 1.0 (all) and the minimum is 0.0 (none)"
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
      read(arg(12:),*,iostat=readstat) train_div
      if (readstat .ne. 0) then
         write(*,*) "The format of the -train_div=[value] command seems to be corrupted!"
         write(*,*)
         stop
      end if
   end if
end do
if (train_div .lt. 0.0 .or. train_div .gt. 1.0) then
   write(*,*) "Please give the desired trainset diversity with the keyword -train_div=[value]"
   write(*,*) " where the maximum diversity is 1.0 and the minimum is 0.0"
   write(*,*)
   stop
end if

!
!     The number of neighborhood classes
!
neigh_classes=0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:15))  .eq. "-neigh_classes=") then
      read(arg(16:),*,iostat=readstat) neigh_classes
      if (readstat .ne. 0) then
         write(*,*) "The format of the -neigh_classes=[number] command seems to be corrupted!"
         write(*,*)
         stop
      end if
   end if
end do
if (neigh_classes .lt. 1) then
   write(*,*) "Please give the number of neighborhood classes based on their element "
   write(*,*) " diversities with the keyword -neigh_classes=[number], where the number"
   write(*,*) " of classes must be 1 (no subdividing) or larger."
   write(*,*)
   stop
end if

!
!     If the number of basis functions for the neighborhood classes shall be scaling
!      linear or square root with the number of contained atoms
!
bas_scale="none"
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:11))  .eq. "-bas_scale=") then
      read(arg(12:),*,iostat=readstat) bas_scale
      if (readstat .ne. 0) then
         write(*,*) "The format of the -bas_scale=(linear or root) command seems to be corrupted!"
         write(*,*)
         stop
      end if
   end if
end do
if ((bas_scale .ne. "linear") .and. (bas_scale .ne. "root")) then
   write(*,*) "Please decide if the number of basis functions per neighborhood class shall "
   write(*,*) " scale linearly or square root wise with the number of contained atoms."
   write(*,*) " Linear scaling: -bas_scale=linear"
   write(*,*) " Square root scaling: -bas_scale=root"
   write(*,*)
   stop
end if

!
!     Number of basis functions allocated to minimum into different neighbor bins 
!
neigh_min_bas=int(nbasis*train_div)


write(*,'(a,i7)') " Number of desired local reference configuraitons per element:",nbasis
write(*,'(a,i7)') " Number of configurations allocated uniformly for neighbor-bins:",neigh_min_bas

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
open(unit=56,file="ML_AB",status="old",iostat=readstat)
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
read(56,*,iostat=readstat)
if (readstat .ne. 0) then
   write(*,*) "The file ML_AB seems to be empty!"
   stop
end if
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
!     The nelems-dimensional increment array
allocate(incs(nelems))
!     The configuration number of the element
allocate(confnum_all(natoms_sum))
!     Boolean mask, if atom was already allocated to final basis
allocate(atom_used(natoms_sum))
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
!     Number of basis functions per environment
allocate(neigh_bas(nelems,max_environ))
!     Number of atoms with gradient components in certain range
allocate(grad_histo(100))
!     The sum of all neighbors around all atoms
allocate(neighnum_global(natoms_sum))


!     The FINAL choice of basis functions!
allocate(final_choice(nbasis,nelems))





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
!      Precompute it for each atom and normalize its integral!
!
      normfac=0.d0
      do k=1,ngrid
         x_act=k*dx
         do l=2,inc2-1
            rdf_all(k,inc)=rdf_all(k,inc)+exp(-alpha_r*(x_act-dist_list(l))**2)
         end do
         normfac=normfac+rdf_all(k,inc)*rdf_all(k,inc)*dx
      end do 
!
!     Normalize its integral to become 1
!
      rdf_all(:,inc)=rdf_all(:,inc)/sqrt(normfac)

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
      normfac=0.d0 
      do k=1,ngrid
         x_act=k*da        
         do l=1,inc3-1
            adf_all(k,inc)=adf_all(k,inc)+exp(-alpha_a*(x_act-angles(l))**2)
         end do
         normfac=normfac+adf_all(k,inc)*adf_all(k,inc)*da
      end do
!
!     Normalize its integral to become 1
!
      adf_all(:,inc)=adf_all(:,inc)/sqrt(normfac)

   end do
end do
!
!     All atoms are still available
!
atom_used=.true.
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
!    A: GRADIENT EXTREMA PRESELECTION:
!    Determine histogram of gradient norms for all atoms, allocate them into
!    bins according to their gradient norms
!    Choose the nbasis*grad_frac atoms with the largest gradient norms
!
!    The number of atoms preselected due to large gradient components
!
grad_pre=int(nbasis*grad_frac)
allocate(atsel_grad(grad_pre,nelems))
allocate(atsel_grad_val(grad_pre,nelems))


grad_min=minval(gradnorm_all)
grad_max=maxval(gradnorm_all)
grad_step=(grad_max-grad_min)/100.d0
grad_histo=0

outer2: do i=1,natoms_sum
   inner2: do j=1,100
      if (gradnorm_all(i) .lt. j*grad_step) then
         grad_histo(j)=grad_histo(j)+1
         exit inner2
      end if
   end do inner2
end do outer2
!
!    Write distributions of gradient norms to file
!
open(unit=60,file="grad_histo.dat",status="replace")
write(60,*) "#Here, the number of atoms with gradient norms in a certain"
write(60,*) "#range (x-coordinate) are listed."
do i=1,100
   write(60,*) i*grad_step,grad_histo(i)

end do
close(60)

!
!    Sort the grad_pre atoms with the largest gradients into an array with
!    their atom indices, do a max-heap search
!  
!    Stored in array: atsel_grad
!

incs=0
atsel_grad=0
atsel_grad_val=0.d0
do
   inc=maxloc(gradnorm_all,dim=1)
   do i=1,nelems
      if (ind_all(inc) .eq. ind_list(i)) then
         if (incs(i) .lt. grad_pre) then
            incs(i)=incs(i)+1
            atsel_grad(incs(i),i)=inc
            atom_used(inc)=.false. 
!    Store atom indices in final basis set!
            final_choice(incs(i),i)=inc
            atsel_grad_val(incs(i),i)=gradnorm_all(inc)
         end if
      end if
   end do
   if (sum(incs) .eq. nelems*grad_pre) then
      exit
   end if
   gradnorm_all(inc)=0.d0
end do
!    
!    B: THE LOCAL NEIGHBORS: number and diversity, minimum number!
!
neigh_global=0
neighnum_global=0
do i=1,natoms_sum
   if (atom_used(i)) then
      do j=1,nelems
         if (ind_all(i) .eq. ind_list(j)) then
            inc=sum(num_around(:,i))
            neigh_global(j,inc)=neigh_global(j,inc)+1
            neighnum_global(i)=neighnum_global(i)+inc
         end if
      end do     
   end if      
end do

!
!    Fill diversity-dependent minimum number of allocated basis functions
!      into different neighbor bins (each bin one, until full or all atoms
!      allocated, line per line)
!

neigh_bas=0
do_elems: do i=1,nelems
   inc=0
   do_column:  do
      do_environ: do j=1,max_environ
         inc=inc+1
         if (inc .gt. neigh_min_bas) exit do_column
         if (neigh_bas(i,j) .lt. neigh_global(i,j)) then
            neigh_bas(i,j)=neigh_bas(i,j)+1
            do_atoms: do k=1,natoms_sum
!
!      Allocate chosen atom to final basis set
!
               if (atom_used(k)) then
                  if (neighnum_global(k) .eq. j) then
                     atom_used(k) = .false.
                     final_choice(inc,i)=k
                     exit do_atoms
                  end if
               end if
            end do do_atoms
         end if
      end do do_environ
   end do do_column
end do do_elems

!
!    Write distribution of neighbor numbers to file
!
open(unit=59,file="neighbor_nums.dat",status="replace")
write(59,*) "#Here, the number of neighbors around atoms in ML_AB are"
write(59,*) "#summed and packed into bins, where the number of atoms with "
write(59,*) "#with an environment of certain size (total number) is listed,"
write(59,*) "#sorted by element."
do i=1,max_environ
   write(59,*) i,real(neigh_global(:,i)),real(neigh_bas(:,i))

end do
close(59)


!    
!     C: THE LOCAL NEIGHBORS: number and diversity, final allocation via 
!      k-means cluster analysis!
!      Outer loop over elements: separate calculation for each element!
!
do i=1,nelems
!
!     Number of available basis functions:
!
   write(*,*) "available:",grad_pre,neigh_min_bas
   inc=0
   inc3=0
   do j=1,max_environ
      if ((neigh_global(i,j)-neigh_bas(i,j)) .gt. 0) then
         inc=inc+1
      end if
   end do
!
!    Allocate array with atom classes based on neighborhood size 
!     and diversity     
!
   if (allocated(basis_sizes)) deallocate(basis_sizes)
   if (allocated(basis_classes)) deallocate(basis_classes)
   if (allocated(neigh_div1)) deallocate(neigh_div1)
   if (allocated(neigh_div2)) deallocate(neigh_div2)
!
!     Add an extra dimension (+1) as remainder for all classes with less 
!     then 10 atoms
!
   allocate(basis_sizes(inc*neigh_classes+1))
   allocate(basis_classes(maxval(neigh_global-neigh_bas),inc*neigh_classes+1))
   allocate(neigh_div1(2,maxval(neigh_global-neigh_bas)))
   allocate(neigh_div2(2,maxval(neigh_global-neigh_bas)))
   remain_class=inc*neigh_classes+1
   inc_remain=0
   basis_sizes=0
   basis_classes=0
!
!    Subdivide the number of neighbor bins into regions based on neighborhood
!    diversity (how many of the different elements are within the neighborhood)
!
!    Fill local array with product of all numbers of elements in the surrounding
!
   inc=0
   do j=1,max_environ
      neigh_div1=1
      neigh_div2=1
      inc=0
      do k=1,natoms_sum
         if (atom_used(k) .and. (ind_all(k) .eq. ind_list(i))) then
            if (neighnum_global(k) .eq. j) then
               inc=inc+1
               neigh_div1(1,inc)=k
               do l=1,nelems
                  if (num_around(l,k) .gt. 0) then
                     neigh_div1(2,inc)=neigh_div1(2,inc)*num_around(l,k)
                  end if
               end do
            end if     
         end if 
      end do
!
!    Second, sort the array with the surroundings
!
      do k=1,inc
         inc2=maxloc(neigh_div1(2,1:inc),dim=1)
         neigh_div2(:,k)=neigh_div1(:,inc2)
         neigh_div1(:,inc2)=1 
      end do
!
!    Third, allocate the atoms into the neighborhood class, and with this 
!    into the final cluster for the k-means clustering!
! 
      if (inc .gt. 0) then
         inc3=inc3+1
         inc4=0
         do k=1,neigh_classes-1
            if ((inc/neigh_classes) .lt. 10) then
               basis_sizes(remain_class)=basis_sizes(remain_class)+inc/neigh_classes
            else
               basis_sizes((inc3-1)*neigh_classes+k)=inc/neigh_classes  
            end if
            if ((inc/neigh_classes) .lt. 10) then
               do l=1,inc/neigh_classes
                  inc4=inc4+1
                  inc_remain=inc_remain+1
                  basis_classes(inc_remain,remain_class)=neigh_div2(1,inc4)
               end do
            else 
               do l=1,inc/neigh_classes
                  inc4=inc4+1
                  basis_classes(l,(inc3-1)*neigh_classes+k)=neigh_div2(1,inc4)
               end do 
            end if
         end do  
         if ((inc-inc4) .lt. 10) then 
            basis_sizes(remain_class)=basis_sizes(remain_class)+inc-inc4
         else
            basis_sizes(inc3*neigh_classes)=inc-inc4
         end if
         if ((inc-inc4) .lt. 10) then
            do l=1,inc-inc4
               inc4=inc4+1
               inc_remain=inc_remain+1
               basis_classes(inc_remain,remain_class)=neigh_div2(1,inc4)
            end do            
         else
            do l=1,basis_sizes(inc3*neigh_classes)
               inc4=inc4+1
               basis_classes(l,inc3*neigh_classes)=neigh_div2(1,inc4)
            end do
         end if
      end if
   end do
!
!    Fourth, determine the number of available spots in the final basis set
!     for the different local environment diversity classes  
!
!     Total number of available spots
!
   spots_avail=nbasis-grad_pre-neigh_min_bas
!
!     Array with number of spots per neighbor subclass
!
   allocate(k_number(size(basis_sizes)))
!
   do j=1,size(basis_sizes)
      if (bas_scale .eq. "linear") then
         inc=inc+basis_sizes(j)
      else if (bas_scale .eq. "root") then
         inc=inc+int(sqrt(real(basis_sizes(j))))
      end if
   end do

!
!     Determine the numbers
!   
   do j=1,size(basis_sizes)
!
!     Linear scaling with number of atoms
!
      if (bas_scale .eq. "linear") then
         k_number(j)=nint(real(basis_sizes(j))/real(inc)*real(spots_avail))

!
!     Square root scaling with number of atoms
!
      else if (bas_scale .eq. "root") then
         k_number(j)=nint(sqrt(real(basis_sizes(j)))/real(inc)*real(spots_avail))

      end if
   end do
!
!     If the tota number of basis functions is too high or too low, increase 
!     or decrease them until the exact number is reached
!
   if (sum(k_number(1:size(basis_sizes))) .gt. spots_avail) then
      inc_remain=sum(k_number(1:size(basis_sizes)))-spots_avail
      do j=1,size(basis_sizes)
         if (k_number(j) .gt. 1) then
            k_number(j) = k_number(j)-1
            inc_remain=inc_remain-1
            if (inc_remain .lt. 1) exit
         end if
      end do
   end if
   if (sum(k_number(1:size(basis_sizes))) .lt. spots_avail) then
      inc_remain=spots_avail-sum(k_number(1:size(basis_sizes)))
      do j=1,size(basis_sizes)
         if (k_number(j) .ge. 1) then
            k_number(j) = k_number(j)+1
            inc_remain=inc_remain+1
            if (inc_remain .gt. -1) exit
         end if
      end do
   end if
!
!     Now perform the k-means clustering for each diversity-neighbor 
!     number group and obtain final choice of atoms for basis set
!
   do j=1,size(basis_sizes)
      if (allocated(mat_overlap)) deallocate(mat_overlap)
      allocate(mat_overlap(basis_sizes(j),basis_sizes(j)))
      if (allocated(hierarchy)) deallocate(hierarchy)
      allocate(hierarchy(basis_sizes(j),basis_sizes(j)))
      mat_overlap=1.1d0
      if (k_number(j) .gt. 0) then
         if (j .eq. 13) write(*,*) k_number(j)
         do k=1,basis_sizes(j)
            do l=1,k-1
!
!     Calculate radial distribution function (RDF) and angular distribution
!      function (ADF) overlaps of all atoms in the cluster
!      Calculate 1-overlap, in order to get the lowest values for the 
!      largest overlaps (diagonal: zero)
!
               rdf_overlap=0.d0
               do m=1,ngrid
                  rdf_overlap=rdf_overlap+rdf_all(m,basis_classes(k,j))* &
                              & rdf_all(m,basis_classes(l,j))*dx
               end do
               mat_overlap(k,l)=1.d0-rdf_overlap
            end do
         end do 
!
!     Perform the hierarchial clustering of the current RDF overlap matrix
!
!     The initial hierarchy: all indices as usual
!
         do l=1,basis_sizes(j)
            hierarchy(l,1)=l
         end do
         write(*,*) hierarchy(:,1)
         inc2=basis_sizes(j)+1
         do k=2,basis_sizes(j)
            do l=1,basis_sizes(j)
               hierarchy(l,k)=hierarchy(l,k-1)
            end do
!
!     Determine indices with lowest current value and update cluster list
!
            mat_coord=minloc(mat_overlap) 
            hierarchy(mat_coord(1),k)=inc2
            hierarchy(mat_coord(2),k)=inc2
            rdum=mat_overlap(mat_coord(1),mat_coord(2))
!
!     If the current index is already a cluster, set all of its atoms to the 
!      new cluster
!
            do l=1,basis_sizes(j)
               if (hierarchy(l,k) .eq. hierarchy(mat_coord(1),k-1)) then
                  hierarchy(l,k) = inc2
               end if
               if (hierarchy(l,k) .eq. hierarchy(mat_coord(2),k-1)) then
                  hierarchy(l,k) = inc2
               end if
            end do
      
!
!     Update all elements in the row and columns of the minloc values to 
!     the larger values of the respective column/row
!
            do l=1,basis_sizes(j)
               if ((mat_overlap(l,mat_coord(1)) .lt. 1.d0) .and. & 
                    & (mat_overlap(l,mat_coord(2)) .lt. 1.d0)) then
                  if (mat_overlap(l,mat_coord(1)) .gt. mat_overlap(l,mat_coord(2))) then
                     mat_overlap(l,mat_coord(2)) = mat_overlap(l,mat_coord(1))
                  else
                     mat_overlap(l,mat_coord(1)) = mat_overlap(l,mat_coord(2))
                  end if 
               end if
            end do

            do l=1,basis_sizes(j)
               if ((mat_overlap(mat_coord(1),l) .lt. 1.d0) .and. &
                    & (mat_overlap(mat_coord(2),l) .lt. 1.d0)) then
                  if (mat_overlap(mat_coord(1),l) .gt. mat_overlap(mat_coord(2),l)) then
                     mat_overlap(mat_coord(2),l) = mat_overlap(mat_coord(1),l)
                  else
                     mat_overlap(mat_coord(1),l) = mat_overlap(mat_coord(2),l)
                  end if
               end if
            end do
!
!     Set all values that appear twice or more in the merged cluster to 1.0
!
            inc3=1
            do l=1,basis_sizes(j)
               if (hierarchy(l,k) .eq. inc2) then
                  if (inc3 .gt. 1) then
                     mat_overlap(:,l) = 1.d0
                  end if
                  inc3=inc3+1
               end if
            end do
            inc3=1
            do l=1,basis_sizes(j)
               if (hierarchy(l,k) .eq. inc2) then
                  if (inc3 .gt. 1) then
                     mat_overlap(l,:) = 1.d0
                  end if
                  inc3=inc3+1
               end if
            end do

            do l=1,basis_sizes(j)
            write(199,'(12f10.6)') mat_overlap(l,:)
            end do
            write(199,*)

            inc2=inc2+1       
            write(*,'(12i5,f10.6)') hierarchy(:,k),rdum
         end do
!
!     Diagonalize overlap matrix
!
!         if (allocated(A_mat)) deallocate(A_mat)
!         if (allocated(W)) deallocate(W)
!         if (allocated(WORK)) deallocate(WORK)
!         JOBZ='V' !eigenvalues and eigenvectors(U)
!         UPLO='U' !upper triangle of a
!         Nn=basis_sizes(j)
!         LDA=Nn
!         INFO=0
!         LWORK=Nn*Nn-1
!         allocate(A_mat(Nn,Nn))
!         allocate(W(Nn))
!         allocate(WORK(LWORK))
!         A_mat=mat_overlap
!         call DSYEV(JOBZ,UPLO,Nn,A_mat,LDA,W,WORK,LWORK,INFO)
!
!     Aanalyze eigenvalues and eigenvectors
!
!         do k=1,basis_sizes(j)
!            write(*,*) "value",k,W(k)
!            do l=1,basis_sizes(j)
!               write(*,*) l,abs(A_mat(k,l))
!            end do
!         end do
!         A_mat=abs(A_mat)
!
!     Calculate overlap sum of largest eigenvector components
!
!         if (allocated(smallest)) deallocate(smallest)
!         allocate(smallest(k_number(k)))         
!         do k=1,basis_sizes(j)
!            do l=1,k_number(j)
!               smallest(l)=maxloc(A_mat(:,k),dim=1) 
!               A_mat(smallest(l),k)=0.0
!            end do
!            write(*,*) k,smallest
!            rdum=0.d0
!            do l=1,k_number(j)
!               do m=1,k_number(j)
!                  rdum=rdum+mat_overlap(smallest(l),smallest(m))
!               end do
!            end do
!            write(*,*) "summ",rdum
!
!         end do
      end if
      if (j .eq. 13) stop "GHüuipg"
      write(*,*) j,basis_sizes(j)
   end do


   stop "Hiühi"
end do





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
