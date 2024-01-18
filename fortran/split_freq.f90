!
!    split_freq: Performs VASP frequency calculations for large 
!      systems where the calculation cannot be done within the 
!      limited walltime of a calculation cluster (and unfortunately
!      those calculations cannot be restarted.      
!      Mode A: Takes a given frequency input and divides it into 
!      a desired number of chunks.      
!      Mode B: Collects the output of all calculation chunks and 
!      constructs the total Hessian matrix and from it the frequencies,
!      normal modes and intensities.      
!    Part of VASP4CLINT
!     Julien Steffen, 2024 (julien.steffen@fau.de) 
!      

program split_freq
implicit none 
integer::i,j,k
logical::setup,eval,file_exists
logical,allocatable::frozen_split(:,:)
character(len=120)::arg,adum
character(len=50)::select_string,coord_string,foldername
integer::nchunks,readstat,ind
integer::natoms,el_num,num_moved,workload
character(len=2),allocatable::el_init(:),el_names(:)
character(len=1),allocatable::selective(:,:)
integer,allocatable::el_numbers(:)
real(kind=8)::cell_scale
real(kind=8)::a_vec(3),b_vec(3),c_vec(3)
real(kind=8),allocatable::xyz(:,:)

write(*,*)
write(*,*) "PROGRAM split_freq: performs VASP frequency calculations"
write(*,*) " for large systems where the calculation cannot be done "
write(*,*) " within the limited walltime of a calculation cluster"
write(*,*) " (and unfortunately, these calculations cannot be restated)"
write(*,*) "Two modes are available:"
write(*,*) "-Mode A: Setup: Takes a given frequency input and divides it"
write(*,*) "  into a desired number of chunks. Must be started in a folder"
write(*,*) "  with full VASP input (INCAR, KPOINTS, POTCAR, POSCAR)"
write(*,*) "  Then, a number of subfolders (chunk1, chunk2, ...) is generated"
write(*,*) "  and you must start the calculations in them invidually."
write(*,*) "  usage: split_freq -setup -chunks=[number of chunks]"
write(*,*) "-Mode B: Evaluation: Takes the calculation output within the "
write(*,*) "  chunk folders and sets them together and reconstructs the "
write(*,*) "  global Hessian matrix and frequencies from them. Must be "
write(*,*) "  started in the folder where all chunk folders are located,"
write(*,*) "  together with a POSCAR file, where all moved atoms are 'T T T'"
write(*,*) "  usage: split_freq -eval "

setup=.false.
eval=.false.
!
!     Determine calculation mode
!
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:6))  .eq. "-setup") then
      setup = .true.
      write(*,*)
      write(*,*) "Mode A was chosen, the setup will be done."    
      write(*,*) 
   end if
end do

do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:6))  .eq. "-eval") then
      if (setup) then
         write(*,*)
         write(*,*) "Please choose either Mode A or B, not both!"
         write(*,*)
         stop
      else
         eval = .true.     
         write(*,*)
         write(*,*) "Mode B was chosen, the evaluation will be done."    
         write(*,*) 
      end if   
   end if
end do

!
!    Abort if no mode has been chosen
!
if ((.not. setup) .and. (.not. eval)) then
   write(*,*)     
   write(*,*) "Please choose either Mode A or Mode B!"
   write(*,*)
   stop
end if        
!
!    If the setup mode has been chosen, read in the number of chunks
!
nchunks=0
if (setup) then
   do i = 1, command_argument_count()
      call get_command_argument(i, arg)
      if (trim(arg(1:8))  .eq. "-chunks=") then
         read(arg(9:),*,iostat=readstat) nchunks
         if (readstat .ne. 0) then
            write(*,*) "The format of the -chunks=x command seems to be corrupted!"
            write(*,*)
            stop
         end if     
      end if
   end do
   if (nchunks .lt. 2) then
      write(*,*) "Please give at least two calculation chunks!"
      write(*,*)
      stop
   end if        
end if

!
!    MODE A -------------------
!
if (setup) then
!
!    First, check if all VASP input files are present
!
   inquire(file="INCAR",exist=file_exists) 
   if (.not. file_exists) then
      write(*,*) "Please provide an INCAR file!"
      stop
   end if
   inquire(file="KPOINTS",exist=file_exists)
   if (.not. file_exists) then
      write(*,*) "Please provide a KPOINTS file!"
      stop
   end if
   inquire(file="POSCAR",exist=file_exists)
   if (.not. file_exists) then
      write(*,*) "Please provide a POSCAR file!"
      stop
   end if
   inquire(file="POTCAR",exist=file_exists)
   if (.not. file_exists) then
      write(*,*) "Please provide a POTCAR file!"
      stop
   end if
!
!    Open the POSCAR file and read in its content
!
   open(unit=56,file="POSCAR",status="old",iostat=readstat)
   read(56,*)
   read(56,*) cell_scale
   read(56,*) a_vec(:)
   read(56,*) b_vec(:)
   read(56,*) c_vec(:)
   read(56,'(a)') adum
!
!    Determine number of different elements in file 
!
   do i=10,1,-1
      allocate(el_init(i))
      read(adum,*,iostat=readstat) el_init(:)
      if (readstat .eq. 0) then
         el_num=i
         allocate(el_names(el_num))
         el_names(1:el_num)=el_init(1:el_num)
         exit
      end if
      deallocate(el_init)
   end do

   allocate(el_numbers(el_num))
   read(56,*) el_numbers
   natoms=sum(el_numbers)
   read(56,'(a)') select_string
   read(56,'(a)') coord_string
   allocate(xyz(3,natoms))
   allocate(selective(3,natoms))
   do i=1,natoms
      read(56,*) xyz(:,i),selective(:,i)
   end do
   close(56)

!
!     Determine total number of non-frozen atoms
!     and asssign them to the different chunks of the calculation
!
   allocate(frozen_split(natoms,nchunks+1))
   frozen_split=.false.
   num_moved=0
   do i=1,natoms
      if (selective(1,i) .eq. "T" .or. selective(2,i) .eq. "T" .or. &
                  & selective(3,i) .eq. "T") then      
         frozen_split(i,1) = .true.
         num_moved = num_moved + 1
      end if             
   end do
!
!     Number of atoms per chunk (the last will get less work)
!
   workload=ceiling(real(num_moved)/real(nchunks))
   do i=1,nchunks-1
      do j=1,workload
         do k=1,natoms
            if (frozen_split(k,1)) then
               frozen_split(k,i+1) = .true.
               frozen_split(k,1) = .false.
               exit
            end if        
         end do
      end do
   end do
!
!     The last chunk
!
   do j=1,num_moved-(nchunks-1)*workload
      do k=1,natoms
         if (frozen_split(k,1)) then
            frozen_split(k,nchunks+1) = .true.
            frozen_split(k,1) = .false.
            exit
         end if        
      end do
   end do

!
!     Now write the input for the nchunks partial frequency calculations
!
   do i=1,nchunks
      if (i .lt. 10) then
         write(foldername,'(a,i1)') "chunk",i
      else if (i .lt. 100) then
         write(foldername,'(a,i2)') "chunk",i
      else 
         write(foldername,'(a,i3)') "chunk",i
      end if       
      call system("mkdir "//trim(foldername))
      call system("cp KPOINTS "//trim(foldername))
      call system("cp POTCAR "//trim(foldername))
      call system("cp INCAR "//trim(foldername))
      open(unit=38,file=trim(foldername)//"/POSCAR",status="replace")
      write(38,*) "Frequency chunk calculation written by split_freq"
      write(38,*) cell_scale
      write(38,*) a_vec(:)
      write(38,*) b_vec(:)
      write(38,*) c_vec(:)
      do j=1,el_num
         write(38,'(a,a)',advance="no") "  ",el_names(j) 
      end do
      write(38,*)
      write(38,*) el_numbers
      write(38,*) select_string
      write(38,*) coord_string
      do j=1,natoms
         if (frozen_split(j,i+1)) then
            write(38,'(3f20.13,a)') xyz(:,j),"    T T T "
         else 
            write(38,'(3f20.13,a)') xyz(:,j),"    F F F "
         end if        
      end do
      close(38)
   end do
 

end if        
!
!    MODE B -------------------
!
!    Read in the calculation results from the different chunks
!    and combine them to obtain the full hessian matrix, normal modes
!    and frequencies
!
if (eval) then
!
!    Determine number of chunks: existence of highest folder
!
   ind=0
   do
      if (ind .lt. 10) then
         write(foldername,'(a,i1)') "chunk",ind
      else if (ind .lt. 100) then
         write(foldername,'(a,i2)') "chunk",ind
      else
         write(foldername,'(a,i3)') "chunk",ind
      end if
      inquire(file=foldername//"/vasprun.cml",exist=file_exists)
      if (.not. file_exists) exit
   end do

!
!    First read in the force vectors of all elongations (apart from the unchanged beginning)
!
   do i=1,nchunks 
      if (i .lt. 10) then
         write(foldername,'(a,i1)') "chunk",i
      else if (i .lt. 100) then
         write(foldername,'(a,i2)') "chunk",i
      else
         write(foldername,'(a,i3)') "chunk",i
      end if

      call chdir(foldername)     



   end do




end if        




end program split_freq      
