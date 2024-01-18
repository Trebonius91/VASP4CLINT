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
integer::i,j,k,l
logical::setup,eval,file_exists
logical,allocatable::frozen_split(:,:),frozen_all(:)
character(len=120)::arg,adum,line
character(len=50)::select_string,coord_string,foldername,string
integer::nchunks,readstat,ind
integer::natoms,el_num,num_moved,workload
integer::force_act,hit_act
character(len=2),allocatable::el_init(:),el_names(:),at_names(:)
character(len=1),allocatable::selective(:,:)
integer,allocatable::el_numbers(:),chunk_lens(:),at_moved(:)
integer,allocatable::at_index(:)
real(kind=8)::cell_scale
real(kind=8)::potim
real(kind=8)::a_vec(3),b_vec(3),c_vec(3)
real(kind=8),allocatable::xyz(:,:),force_vecs(:,:,:),hess(:,:)
real(kind=8),allocatable::el_mass(:),coord_mass(:),at_mass(:)
real(kind=8),allocatable::n_vib(:,:,:)
character(len=1)::JOBZ,UPLO
integer::M,Nn,LDA,LDU,LDV,LWORK,INFO
real(kind=8),dimension(:),allocatable::WORK,W
real(kind=8),dimension(:,:),allocatable::A_mat
real(kind=8),dimension(:,:),allocatable::U

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
   if (readstat .ne. 0) then
      write(*,*) "The POSCAR file could not been found!"
      stop
   end if        
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
   ind=1
   do
      if (ind .lt. 10) then
         write(foldername,'(a,i1)') "chunk",ind
      else if (ind .lt. 100) then
         write(foldername,'(a,i2)') "chunk",ind
      else
         write(foldername,'(a,i3)') "chunk",ind
      end if
      inquire(file=trim(foldername)//"/vasprun.xml",exist=file_exists)
      if (.not. file_exists) exit
      ind=ind+1
   end do
   nchunks=ind-1
   if (nchunks .lt. 2) then
      write(*,*) "There need to be at least two chunk folders with vasprun.xml files!"
      stop 
   end if  
   write(*,'(a,i4,a)') " In total, ",nchunks," frequency calculation chunks were performed." 
!
!    Open the POSCAR file to determine the total number of atoms and which atoms were 
!    moved during the frequency calculation
!
      
   open(unit=56,file="POSCAR",status="old",iostat=readstat)
   if (readstat .ne. 0) then
      write(*,*) "The POSCAR file could not been found!"
      stop
   end if        
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
   allocate(el_mass(el_num))
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
!     Convert coordinates to cartesian coordinates if needed
!
   if (index(coord_string,'Direct') .ne. 0 .or. index(coord_string,'direct') .ne. 0) then
      do i=1,natoms
         xyz(:,i)=xyz(1,i)*a_vec(:)+xyz(2,i)*b_vec(:)+xyz(3,i)*c_vec(:)
      end do
   end if        

!
!     Determine total number of non-frozen atoms
!     and asssign them to the different chunks of the calculation
!
   allocate(chunk_lens(nchunks))
   allocate(frozen_split(natoms,nchunks+1))
   allocate(frozen_all(natoms))
   frozen_split=.false.
   frozen_all=.false.
   num_moved=0
   do i=1,natoms
      if (selective(1,i) .eq. "T" .or. selective(2,i) .eq. "T" .or. &
                  & selective(3,i) .eq. "T") then
         frozen_split(i,1) = .true.
         frozen_all(i)=.true.
         num_moved = num_moved + 1
      end if
   end do
   allocate(at_moved(num_moved))
   ind=1
   do i=1,natoms
      if (frozen_split(i,1)) then
         at_moved(ind) = i
         ind=ind+1
      end if        
   end do
   write(*,'(a,i8,a)') " The system contains ",natoms," atoms."
   write(*,'(a,i8,a)') " Of them, ",num_moved," atoms were activated for numerical frequencies."
!
!     Number of atoms per chunk (the last will get less work)
!
   workload=ceiling(real(num_moved)/real(nchunks))
   do i=1,nchunks-1
      chunk_lens(i)=workload
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
   chunk_lens(nchunks)=num_moved-(nchunks-1)*workload
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
!    First read in the force vectors of all elongations (apart from the unchanged beginning)
!
   force_act=0
   allocate(force_vecs(3,natoms,6*num_moved))
   allocate(hess(3*num_moved,3*num_moved))

   write(*,*)
   write(*,*) "Set all chunks together and construct the Hessian matrix..."

   do i=1,nchunks 
      if (i .lt. 10) then
         write(foldername,'(a,i1)') "chunk",i
      else if (i .lt. 100) then
         write(foldername,'(a,i2)') "chunk",i
      else
         write(foldername,'(a,i3)') "chunk",i
      end if

      call chdir(foldername)   
!
!    Open and read in the vasprun.xml file
!
      open(unit=78,file="vasprun.xml",iostat=readstat)
      if (readstat .ne. 0) then
         write(*,*) "The vasprun.xml file in folder ",trim(foldername)," could not been found!"
         stop
      end if 
      hit_act=0      
      do
         read(78,'(a)',iostat=readstat) line
         if (readstat .ne. 0) exit
         if (index(line,'forces') .ne. 0) then
            hit_act=hit_act+1
            if (hit_act .gt. 1) then     
               force_act=force_act+1    
               do j=1,natoms
                  read(78,*) adum,force_vecs(:,j,force_act)
               end do
            end if   
         end if   
!
!    The elongation during the numerical calculation
!
         if (index(line,'POTIM') .ne. 0) then
            read(line,*) adum,adum,string
            read(string(1:8),*) potim 
         end if        
!
!    The masses of all atoms 
!
         if (index(line,'atomspertype') .ne. 0) then
            if (el_mass(1) .lt. 0.01) then         
               read(78,*)
               read(78,*)
               read(78,*)
               read(78,*)
               read(78,*)
               do k=1,el_num
                  read(78,'(a)') line
                  read(line(34:46),*) el_mass(k)
               end do
               allocate(coord_mass(3*num_moved),at_mass(natoms))
               allocate(at_names(natoms))
               allocate(at_index(natoms))
               ind=1
               do k=1,el_num
                  do l=1,el_numbers(k)
                     at_mass(ind)=el_mass(k)
                     at_names(ind)=el_names(k)
                     call elem(at_names(ind),at_index(ind))
                     ind=ind+1
                  end do                  
               end do

               do k=1,num_moved
                  coord_mass((k-1)*3+1)=at_mass(at_moved(k))
                  coord_mass((k-1)*3+2)=at_mass(at_moved(k))
                  coord_mass((k-1)*3+3)=at_mass(at_moved(k))
               end do
            end if        
         end if


      end do
    !  force_act=force_act-1
      close(78)

      call chdir("..")
      
   end do

   do i=1,3*num_moved
   !   write(*,*) "mass",i,coord_mass(i)
   end do
!
!    Now calculate the full Hessian numerically
!
   do i=1,3*num_moved
      do j=1,num_moved
         do k=1,3
            hess(i,(j-1)*3+k) = -(force_vecs(k,at_moved(j),(i-1)*2+1)- &
                     & force_vecs(k,at_moved(j),(i-1)*2+2))/(2d0*potim)
         end do
      end do
   end do
!
!    Mass-scale the hessian for correct frequencies
!
   do i=1,3*num_moved
      do j=1,3*num_moved
         hess(i,j)=hess(i,j)/sqrt(coord_mass(i)*coord_mass(j))
      end do
   end do
   write(*,*) "done!"
   write(*,*) "Diagonalize the Hessian and calculate normal modes..."
!
!    Diagonalize the hessian
!
!   therafter: W=frequencies, A=normal coordinates
!
   JOBZ='V' !eigenvalues and eigenvectors(U)
   UPLO='U' !upper triangle of a
   Nn=3*num_moved
   LDA=Nn
   INFO=0
   LWORK=Nn*Nn-1
   allocate(A_mat(Nn,Nn))
   allocate(W(Nn))
   allocate(WORK(LWORK))
   A_mat=hess
   call DSYEV(JOBZ,UPLO,Nn,A_mat,LDA,W,WORK,LWORK,INFO)

   if (info .ne. 0) then
      write(*,*) "The diagonalization of the hessian matrix with Lapack returned"
      write(*,*) " an error code!"
      stop
   end if
!
!     Obtain frequencies with correct unit etc.
!
   do i = 1, 3*num_moved
      W(i)=sign(sqrt(abs(W(i))),W(i))
   end do
!
!     Do conversion according to Mopac manual: 
!     http://openmopac.net/manual/Hessian_Matrix.html
!
!     1: from eV/Ang^2 to kcal/(mol*Ang^2)
!
   W=W*sqrt(23.0609d0)
!
!     2: from kcal/(mol*Ang^2) to millidynes/Ang (and to Newton/m/kg)
!
   W=W*sqrt(1E8*4184d0/(1E-10*6.023E23))
!
!     3: from millidynes/Ang to cm^-1
!
   W=W*1302.79d0
!
!     Calculate normal mode vibration vectors
!
   allocate(n_vib(3,num_moved,Nn))
   do i=1,3*num_moved
      do j=1,num_moved
         do k=1,3
            n_vib(k,j,i)=A_mat((j-1)*3+k,i)/sqrt(coord_mass((j-1)*3+k))
         end do
      end do
   end do
   write(*,*) "done!"

!
!     Write molden format output file
!
   open(unit=49,file="combined_full.molden",status="replace")
   write(49,*) "[Molden Format]"
   write(49,*) "[Atoms] Angs"
   do i=1,natoms
      write(49,'(a,a,i7,i4,3f12.6)') " ",at_names(i),i,at_index(i),xyz(:,i)
   end do
   write(49,*) "[FREQ]"
   do i=1,3*num_moved
      write(49,'(f10.4)') W(i)
   end do
   write(49,*) "[INT]"
   do i=1,3*num_moved
      write(49,*) 1.0d0
   end do
   write(49,*) "[FR-COORD]"
   do i=1,natoms
      write(49,'(a,a,3f11.6)') " ",at_names(i),xyz(:,i)/0.52917721092d0
   end do
   write(49,*) "[FR-NORM-COORD]"
   do i=1,3*num_moved
      write(49,'(a,i7)') "vibration",i
      ind=1
      do j=1,natoms
         if (frozen_all(j)) then
            write(49,'(3f11.6)') n_vib(:,ind,i)
            ind=ind+1
         else
            write(49,'(3f11.6)') 0.0d0, 0.0d0, 0.0d0  
         end if    
      end do
   end do
   close(49)
   write(*,*)
   write(*,*) "File 'combined_full.molden' for mode visualization written!"
end if        



write(*,*)
write(*,*) "Program split_freq exited normally."
write(*,*)

end program split_freq    


!
!     subroutine elem: read in character with element symbol and
!       give out the number
!
!     part of QMDFF
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

