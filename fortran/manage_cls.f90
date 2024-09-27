!
!    manage_cls: set up and evaluate core level shift (CLS)
!      calculations for structures with several atoms
!    Part of VASP4CLINT
!     Julien Steffen, 2024 (julien.steffen@fau.de)
!
program manage_cls
implicit none 
integer::i,j,k
character(len=2),allocatable::el_list_read(:)
character(len=2),allocatable::el_list(:)
character(len=80)::arg,adum,select_string
character(len=80)::coord_string
character(len=40)::foldername
character(len=120)::a120
character(len=120),allocatable::potcar_block(:)
integer::readstat,el_num
integer::nelems_choice,counter,block_len
logical::cls_elements,el_all
logical::mode_eval,mode_setup
logical::use_slurm
integer,allocatable::list_active(:)
character(len=2),allocatable::el_active(:)
real(kind=8)::cell_scale,a_vec(3),b_vec(3),c_vec(3)
character(len=2),allocatable::el_init(:)
character(len=2),allocatable::el_names(:)
character(len=2),allocatable::at_names(:)
integer,allocatable::el_numbers(:)
integer::natoms,num_active
real(kind=8),allocatable::xyz(:,:)



write(*,*)
write(*,*) "PROGRAM manage_cls: Management of VASP core level shift (CLS)"
write(*,*) " calculations of structures with several atoms."
write(*,*) "In the setup mode, the input files for a new calculation are"
write(*,*) " prepared, in choosing the atoms whose final state energies "
write(*,*) " shall be calculated. A POSCAR file with the structure as well"
write(*,*) " as a INCAR (for final state calc.), KPOINTS and POTCAR file"
write(*,*) " need to be located in the folder."
write(*,*) "In the evaluation mode, the prepared and finished calculations"
write(*,*) " are evaluated and the results collected such that one plot or"
write(*,*) " picture can be produced for the whole structure."
write(*,*) " The evaluation must be called in the same folder as the setup."
write(*,*) "The following command line arguments can/must be given:"
write(*,*) " -setup : chooses the setup mode for a new calculation."
write(*,*) " -eval : chooses the evaluation mode for a done calculation."
write(*,*) " -element=[list of elements]: which element shall be analzed."
write(*,*) "   For each atom of the chosen element in the system, a "
write(*,*) "   separate final state energy calculation will be done!"
write(*,*) "   Choose 'all', if all elements/atoms shall be calculated."
write(*,*) " -slurm : assumes that slurm is used as queue manager. A file"
write(*,*) "   named slurm_script will be copied into each folder."


mode_setup=.false.
mode_eval=.false.
!
!     Determine calculation mode
!
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:6))  .eq. "-setup") then
      mode_setup = .true.
      write(*,*)
      write(*,*) "Mode A was chosen, the setup will be done."
      write(*,*)
   end if
end do

do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:6))  .eq. "-eval") then
      if (mode_setup) then
         write(*,*)
         write(*,*) "Please choose either Mode -eval or -setup, not both!"
         write(*,*)
         stop
      else
         mode_eval = .true.
         write(*,*)
         write(*,*) "Mode B was chosen, the evaluation will be done."
         write(*,*)
      end if
   end if
end do

!
!    Abort if no mode has been chosen
!
if ((.not. mode_setup) .and. (.not. mode_eval)) then
   write(*,*)
   write(*,*) "Please choose either Mode -setup or Mode -eval!"
   write(*,*)
   stop
end if


!
!    The elements whose core level energies shall be calculated
!
nelems_choice=0
cls_elements=.false.
el_all=.false.
allocate(el_list_read(10))
el_list_read="XX"
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:9))  .eq. "-element=") then
      read(arg(10:),*,iostat=readstat) el_list_read
      cls_elements=.true.
   end if
end do

if (cls_elements) then
   nelems_choice=0
   do i=1,10
      if (el_list_read(i) .eq. "XX") exit
!
!    If all elements are chosen (all)
!
      if (el_list_read(i) .eq. "al") then
         el_all = .true.
         exit
      end if
      nelems_choice=nelems_choice+1
   end do

   allocate(el_list(nelems_choice))
   do i=1,nelems_choice
      el_list(i)=el_list_read(i)
   end do
   if (el_all) then
      write(*,*) "All atoms in the system were chosen for the DOS evalulation."
   else
      write(*,*) "The following elements were chosen for the DOS evaluation:"
      do i=1,nelems_choice
         write(*,*) " - ",el_list(i)
      end do
   end if
end if
if (nelems_choice .lt. 1) then
   write(*,*) "Please give at least one element with the -element command!"
   stop
end if        

use_slurm=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:6))  .eq. "-slurm") then
      write(*,*) "The calculations will be started with slurm_scripts."     
      use_slurm=.true.
   end if
end do


!
!    Open the POSCAR file and read in its content
!    This is needed for both modes!
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
do i=1,natoms
   read(56,*) xyz(:,i)
end do
close(56)
!
!    Define individual elements for each atom
!

allocate(at_names(natoms))
counter = 1
do i=1,el_num
   do j=1,el_numbers(i)
      at_names(counter) = el_names(i)
      counter = counter +1
   end do
end do



!
!    MODE A (setup)
!

if (mode_setup) then
   allocate(list_active(natoms))
   allocate(el_active(natoms)) 
!
!    Determine the number and indices of the atoms that 
!    shall be calculated
!
   write(*,*) at_names
   counter=1
   do i=1,nelems_choice
      do j=1,natoms
         if (at_names(j) .eq. el_list(i)) then
            list_active(counter)=j
            el_active(counter)=el_list(i)
            counter=counter+1
         end if
      end do
   end do
   num_active=counter-1

   open(unit=57,file="active_list.dat",status="replace")
   write(57,*) "# This file contains all atoms chosen for final state CLS calculations"
   do i=1,num_active
      write(57,*) list_active(i)
   end do
   write(*,*) "File active_list.dat written."
!
!    Now, setup a calculation folder for each selected atom!
!    In the respective POSCAR, the atom is printed in the last 
!    line
!    Further, a new species is added to the POTCAR file, resembling
!    the element of the atom in the last line
!
   allocate(potcar_block(8000))
   do i=1,num_active
      write(foldername,*) list_active(i)
      call system("mkdir " // adjustl(trim(foldername)))
      call chdir(adjustl(trim(foldername)))
      call system("cp ../KPOINTS .")
      call system("cp ../INCAR .")
      call system("cp ../POTCAR .")
      if (use_slurm) then
         call system("cp ../slurm_script .")
      end if        
      open(unit=45,file="POSCAR",status="replace")
      write(45,*) "Final state calculation for atom ",list_active(i)
      write(45,*) cell_scale
      write(45,*) a_vec(:)
      write(45,*) b_vec(:)
      write(45,*) c_vec(:)
      do j=1,el_num 
         write(45,'(a3,a2)',advance="no") "   ",el_names(j)
      end do
      write(45,'(a3,a2)',advance="no") "   ",el_active(i)
      write(45,*)
      do j=1,el_num
         if (el_names(j) .eq. el_active(i)) then
            write(45,'(a3,i8)',advance="no") "   ",el_numbers(j)-1 
         else
            write(45,'(a3,i8)',advance="no") "   ",el_numbers(j)
         end if        
      end do
      write(45,'(a3,i8)',advance="no") "   ",1
      write(45,*)
      write(45,*) coord_string
!
!     Write the coordinate section, move the line of the active atom
!     to the last line
! 
      do j=1,natoms
         if (j .ne. list_active(i)) then
            write(45,*) xyz(:,j)
         end if        
      end do
      write(45,*) xyz(:,list_active(i))
      close(45)
!
!     Write the entry for the POTCAR file at the last line
!
      block_len=0
      open(unit=67,file="POTCAR",status="old")
      do 
         read(67,'(a)',iostat=readstat) a120
         if (readstat .ne. 0) exit
         if (index(a120,'PAW_PBE '//trim(el_active(i))) .ne. 0) then

            counter=1
            potcar_block(counter)=a120
            do 
               counter=counter+1
               read(67,'(a120)') a120
               potcar_block(counter)=a120
               if (index(a120,'End of Dataset') .ne. 0) exit
            end do      
            block_len=counter
         end if
         if (block_len .ne. 0) exit
      end do
      close(67)
      open(unit=68,file="POTCAR",status="old",access="append")
      do j=1,block_len
         write(68,'(a120)') potcar_block(j)
      end do
      close(68)
!
!     Start the calculation if the slurm_script usage is activated
!
      if (use_slurm) then
         call system("sbatch slurm_script")
         call system("sleep 0.5")
      end if        

      call chdir("..")

   end do
end if        








end program manage_cls        
