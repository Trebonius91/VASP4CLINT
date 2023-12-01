!
!    modify_xdatcar: Perform different tasks on a XDATCAR file,
!      such as translating or multiplying the unit cell, picking
!      certain structures or write an xyz trajectory file
!    Part of VASP4CLINT
!     Julien Steffen, 2023 (julien.steffen@fau.de)
!


program modify_xdatcar
implicit none 
integer::i,j,k
integer::readstat,openstat,counter,endl
integer::natoms,nelems,xdat_lines,nframes
real(kind=8)::a_vec(3),b_vec(3),c_vec(3)
character(len=2),allocatable::el_names_read(:),el_names(:)
character(len=2),allocatable::at_names(:)
integer,allocatable::el_nums(:)
real(kind=8)::act_num(3),factor
real(kind=8),allocatable::xyz(:,:,:)
real(kind=8)::xlen,ylen,zlen
real(kind=8)::shift_vec(3)
integer::multiply_vec(3)
logical::eval_stat(10)
logical::shift_cell,multiply_cell
character(len=120)::a120,cdum,arg
character(len=220)::a220
character(len=50)::atest

write(*,*) "PROGRAM modify_xdatcar: Modification of XDATCAR trajectory"
write(*,*) " files for arbitrary systems."
write(*,*) "The file XDATCAR must be present."
write(*,*) "The following calculations can be done, called by one of the"
write(*,*) " listed command line arguments:"
write(*,*) " -shift=a,b,c : Shift each frame of the trajectory by a vector"
write(*,*) "    given in direct coordinates. Example: -shift=0.1,0.0,0.2"
write(*,*) " -multiply=a,b,c : Multiply the unit cell of each frame by some"
write(*,*) "    replications in each of the coordinate directions. Integers"
write(*,*) "    must be given as arguments. Example: -multiply=2,2,1"
write(*,*) " -pick_frame=num : Pick one frame of the XDATCAR file and print"
write(*,*) "    it to a POSCAR file (POSCAR_pick). Example: -pick_frame=283"
write(*,*) " -print_xyz : Print each frame to a xyz trajectory: xdatcat.xyz"
write(*,*) "One or several of these commands can be executed. The ordering"
write(*,*) " will be the same as the listing of the commands."
!
!    PART A: Read in the command line arguments !!!!!!!!!!!!!!!!
!
!    The coordinate shift vector
!
shift_vec=0.0d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:7))  .eq. "-shift=") then
      read(arg(18:),*,iostat=readstat) shift_vec
      shift_cell=.true.
      if (readstat .ne. 0) then
         stop "Check the command -shift=..., something went wrong!"
      end if        
   end if
end do
!
!    The multiplication of unit cells
!
multiply_vec=1
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:10))  .eq. "-multiply=") then
      read(arg(18:),*,iostat=readstat) multiply_vec
      multiply_cell=.true.
      if (readstat .ne. 0) then
         stop "Check the command -multiply=..., something went wrong!"
      end if
   end if
end do



!
!    PART B: Read in the XDATCAR file !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    First, determine the number of lines in the XDATCAR file
!
call system("wc -l XDATCAR > xdat_length")
open(unit=45,file="xdat_length",status="old")
read(45,*) xdat_lines
close(45)

open(unit=14,file="XDATCAR",status="old",iostat=openstat)
if (openstat .ne. 0) then
   stop "ERROR! The file 'XDATCAR' could not been found!"
end if
read(14,*)
read(14,*) factor   ! the global geometry conversion factor
!
!    Read in the cell dimensions
!
read(14,*) a_vec(1),a_vec(2),a_vec(3)
read(14,*) b_vec(1),b_vec(2),b_vec(3)
read(14,*) c_vec(1),c_vec(2),c_vec(3)
!
!    Read in the elements
!
allocate(el_names_read(10))
el_names_read="XX"
read(14,'(a)') cdum
read(cdum,*,iostat=readstat) el_names_read
nelems=0
do i=1,10
   if (el_names_read(i) .eq. "XX") exit
   nelems=nelems+1
end do
allocate(el_names(nelems),el_nums(nelems))

do i=1,nelems
   el_names(i)=el_names_read(i)
end do
read(14,*) el_nums

!
!    Define the element symbols for all atoms in the system
!
natoms = sum(el_nums)
allocate(at_names(natoms))

counter = 1
do i=1,nelems
   do j=1,el_nums(i)
      at_names(counter) = el_names(i)
      counter = counter +1
   end do
end do

nframes = (xdat_lines - 7)/(natoms+1)

allocate(xyz(3,natoms,nframes))
!
!    Read in the coordinates of the frames and correct for image flags
!
eval_stat = .false.
write(*,*)
write(*,*) "Read in the trajectory from XDATCAR..."
do i=1,nframes
!
!    Every 10% of the read in, give a status update
!
   do j=1,10
      if (real(i)/real(nframes) .gt. real(j)*0.1d0) then
         if (.not. eval_stat(j)) then
            write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
            eval_stat(j) = .true.
         end if
      end if
   end do
   read(14,*)
   do j=1,natoms
      read(14,'(a)') a120
      read(a120,*,iostat=readstat) act_num(:)
!
!   In long trajectories, problems might arise with double-digit negative numbers!
!   (no space between them)   Insert a space again!
!
      if (readstat .ne. 0) then
         a220 = ""
         endl = 1
         do k=1,119
            atest=a120(k:k+1)
            if (index(atest,"-") .ne. 0) then
               write(a220(1:endl+4),'(a,a)') a220(1:endl),"   -"
               endl = endl + 4
            else
               write(a220(1:endl+1),'(a,a)') a220(1:endl),atest
               endl = endl + 1
            end if
         end do
    !  write(*,*) "Problem in XDATCAR file!"
         read(a220,*,iostat=readstat) act_num(:)
      end if

      do k=1,3
         do
            if (act_num(k) > 1d0) then
               act_num(k) = act_num(k) - 1d0
            else
               exit
            end if
         end do
         do
            if (act_num(k) < 0d0) then
!
!    Special case: move atoms near the lower border to negative values
!
!
               if (act_num(k) >= -0.2d0) then
                  exit
               end if
               act_num(k) = act_num(k) + 1d0
            else
               exit
            end if
         end do
!
!     We assume that the bulk is located in the lower half of the simulation
!     cell. If atoms go through the lower x-y surface z-values near 1,
!     move them to values close below zero for better appearance
!

         if (act_num(k) .gt. 0.9d0) then
            act_num(k) = act_num(k)-1.d0
         end if

         if (k .eq. 1) then
            xyz(k,j,i) = act_num(k)*xlen
         end if
         if (k .eq. 2) then
            xyz(k,j,i) = act_num(k)*ylen
         end if
         if (k .eq. 3) then
            xyz(k,j,i) = act_num(k)*zlen
         end if
      end do
   end do
end do
write(*,*) " completed!"
close(14)

write(*,*)
write(*,*) "---------- SETTINGS ---------------------------"
write(*,*) "Number of atoms in the system:",natoms
write(*,*) "Number of frames in the trajectory:",nframes
write(*,*) "-----------------------------------------------"

end program modify_xdatcar        
