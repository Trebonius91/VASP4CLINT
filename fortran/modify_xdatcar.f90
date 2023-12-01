!
!    modify_xdatcar: Perform different tasks on a XDATCAR file,
!      such as translating or multiplying the unit cell, picking
!      certain structures or write an xyz trajectory file
!

program modify_xdatcar
implicit none 
integer::i,j,k
integer::readstat,openstat
real(kind=8)::a_vec(3),b_vec(3),c_vec(3)


!
!    PART A: Read in the XDATCAR file !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!


!
!    First, determine the number of lines in the XDATCAR file
!
if (use_reaxff) then
   call system("wc -l dump.xyz > dump_length")
   open(unit=45,file="dump_length",status="old")
   read(45,*) dump_lines
   close(45)
else
   call system("wc -l XDATCAR > xdat_length")
   open(unit=45,file="xdat_length",status="old")
   read(45,*) xdat_lines
   close(45)
end if

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
   if (el_names_read(i) .eq. cls_element) then
      cls_elem_ind=nelems
   end if
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
!
!    If the XDATCAR file shall not be read in, skip the rest
!
if (skip_xdat) then
   close(14)
   goto 33
end if

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


end program modify_xdatcar        
