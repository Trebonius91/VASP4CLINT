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
integer::inc  ! incremented index
integer::readstat
integer::conf_num  ! number of configurations (full structures)
integer::act_conf   ! the current configuration number
integer::natoms_max  ! maximum number of atoms per configuration
integer::natoms_sum  ! sum of all atoms within all configurations
integer::nelems  ! current number of elements
integer::nelems_all  ! total number of different elements
integer::nelems_glob  ! number of elements in current ML_AB
integer::ind_act ! current element index
integer,allocatable::natoms(:)  ! atom numbers of configurations
real(kind=8),allocatable::xyz(:,:,:) ! the geometries
real(kind=8),allocatable::xyz_dir(:,:,:) ! the geometries (direct coords)
real(kind=8),allocatable::cells(:,:,:) ! the unit cells
integer,allocatable::inds(:,:)  ! the element indices (core charges)
integer::ind_list(50)   ! list of different element indices
character(len=2),allocatable::el_list(:)  ! current list of elements
character(len=2)::el_list_glob(50)  ! global list of elements 
real(kind=8),allocatable::energies(:)  ! the energies 
real(kind=8),allocatable::grads(:,:,:)  ! the gradient vectors
real(kind=8)::cell_inv(3,3)  ! inverted coordinate unit cell
real(kind=8)::dist_act   ! scalar distance between two atoms
real(kind=8)::cutoff    ! the distance cutoff during the ML-FF learning
real(kind=8)::dist_vec(3)  ! distance vector between two atoms
integer,allocatable::el_nums(:)  ! current numbers of elements
integer,allocatable::confnum_all(:)  ! numbers of configurations 
integer,allocatable::nat_all(:) ! number of atom in configurations
integer,allocatable::ind_all(:)  ! element indices of all atoms
integer,allocatable::num_around(:,:) ! number of atoms around the atom

real(kind=8),allocatable::gradnorm_all(:)  ! gradient forms for all atoms

character(len=130)::a130


write(*,*) 
write(*,*) "PROGRAM mlff_select: selection of an effective set of"
write(*,*) " local reference configurations (basis sets) for a VASP"
write(*,*) " machine-learning force field from a given ML-AB file."
write(*,*) "The procedure is similar to the VASP command ML_MODE = select,"
write(*,*) " but is faster (since purely based on geometrical comparisons)"
write(*,*) " and has more options to be adjusted (global approach)"
write(*,*)
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
!     Open the ML_AB file and check if it's there
!
open(unit=56,file="ML_AB",iostat=readstat)
if (readstat .ne. 0) then
   write(*,*) "The file ML_AB could not been found!"
   stop
end if
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
              & xyz(1,i,act_conf)*cell_inv(2,:)+xyz(3,i,act_conf)*cell_inv(3,:)
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
num_around=0
write(*,*) "Total number of atoms (possible basis functions):",natoms_sum

inc=0
do i=1,conf_num
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
!     All distances below cutoff
!
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
            dist_vec(:)=dist_vec(1)*cells(1,i,act_conf)+dist_vec(2)*cells(2,i,act_conf)+ &
                       & dist_vec(3)*cells(3,i,act_conf)
!
!     Calculate the distance
!
            dist_act=sqrt(dot_product(dist_vec,dist_vec))
!
!     Determine if atom is within the cutoff. If yes, add it to the list of 
!      surrounding atoms
!
            if (dist_act .lt. cutoff) then
               write(*,*) "dist",dist_act
               do l=1,nelems_all
                  if (inds(k,i) .eq. ind_list(l)) then
                     num_around(l,inc)=num_around(l,inc)+1
                  end if
               end do
            end if
         end if
      end do
      write(*,*) num_around(:,inc)
      stop "Gupgu"

   end do
end do



end program mlff_select

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
