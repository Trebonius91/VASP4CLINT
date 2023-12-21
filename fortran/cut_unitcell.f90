!
!    cut_unitcell: cut an arbitrary unit cell from a 
!       given large surface cell
!    Part of VASP4CLINT
!     Julien Steffen, 2023 (julien.steffen@fau.de)
!

program cut_unitcell
implicit none
integer::i,j,k 
real(kind=8)::coord_scale
real(kind=8)::a_vec(3),b_vec(3),c_vec(3)
character(len=2)::el_name
integer::natoms,natoms_chosen
integer::ind1,ind2,ind3,ind4
character(len=1),allocatable::active(:,:)
real(kind=8),allocatable::xyz_ref(:,:)
real(kind=8)::pos1(2),pos2(2),pos3(2),pos4(2)
real(kind=8)::vec1(2),vec2(2),origin(2)
real(kind=8)::vec_final(2,2),pos_act(2)
integer::vec_inds(6,2)
real(kind=8)::angle,arg
real(kind=8)::pi
real(kind=8)::par_x,par_y,x_act,y_act
real(kind=8)::x_check,y_check,dist
real(kind=8)::xyz_chosen(1000,3)
character(len=1)::active_chosen(1000,3)

pi=3.14159265359d0

write(*,*) "This program is able to cut an (almost) arbitrary unit cell"
write(*,*) " from a large given surface."
write(*,*) "Two files must be given:"
write(*,*) " a) POSCAR_surf, containing the large reference unit cell"
write(*,*) " b) atom_inds.txt, containing three atom numbers, each in a line"
write(*,*) "    These indices should be chosen in advance, e.g., with VESTA"
write(*,*) "The unit cell vectors will be set up in the following way:"
write(*,*) " vector1=atom2-atom1"
write(*,*) " vector2=atom3-atom1"
!
!     First, read in the given POSCAR of the large surface
!

open(unit=16,file="POSCAR_surf",status="old")
read(16,*)
read(16,*) coord_scale
read(16,*) a_vec(:)
read(16,*) b_vec(:)
read(16,*) c_vec(:)
read(16,*) el_name
read(16,*) natoms
read(16,*) 
read(16,*)
allocate(xyz_ref(natoms,3))
allocate(active(natoms,3))
do i=1,natoms
   read(16,*) xyz_ref(i,:),active(i,:)
end do
close(16)


!
!     Second, read indices of chosen atoms from file atom_inds.txt
!

open(unit=17,file="atom_inds.txt",status="old")
read(17,*) ind1
read(17,*) ind2
read(17,*) ind3
close(17)

!
!     Now locate all atoms given in atom_inds.txt     
!

pos1=xyz_ref(ind1,1:2)
pos2=xyz_ref(ind2,1:2)
pos3=xyz_ref(ind3,1:2)


!
!     Calculate the unit cell vectors 
!     v1: atom2-atom1
!     v2: atom3-atom1
!     Origin: atom1
!   

vec1=pos2-pos1
vec2=pos3-pos1

origin=pos1

write(*,*) "Unit cell vector 1: ", vec1
write(*,*) "Unit cell vector 2: ", vec2
write(*,*) "origin",origin
!
!     Determine all atoms that are in the defined unit cell
!     Solve linear system of equations for each atom in the system
!     Subtract position of atom 1 in the unit cell to shift them to the 
!     right position
!
natoms_chosen=0
outer: do i=1,natoms
   x_act=xyz_ref(i,1)
   y_act=xyz_ref(i,2)
   par_x=(x_act*vec2(2)-vec2(1)*y_act)/(vec1(1)*vec2(2)-vec2(1)*vec1(2))
   par_y=(y_act*vec1(1)-vec1(2)*x_act)/(vec1(1)*vec2(2)-vec2(1)*vec1(2))
   if ((par_x .ge. 0.d0) .and. (par_x .lt. 1.d0) .and. (par_y .ge. 0.d0) &
              &  .and. (par_y .lt. 1.d0)) then
!
!     Check if another atom within the cell would be (nearly) reproced by
!     shifting the new atom with the unit cell vectors
!
!      x_check=x_act+vec1(1)
!      y_check=y_act+vec1(2)
!      inner1: do j=1,natoms_chosen
!         dist=sqrt((x_check-xyz_chosen(j,1))**+(y_check-xyz_chosen(j,2))**2)
!         if (dist .lt. 1E-6) cycle outer
!      end do inner1
!      x_check=x_act-vec1(1)
!      y_check=y_act-vec1(2)
!      inner2: do j=1,natoms_chosen
!         dist=sqrt((x_check-xyz_chosen(j,1))**+(y_check-xyz_chosen(j,2))**2)
!         if (dist .lt. 1E-6) cycle outer
!      end do inner2
!      x_check=x_act+vec2(1)
!      y_check=y_act+vec2(2)
!      inner3: do j=1,natoms_chosen
!         dist=sqrt((x_check-xyz_chosen(j,1))**+(y_check-xyz_chosen(j,2))**2)
!         if (dist .lt. 1E-6) cycle outer
!      end do inner3
!      x_check=x_act-vec2(1)
!      y_check=y_act-vec2(2)
!      inner4: do j=1,natoms_chosen
!         dist=sqrt((x_check-xyz_chosen(j,1))**+(y_check-xyz_chosen(j,2))**2)
!         if (dist .lt. 1E-6) cycle outer
!      end do inner4
!
      natoms_chosen=natoms_chosen+1
      xyz_chosen(natoms_chosen,:)=xyz_ref(i,:)
      active_chosen(natoms_chosen,:)=active(i,:)
   end if
end do outer
!
!     Write resulting POSCAR with new unit cell
!
open(unit=27,file="POSCAR",status="replace")
write(27,*) "New unit cell, build by cut_unitcell.f90"
write(27,*) coord_scale
write(27,*) vec1(1),vec1(2),0.d0
write(27,*) vec2(1),vec2(2),0.d0
write(27,*) 0.d0,0.d0,c_vec(3)
write(27,*) el_name
write(27,*) natoms_chosen
write(27,*) "Selective "
write(27,*) "Cartesian "
do i=1,natoms_chosen
   write(27,*) xyz_chosen(i,:),active_chosen(i,1), "  ",active_chosen(i,2),&
                & "   ",active_chosen(i,3)
end do


close(27)


end program cut_unitcell
