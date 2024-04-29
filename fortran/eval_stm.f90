!
!    eval_stm: evaluates the results of a STM calculation 
!      (partial charge density, to be evaluated by the Tersoff-Hamann
!      approach). A PARCHG file is read in and depending on the 
!      user keywords, a constant current or constant height STM
!      picture is evaluated. Both a bitmap file and a plottable 
!      data file with a gnuplot file to plot it are given as 
!      results.
!    Part of VASP4CLINT
!     Julien Steffen, 2024 (julien.steffen@fau.de)
!

program eval_stm
implicit none 
integer::readstat
integer::i,j,k,l,m
integer::k_act,l_act,m_act
integer::natoms,nelems
integer::grida,gridb,gridc
integer::nx,ny,nz
integer::stm_gridx,stm_gridy
integer::near_x,near_y,near_z
integer::include_x,include_y,include_z
real(kind=8)::scale
real(kind=8)::a_fac,cutoff,g_width
real(kind=8)::x_len,y_len,z_len
real(kind=8)::stm_height,isos_dens
real(kind=8)::dist_act
real(kind=8)::weight_sum,weight
real(kind=8)::coord_mat(3,3)
real(kind=8)::act_pos(3)
real(kind=8)::vec_act(3)
real(kind=8)::diff_vec(3)
real(kind=8),allocatable::xyz(:,:)
real(kind=8),allocatable::cdens(:,:,:)
real(kind=8),allocatable::c_coord(:,:,:,:)
real(kind=8),allocatable::stm_dat(:,:)
real(kind=8),allocatable::stm_pos(:,:,:)
character(len=120)::cdum,arg
character(len=30)::stm_mode
character(len=2),allocatable::el_names_read(:),el_names(:)
character(len=2),allocatable::at_names(:)
integer,allocatable::el_nums(:)

!
!    Read command line arguments
!
!    If the constant height or the constant current mode 
!      shall be used
!
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:6))  .eq. "-mode=") then
      read(arg(7:),*,iostat=readstat) stm_mode
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -mode=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
if (trim(stm_mode) .eq. "height") then
   write(*,*) "The constant height mode will be evaluated."
else if (trim(stm_mode) .eq. "current") then
   write(*,*) "The constant current mode will be evaluated."
else 
   write(*,*) "Please give either the constant height mode (-mode=height) "
   write(*,*) " or the constant current mode (-mode=current)!"
   stop
end if        
!
!    If constant height mode: read in the height (z/c coordinate)
!
stm_height=-1.d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:5))  .eq. "-pos=") then
      read(arg(6:),*,iostat=readstat) stm_height
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -pos=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
if (trim(stm_mode) .eq. "height") then
   if (stm_height .lt. 0.d0) then
      write(*,*) "Please give the height (z-coordinate) at which the STM"
      write(*,*) " tip shall scan the surface!"
      stop
   else
      write(*,'(a,f13.6,a)') " The STM tip will scan the surface at z=", &
                      & stm_height," Angstroms" 
   end if        
end if        
!
!    If constant current mode: read in the isos density 
!
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:6))  .eq. "-dens=") then
      read(arg(7:),*,iostat=readstat) isos_dens
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -dens=..., something went wrong!"
         write(*,*)
      end if
   end if
end do

!
!     The number of grid points for STM printout along the x axis 
!
!     Default
!
stm_gridx=100
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:8))  .eq. "-grid_x=") then
      read(arg(9:),*,iostat=readstat) stm_gridx
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -grid_x=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
!
!     The number of grid points for STM printout along the y axis
!
!     Default
!
stm_gridy=100
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:8))  .eq. "-grid_y=") then
      read(arg(9:),*,iostat=readstat) stm_gridy
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -grid_y=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
!
!     The width at half maximum (FWHM) of the Gaussian used for 
!       smoothing the picture
!
!     Default
!
g_width=0.3d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:13))  .eq. "-gauss_width=") then
      read(arg(14:),*,iostat=readstat) g_width
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -gauss_width=..., something went wrong!"
         write(*,*)
      end if
   end if
end do

!
!    Determine the properties and the function of the smoothing Gaussian
!
a_fac=dlog(2.d0)/((g_width/2d0)**2)
cutoff=sqrt(dlog(100d0)/a_fac)

write(*,*) cutoff
!
!    Open file PARCHG and read in relevant data
!
open(unit=15,file="PARCHG",status="old",iostat=readstat)
if (readstat .ne. 0) then
   write(*,*) "The file PARCHG is not there!"
   stop
end if        
!
!    Read the unit cell definition
!
read(15,*)
read(15,*) scale
read(15,*) coord_mat(1,:)
read(15,*) coord_mat(2,:)
read(15,*) coord_mat(3,:)
!
!    The lengths of coordinate axes
!
x_len=sqrt(dot_product(coord_mat(1,:),coord_mat(1,:)))
y_len=sqrt(dot_product(coord_mat(2,:),coord_mat(2,:)))
z_len=sqrt(dot_product(coord_mat(3,:),coord_mat(3,:)))
!
!    Read in the elements
!
allocate(el_names_read(10))
el_names_read="XX"
read(15,'(a)') cdum
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
read(15,*) el_nums

!
!    Define the element symbols for all atoms in the system
!
natoms = sum(el_nums)
allocate(at_names(natoms))
k=0
do i=1,nelems
   do j=1,el_nums(i)
      k=k+1
      at_names(k)=el_names(i)
   end do
end do
!
!    Read in the atomic coordinates
!
allocate(xyz(3,natoms))
read(15,*)
do i=1,natoms
   read(15,*) xyz(:,i)
end do
!
!    Read in the number of grid points
!
read(15,*)
read(15,*) grida,gridb,gridc
allocate(cdens(grida,gridb,gridc))
!
!    Read in the partial charge density 
!    Use compactified read command for array
!
read(15,*) (((cdens(nx,ny,nz),nx=1,grida),ny=1,gridb),nz=1,gridc)
close(15)

!
!    Calculate x,y,z cartesian coordinates for all grid points, for the 
!    subsequent determination of distances to the current sampling point
!
allocate(c_coord(3,grida,gridb,gridc))
do i=1,gridc
   do j=1,gridb
      do k=1,grida
         vec_act(1)=real((k-1))/real(grida)
         vec_act(2)=real((j-1))/real(gridb)
         vec_act(3)=real((i-1))/real(gridc)
         c_coord(:,k,j,i)=matmul(vec_act,coord_mat)
      end do
   end do
end do    
!
!    Allocate data array for 2D STM picture
!
allocate(stm_dat(stm_gridy,stm_gridx))
allocate(stm_pos(2,stm_gridy,stm_gridx))
stm_dat=0.d0
!
!    Determine number of grid points to be included into the Gaussian
!    smoothin region
!
include_x=nint(real(cutoff)/real(x_len)*real(grida))
include_y=nint(real(cutoff)/real(y_len)*real(gridb))
include_z=nint(real(cutoff)/real(z_len)*real(gridc))
!
!    A: CONSTANT HEIGHT MODE
!
if (stm_mode .eq. "height") then
!
!    Desired height/z-coordinate in internal coordinates
!
  
   write(*,*) "pos",act_pos(3)
   write(*,*) include_x,include_y,include_z
!
!    Loop over all pricture gridpoints in the x-y (or a-b) plane and calculate 
!    the Gaussian-weighted interpolation of nearby grid points
!
   do i=1,stm_gridx
      do j=1,stm_gridy
         act_pos(1)=real((i-1))/real(stm_gridx)
         act_pos(2)=real((j-1))/real(stm_gridy)
         act_pos(3)=stm_height/z_len
!
!    Determine nearest point in density grid
!
         near_x=nint(act_pos(1)*grida)
         near_y=nint(act_pos(2)*gridb)
         near_z=nint(act_pos(3)*gridc)

         act_pos(:)=matmul(act_pos,coord_mat)
!
!    Now loop over all other grid points within the proposed cutoff of the 
!    Gaussian smearing function
!
         weight_sum=0
         do k=near_x-include_x,near_x+include_x
            if (k .lt. 1) then
               k_act=grida+k
            else if (k .gt. grida) then
               k_act=k-grida
            else     
               k_act=k
            end if     
            do l=near_y-include_y,near_y+include_y
               if (l .lt. 1) then
                  l_act=gridb+l
               else if (l .gt. gridb) then
                  l_act=l-gridb     
               else 
                  l_act=l
               end if     
               do m=near_z-include_z,near_z+include_z
                  if (m .lt. 1) then
                     m_act=gridc+m
                  else if (m .gt. gridc) then
                     m_act=m-gridc
                  else
                     m_act=m
                  end if   
                  diff_vec(:)=act_pos(:)-c_coord(:,k_act,l_act,m_act)                  
                  dist_act=sqrt(dot_product(diff_vec,diff_vec))
                  weight=exp(-g_width*dist_act*dist_act)
                  weight_sum=weight_sum+weight
                  stm_dat(j,i)=stm_dat(j,i)+weight*cdens(k_act,l_act,m_act)
               end do
            end do
         end do      
         stm_dat(j,i)=stm_dat(j,i)/weight_sum
         stm_pos(:,j,i)=act_pos(1:2)

         !    write(*,*) near_x,near_y,near_z
!
!    Now sum over all grid points within the chosen cutoff for consideration
!    and add the components together
!
      end do
   end do    
end if        

!
!    Write plot file with STM data
!
open(unit=37,file="stm_plot.dat",status="replace")
write(37,*) "# This data has been generated by the program eval_stm."
write(37,*) "# x-coord     y-coord     intenstiy  "
do i=1,stm_gridx
   do j=1,stm_gridy
      write(37,*) stm_pos(1,j,i),stm_pos(2,j,i),stm_dat(j,i)
   end do
   write(37,*)
end do   
close(37)



end program eval_stm        
