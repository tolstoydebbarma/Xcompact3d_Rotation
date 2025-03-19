module param_swirling
	use decomp_2d, only : mytype
	implicit none
	
	real(mytype), save :: alp_p,del_p,del_cf,S_p, len_nz
	
end module param_swirling

module param_swirling_nonozzle
	use decomp_2d, only : mytype
	implicit none
	
	real(mytype), save :: alp_nonz,del_nonz,del_cf_nonz,S_nonz,rao_nonz
	
end module param_swirling_nonozzle

  module param_p
!********************************

! Module to set common parmeters I will be using
  use decomp_2d, only : mytype
  implicit none
  
  character, parameter :: foldName_p*4='out/'    ! Output Folder name
  
  real(mytype), parameter :: convT=5d-1
  
  real(mytype), parameter :: blowUpVel_p=1d1                       ! Velocity above which code is assumed to have blown
  
  real(mytype), parameter :: convVelDiff_p=1d-4                     ! Velocity difference between two successive visualizations, freq of viz
  
  integer, save :: convModulo_p  !@@ time iteration difference between two adjacent convergence velocities
  
  end module param_p
  

  module var_p 
!********************************

! Module to initialize common variables I will be using

 use decomp_2d, only : mytype
 implicit none

 character,save :: dummyString*200                         ! String used for multiple purposes

! Probes
  integer :: totPb_p						  	  !@@7A total no of probes	
  integer, allocatable, save, dimension(:) :: pbXPos_p,pbYPos_p,pbZPos_p  !@@7A Probe positions
  real(mytype),save,allocatable,dimension(:,:,:) :: pbU_p ! Probe save till isave 

  end module var_p
  
!#########################
  subroutine readProbes_p
!#########################

! Subroutine to read and initialize probe location

use decomp_2d, only: mytype
use var_p
use variables, only: nxm,nym,nzm
use param, only: xlx,yly,zlz,irestart,icheckpoint

implicit none

real(mytype) :: xPb,yPb,zPb	! Probe positions
integer :: i

!Probe file: probes.prm read
open(891,file = 'probes_points.prm',status='old')
read(891,*) dummyString
read(891,*) dummyString
read(891,*) dummyString

! No of probes
read(891,*) totPb_p
read(891,*) dummyString

! Allocate that many probes
allocate(pbXPos_p(totPb_p))
allocate(pbYPos_p(totPb_p))
allocate(pbZPos_p(totPb_p))
allocate(pbU_p(totPb_p,3,icheckpoint))

! Create files if starting afresh
if(irestart.eq.0) then
do i=1,totPb_p
   write(dummyString,'(a,i3.3)') 'probes/Pb',i
   open(892,file=trim(dummyString),status='replace')
   close(892)
enddo
endif

! Read positions & calculate indices
do i=1,totPb_p
   read(891,*) xPb,yPb,zPb
! Probe Global indices
   pbxPos_p(i) = nint(xPb/xlx*nxm)+1
   pbyPos_p(i) = nint(yPb/yly*nym)+1+nym/2
   pbzPos_p(i) = nint(zPb/zlz*nzm)+1+nzm/2
enddo

close(891)

!#########################
  end subroutine readProbes_p
!#########################

!#########################
subroutine probes_p(ux,uy,uz,t)
!#########################

! Subroutine to write velocity at given probe position
use var_p
use decomp_2d, only: mytype, xsize,xstart,xend
use param, only: icheckpoint
implicit none

real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
integer, intent(in) :: t
integer :: i,j,k,pbNo

! Expensive way of doing

do pbNo = 1,totPb_p

  if(pbyPos_p(pbNo) .ge. xstart(2) .and. pbyPos_p(pbNo) .le. xend(2) ) then
  if(pbzPos_p(pbNo) .ge. xstart(3) .and. pbzPos_p(pbNo) .le. xend(3) ) then

    i = pbxPos_p(pbNo)
    j = pbyPos_p(pbNo)-xstart(2)+1
    k = pbzPos_p(pbNo)-xstart(3)+1
    pbU_p(pbNo,1,mod(t-1,icheckpoint)+1) = ux(i,j,k)
    pbU_p(pbNo,2,mod(t-1,icheckpoint)+1) = uy(i,j,k)
    pbU_p(pbNo,3,mod(t-1,icheckpoint)+1) = uz(i,j,k)

    if(mod(t,icheckpoint) == 0) then
       write(dummyString,'(a,i3.3)') 'probes/Pb',pbNo
       open(892,file=trim(dummyString),status='unknown',action='write', &
       form='formatted',position="append")
       do i = 1,icheckpoint
           write(892,*)  t-icheckpoint+i, pbU_p(pbNo,1,i), pbU_p(pbNo,2,i), pbU_p(pbNo,3,i)
       enddo
       close(892)
    endif

  endif
  endif

enddo

!#########################
end subroutine probes_p
!#########################  