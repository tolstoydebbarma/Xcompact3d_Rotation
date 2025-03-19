  module read_buffer
!********************************

! Module to set common parmeters I will be using
  use decomp_2d, only : mytype, nrank
  use param, only : xlx
  use variables, only : nx
  
  implicit none
  
  integer, parameter :: xPosnBuff = 2	! Length of the buffer
  integer, save :: ixBuffPosn_p ! Grid starting position of buffer
!********************************

contains

!#########################
 subroutine computeBufferPosn
 
 ixBuffPosn_p = nint((xlx-xPosnBuff)/xlx*(nx-1))+1
 
 if (nrank==0) then
	write(*,*) 'ixBuffPosn_p is ', ixBuffPosn_p
 endif
 
 end subroutine computeBufferPosn
!#########################

  end module read_buffer