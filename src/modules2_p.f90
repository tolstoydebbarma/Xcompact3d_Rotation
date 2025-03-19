  module tools_p
!********************************

! Module for various tools for processing data etc
implicit none
contains

subroutine conv_p(itime,dt,cutOff,ux1,uy1,uz1,uxR1_p,uyR1_p,uzR1_p,ta1)
!#########################

use decomp_2d
use decomp_2d_io
use param_p, only: foldName_p
use var_p, only: dummyString
use MPI

implicit none

real(mytype), intent(in) :: cutOff
real(mytype),dimension(:,:,:) :: ux1,uy1,uz1,uxR1_p,uyR1_p,uzR1_p,ta1
integer :: mpiErr,itime
real(mytype) :: maxErr,totErr,mETemp,tETemp,dt

if(cutOff .ge. 0d0) then   
  ta1 = sqrt((ux1-uxR1_p)**2+(uy1-uyR1_p)**2 + (uzR1_p-uz1)**2) ! uxR1_p no need any longer
  
  mETemp = maxval(ta1); tETemp = sum(ta1)
  
  call MPI_REDUCE(mETemp,maxErr,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,mpiErr)
  call MPI_REDUCE(tETemp,totErr,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,mpiErr)
 
  uxR1_p = ux1; uyR1_p = uy1 ; uzR1_p = uz1 ! Replace references
 

  if(nrank == 0) then
      open(163,file=trim(foldName_p)//'maxErr.txt',status='unknown',position='append')
      write(163,*) itime*dt,' ',maxErr
      close(163)
      open(163,file=trim(foldName_p)//'totErr.txt',status='unknown',position='append')
      write(163,*) itime*dt,' ',totErr
      close(163)
      if(maxErr .le. cutOff) then
         write(*,*) 'Code converged with max diff betwn consecutive viz of: ', maxErr
         call decomp_2d_abort(163,'Code converged')
      endif
  endif
endif

!#########################
end subroutine conv_p


!********************************
  end module tools_p
!********************************
