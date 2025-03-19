!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module jet_nonozzle

  use decomp_2d
  use variables
  use param
  !!$use complex_geometry, only: tol

  implicit none

  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_jet_nonozzle, boundary_conditions_jet_nonozzle

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  subroutine: init_pipe
  !!      AUTHOR: Rodrigo Vicente Cruz
  !! DESCRIPTION: Initial laminar conditions in pipe
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !********************************************************************
  !
  subroutine init_jet_nonozzle (ux1,uy1,uz1,phi1)
  !
  !********************************************************************
    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI
    use dbg_schemes, only: exp_prec

	use param_swirling_nonozzle
	
    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
	
	! re5strt problem
	real(mytype) :: ym,zm,r,um,yc,zc
	
	real(mytype) :: ux_p,uy_p,uz_p,uth_p    						!@@6A velocity components
	real(mytype) :: sub_fr_r
	
	! Make changes in ecoule_nonozzle too
	
    integer :: k,j,i,ii,is,code
	
	if (nrank==0) then 
		write(*,*) 'alp_nonz = ', alp_nonz
		write(*,*) 'del_nonz = ', del_nonz
		write(*,*) 'del_cf_nonz = ', del_cf_nonz
		write(*,*) 'S_nonz = ', S_nonz
		write(*,*) 'rao_nonz = ', rao_nonz
	endif
	
	yc = yly / two
    zc = zlz / two

    ux1=zero; uy1=zero; uz1=zero
	sub_fr_r=1d0-2d0*del_nonz
	
	if (iscalar==1) then

       phi1(:,:,:,:) = zero !change as much as you want

    endif
	
    if (iin.ne.0) then
       call system_clock(count=code)
       if (iin.eq.2) code=0
       call random_seed(size = ii)
       call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /))

       call random_number(ux1)
       call random_number(uy1)
       call random_number(uz1)

       do k=1,xsize(3)
          do j=1,xsize(2)
             do i=1,xsize(1)
                ux1(i,j,k)=init_noise*ux1(i,j,k)
                uy1(i,j,k)=init_noise*uy1(i,j,k)
                uz1(i,j,k)=init_noise*uz1(i,j,k)
             enddo
          enddo
       enddo
		
	   ! re5strt problem
	   do k=1,xsize(3)
	      zm=dz*real(xstart(3)-1+k-1,mytype)-zc 
          do j=1,xsize(2)
	         ym=real(j+xstart(2)-1-1,mytype)*dy-yc
		     r=sqrt(ym*ym+zm*zm)
		     !****************** Maxworthy *************************
		     ux_p = 1d0 - 0.5d0*(1d0+erf((r-sub_fr_r)/del_nonz))
		     uth_p = S_nonz*(r/two)*(1d0 - erf((r-sub_fr_r)/del_nonz))
		     
			 if (r.gt.rao_nonz) then
				ux_p = ux_p+erf((r-rao_nonz)/del_cf_nonz)*alp_nonz
			 endif
			 
			 if (r .ne. 0) then
			    uz_p =  uth_p*ym/r 		!@@P Prev error correct, right-hand coord now!
			    uy_p = -uth_p*zm/r 		!@@P y-z to r-th (uy=urcos-uthsin,uz=ursin+uthcos)
		     else
			    uz_p = 0d0
			    uy_p = 0d0
		     endif
		     bxx1(j,k) = ux_p
		     bxy1(j,k) = uy_p
		     bxz1(j,k) = uz_p
		  enddo
	   enddo
	   ! Make changes in ecoule_nonozzle too

       !modulation of the random noise + inlet velocity profile
       do k=1,xsize(3)
		  zm=dz*real(xstart(3)-1+k-1,mytype)-zlz/2d0
          do j=1,xsize(2)
             if (istret.eq.0) ym=(j+xstart(2)-1-1)*dy-yly/2d0
             if (istret.ne.0) ym=yp(j+xstart(2)-1)-yly/2d0
             r=sqrt(ym*ym+zm*zm)
			 um=exp_prec(-zptwo*r*r)
			 do i=1,xsize(1)
                ux1(i,j,k)=um*ux1(i,j,k)+bxx1(j,k) ! + u1
                uy1(i,j,k)=um*uy1(i,j,k)+bxy1(j,k)
                uz1(i,j,k)=um*uz1(i,j,k)+bxz1(j,k)
             enddo
          enddo
       enddo
    endif

#ifdef DEBG
    if (nrank .eq. 0) write(*,*) '# init end ok'
#endif

    return

  end subroutine init_jet_nonozzle
  !********************************************************************
  !
  subroutine boundary_conditions_jet_nonozzle (ux,uy,uz,phi)
  !
  !********************************************************************

    use param
    use variables
    use decomp_2d

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    !
	
	call inflow (ux,uy,uz,phi)
    call outflow (ux,uy,uz,phi)
	
    return
    !
  end subroutine boundary_conditions_jet_nonozzle
  
  subroutine ecoule_nonozzle
	USE param
	use variables
    USE decomp_2d
	
	use param_swirling_nonozzle
	
	implicit none
	
	integer  :: j,k
	
	! re5strt problem
	real(mytype) :: ym,zm,r,yc,zc
	
	real(mytype) :: ux_p,uy_p,uz_p,uth_p    						!@@6A velocity components
	real(mytype) :: sub_fr_r
	
	! Make changes in ecoule_nonozzle too
	
	yc = yly / two
    zc = zlz / two
	
	sub_fr_r=1d0-2d0*del_nonz
	
		! re5strt problem
	   do k=1,xsize(3)
	      zm=dz*real(xstart(3)-1+k-1,mytype)-zc 
          do j=1,xsize(2)
	         ym=real(j+xstart(2)-1-1,mytype)*dy-yc
		     r=sqrt(ym*ym+zm*zm)
		     !****************** Maxworthy *************************
		     ux_p = 1d0 - 0.5d0*(1d0+erf((r-sub_fr_r)/del_nonz))
		     uth_p = S_nonz*(r/two)*(1d0 - erf((r-sub_fr_r)/del_nonz))
		     
			 if (r.gt.rao_nonz) then
				ux_p = ux_p+erf((r-rao_nonz)/del_cf_nonz)*alp_nonz
			 endif
			 
			 if (r .ne. 0) then
			    uz_p =  uth_p*ym/r 		!@@P Prev error correct, right-hand coord now!
			    uy_p = -uth_p*zm/r 		!@@P y-z to r-th (uy=urcos-uthsin,uz=ursin+uthcos)
		     else
			    uz_p = 0d0
			    uy_p = 0d0
		     endif
		     bxx1(j,k) = ux_p
		     bxy1(j,k) = uy_p
		     bxz1(j,k) = uz_p
		  enddo
	   enddo
		! Make changes in ecoule_nonozzle too
	
  end subroutine ecoule_nonozzle
  !********************************************************************
  subroutine inflow (ux,uy,uz,phi)

    USE param
    USE variables
    USE decomp_2d
    USE ibm_param

    implicit none

    integer  :: j,k,is
	real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
	
    call random_number(bxo)
    call random_number(byo)
    call random_number(bzo)
	
	if (irestart==1) then
		call ecoule_nonozzle
	endif	
	
	do k=1,xsize(3)
       do j=1,xsize(2)
		  bxx1(j,k) = bxx1(j,k)+bxo(j,k)*inflow_noise
		  bxy1(j,k) = bxy1(j,k)+bxo(j,k)*inflow_noise
		  bxz1(j,k) = bxz1(j,k)+bxo(j,k)*inflow_noise	
       enddo
    enddo
	
    if (iscalar.eq.1) then
       do is=1, numscalar
          do k=1,xsize(3)
             do j=1,xsize(2)
                phi(1,j,k,is)=cp(is)
             enddo
          enddo
       enddo
    endif

    return
  end subroutine inflow
  !********************************************************************
  
  subroutine outflow (ux,uy,uz,phi)

    USE param
    USE variables
    USE decomp_2d
    USE MPI
    USE ibm_param

    implicit none

    integer :: j,k,code
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx,uxmin,uxmax,uxmin1,uxmax1

    udx=one/dx; udy=one/dy; udz=one/dz; uddx=half/dx; uddy=half/dy; uddz=half/dz

    uxmax=-1609._mytype
    uxmin=1609._mytype
    do k=1,xsize(3)
       do j=1,xsize(2)
          if (ux(nx-1,j,k).gt.uxmax) uxmax=ux(nx-1,j,k)
          if (ux(nx-1,j,k).lt.uxmin) uxmin=ux(nx-1,j,k)
       enddo
    enddo

    call MPI_ALLREDUCE(uxmax,uxmax1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(uxmin,uxmin1,1,real_type,MPI_MIN,MPI_COMM_WORLD,code)

    if (u1 == zero) then
       cx=(half*(uxmax1+uxmin1))*gdt(itr)*udx
    elseif (u1 == one) then
       cx=uxmax1*gdt(itr)*udx
    elseif (u1 == two) then
       cx=u2*gdt(itr)*udx    !works better
    else
       cx=(half*(u1+u2))*gdt(itr)*udx
    endif

    do k=1,xsize(3)
       do j=1,xsize(2)
          bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
          bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
          bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
       enddo
    enddo

    if (iscalar==1) then
       if (u2==zero) then
          cx=(half*(uxmax1+uxmin1))*gdt(itr)*udx
       elseif (u2==one) then
          cx=uxmax1*gdt(itr)*udx
       elseif (u2==two) then
          cx=u2*gdt(itr)*udx    !works better
       else
          stop
       endif

       do k=1,xsize(3)
          do j=1,xsize(2)
             phi(nx,j,k,:)=phi(nx,j,k,:)-cx*(phi(nx,j,k,:)-phi(nx-1,j,k,:))
          enddo
       enddo
    endif

    if (nrank==0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime == ilast)) &
       write(*,*) "Outflow velocity ux nx=n min max=",real(uxmin1,4),real(uxmax1,4)

    return
  end subroutine outflow
  !********************************************************************

end module jet_nonozzle
