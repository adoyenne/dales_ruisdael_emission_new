!> \file modemission.f90
!!  (Anthropogenic) emissions

!>
!!  \author Marco de Bruine and Arseni Doyennel, VU
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2020 Delft University of Technology, Wageningen
!  University, Utrecht University, KNMI
!

module modemission

use modemisdata

implicit none

contains

  subroutine initemission

    use modglobal,    only : i2, j2,kmax, nsv, ifnamopt, fname_options, checknamelisterror
    use modmpi,       only : myid, comm3d, mpi_logical, mpi_integer, mpi_character
    use moddatetime,  only : datex, prevday, nextday

    implicit none
  
    ! Auxiliary variables
    integer :: ierr
    
     ! Declare variables to hold next and previous hours
  	integer :: next_hour, prev_hour

    ! --- Read & broadcast namelist EMISSION -----------------------------------
    namelist/NAMEMISSION/ l_emission, l_points, kemis, svskip, emisnames, svco2ags, svco2sum, nullemiss

	nullemiss(1:10) = -1

    if (myid == 0) then

      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMEMISSION,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMEMISSION')
      write(6, NAMEMISSION)
      close(ifnamopt)

    endif

    call mpi_bcast(l_emission,    1,   mpi_logical,   0, comm3d, ierr)
    call mpi_bcast(l_points,      1,   mpi_logical,   0, comm3d, ierr)
    call mpi_bcast(kemis,         1,   mpi_integer,   0, comm3d, ierr)
    call mpi_bcast(svskip,        1,   mpi_integer,   0, comm3d, ierr)
	call mpi_bcast(emisnames(1:100), 100, mpi_character, 0, comm3d, ierr)
    call mpi_bcast(svco2ags,        1,   mpi_integer,   0, comm3d, ierr)
	call mpi_bcast(svco2sum,        1,   mpi_integer,   0, comm3d, ierr)
	call mpi_bcast(nullemiss(1:10), 10, mpi_integer, 0, comm3d, ierr)
	

    ! --- Local pre-calculations and settings
    if (.not. (l_emission)) return

    ! -- Interaction with AGs   ----------------------------------------------------
     
    allocate(co2fields(nsv))
    co2fields = 0
    co2fields(svskip+1:nsv) = index(emisnames(1:nsv-svskip), "co2") 

    !svco2ags = findloc(emisnames(1:nsv-svskip), value = "co2ags", dim = 1)
    !svco2ags = svco2ags + svskip
 
    if (myid == 0) then
      write(6,*) 'modemission: co2fields (scalar fields with CO2 0=no, 1=yes)'
      write(6,*) co2fields
      write(6,*) 'modemission: svco2ags (scalar field number for AGS emissions)'
      write(6,*) svco2ags
    endif

    if (kemis == -1) kemis = kmax

    ! --- Read emission files for first time step ----------------------------------

    ! Two hourly emission fields are loaded at all times: 
    ! (1) before model time,   t_field < t_model, "in-the-past"
    ! (2) ahead of model time, t_field > t_model, "in-the-future"
    
    allocate(emisfield(i2, j2, kemis, svskip+1:nsv, 2))

    if (datex(5) >= 30) then
      call reademission(    datex(1),   datex(2),   datex(3),   datex(4), emisfield(:,:,:,:,1))
      
      if (datex(4) == 23) then  
        call reademission(nextday(1), nextday(2), nextday(3),          0, emisfield(:,:,:,:,2))
      else
        call reademission(  datex(1),   datex(2),   datex(3), datex(4)+1, emisfield(:,:,:,:,2))
      endif

    else
      call reademission(    datex(1),   datex(2),   datex(3),   datex(4), emisfield(:,:,:,:,2)) 
     
      if (datex(4) == 0) then
        call reademission(prevday(1), prevday(2), prevday(3),         23, emisfield(:,:,:,:,1))
      else
        call reademission(  datex(1),   datex(2),   datex(3), datex(4)-1, emisfield(:,:,:,:,1))
      endif

    endif
	

if (l_points) then
  npoints = 0

  
  ! Check if the current time step is at the half-hour mark
  if (datex(5) >= 30) then
  	! Check if you need to read data for the current time step
  	call inquirepoints(datex(1), datex(2), datex(3), datex(4))
    if (datex(4) == 23) then
      call inquirepoints(nextday(1), nextday(2), nextday(3), 0)
    else
      call inquirepoints(datex(1), datex(2), datex(3), datex(4) + 1)
    endif
  else
  	call inquirepoints(datex(1),   datex(2),   datex(3),   datex(4))
    if (datex(4) == 0) then
      call inquirepoints(prevday(1), prevday(2), prevday(3), 23)
    else
      call inquirepoints(datex(1), datex(2), datex(3), datex(4) - 1)
    endif
  endif

  ! Check if there are points to read
  if (npoints > 0) then
  
  	allocate(point_source_data(npoints, 7,2))
    
    if (datex(5) >= 30) then

    	call readpoints(datex(1), datex(2), datex(3), datex(4), point_source_data(:,:,1))

      	if (datex(4) == 23) then
      	
        	call readpoints(nextday(1), nextday(2), nextday(3), 0, point_source_data(:,:,2))
      	else
        	
        	call readpoints(datex(1), datex(2), datex(3), datex(4) + 1, point_source_data(:,:,2))
      	endif
    else
      	
        call readpoints(datex(1),   datex(2),   datex(3),   datex(4), point_source_data(:,:,2))
        
      	if (datex(4) == 0) then
        	
        	call readpoints(prevday(1), prevday(2), prevday(3), 23, point_source_data(:,:,1))
      	
      	else
        	
        	call readpoints(datex(1), datex(2), datex(3), datex(4) - 1, point_source_data(:,:,1))
      endif
    endif
  endif
endif

  end subroutine initemission

  subroutine reademission(iyear, imonth, iday, ihour, emisfield)

  ! ----------------------------------------------------------------------
  ! Reading of emission files
  ! Multiple/all tracers
  ! ----------------------------------------------------------------------

    use netcdf
    use modmpi,      only : myid, myidx, myidy
    use modglobal,   only : i1, j1, i2, j2, imax, jmax, nsv
    use moddatetime, only : datex

    implicit none

    integer, intent(in)  :: iyear, imonth, iday, ihour     
    real, intent(out)    :: emisfield(i2, j2, kemis, 1+svskip:nsv)

    integer, parameter   :: ndim = 3 
    integer              :: start(ndim), count(ndim)
    integer              :: ncid, varid
    integer              :: isv
 
    character(len=12)    :: sdatetime

    ! Create string from given date
    write(sdatetime, "(I0.4,2I0.2,2I0.2)") iyear, imonth, iday, ihour, 0 

    write(6,"(A18, A12)") "Reading emission: ", sdatetime

    do isv = 1+svskip, nsv

      call check( nf90_open( trim(emisnames(isv-svskip))//'_emis_'//sdatetime//'_3d.nc', NF90_NOWRITE, ncid))

      call check( nf90_inq_varid( ncid, emisnames(isv-svskip), varid) )
      call check( nf90_get_var  ( ncid, varid, emisfield(2:i1,2:j1,1:kemis,isv), & 
                                  start = (/1 + myidx * imax, 1 + myidy * jmax, 1, 1/), &
                                  count = (/imax, jmax, kemis, 1/) ) )
      call check( nf90_close( ncid ) )

    end do

  contains

    subroutine check(status)
      integer, intent(in) :: status
    
      if(status /= nf90_noerr) then 
        print *, trim(nf90_strerror(status))
        stop 'NetCDF error in modemission. See outputfile for more information.'
      end if
    end subroutine check

  end subroutine reademission

  subroutine emission
  ! ----------------------------------------------------------------------
  ! Read appropriate emission fields, interpolate and transfer to svp
  !
  ! NOTES
  ! 1. Emission files (currently) in kg per gridbox per hour!
  !    What results from this routine now is ug/g (same as CAMS input), i.e. we scale for time,
  !    gridbx size and air density AND apply a factor of 1e6.
  ! 
  ! 
  !    Note that svp is tracer tendency in ug g-1 s-1
  !
  ! TODO
  ! 1. MDB Align properly with "non-emitted" tracers, i.e. cloud scalars from e.g.
  ! microphysics/chemistry
  ! ----------------------------------------------------------------------

    use modfields,   only : svp
    use modglobal,   only : i1, j1, ih, jh, nsv, &
                            rdt, rtimee, rk3step, &
                            dzf, dx, dy    
    use modfields,   only : rhof
    use moddatetime, only : datex, nextday
    
    implicit none

    integer         :: k, l, count, n, i
    
    real            :: emistime_s, emistime_e ! Emission timers
    real, parameter :: div3600 = 1./3600.     ! Quick division

	integer, dimension(:), allocatable :: new_nullemiss
    ! Create a logical mask for elements not in nullemiss
  	logical :: is_in_nullemiss(nsv)

 
    if (.not. (l_emission)) return
 
    ! --------------------------------------------------------------------------
    ! Interpolate and apply emission
    ! --------------------------------------------------------------------------
    emistime_s = mod(rtimee +       1800., 3600.)*div3600

    ! MdB NOTE : Better way to do this? Problem is the broadcasting of 1D arrays
    ! rhof and dzf to emisfield. For now, loop over k.

   
   ! Create a logical array indicating if each 'l' is in nullemiss
	is_in_nullemiss = .false.
	
	! Count the number of elements greater than or equal to 1
  count = 0
  do i = 1, size(nullemiss)
    if (nullemiss(i) >= 1) then
      count = count + 1
    end if
  end do

  ! Allocate the new array
  allocate(new_nullemiss(count))

  ! Copy elements greater than or equal to 1 to the new array
  count = 0
  do i = 1, size(nullemiss)
    if (nullemiss(i) >= 1) then
      count = count + 1
      new_nullemiss(count) = nullemiss(i)
    end if
  end do
	
	
	is_in_nullemiss(new_nullemiss) = .true.
	

	do k = 1, kemis
  		do l = svskip + 1, nsv
    		if (.not. is_in_nullemiss(l)) then
    		
      			svp(2:i1, 2:j1, k, l) = svp(2:i1, 2:j1, k, l) + &
        			((1.0 - emistime_s) * emisfield(2:i1, 2:j1, k, l, 1) + &
        			emistime_s * emisfield(2:i1, 2:j1, k, l, 2)) / &
        			(3600.0 * rhof(k) * dzf(k) * dx * dy * 1.0e-6)
        			
   			 end if
  		end do
	end do


    ! -----
    ! Point sources
    ! ----

    !Intra-hour interpolation is applied
    
    if ( l_points .and. (npoints > 0) )  then
      call applypoints
    end if

    ! --------------------------------------------------------------------------
    ! Read emission files when neccesary, i.e. simulation reaches half hour mark
    ! after current timestep
    ! --------------------------------------------------------------------------

    if ( rk3step == 3 ) then
      emistime_e = mod(rtimee + rdt + 1800., 3600.)*div3600

      if ( emistime_e < emistime_s ) then
        ! Transfer data from 'ahead-of-modeltime' field to 'past-modeltime' field
        emisfield(:,:,:,:,1) = emisfield(:,:,:,:,2)
       
        ! Read new 'ahead-of-modeltime' emission field
        if ( datex(4) == 23 ) then
          call reademission(nextday(1), nextday(2), nextday(3),          0, emisfield(:,:,:,:,2))
        else
          call reademission(  datex(1),   datex(2),   datex(3), datex(4)+1, emisfield(:,:,:,:,2))
        endif
      endif

      ! --------------------------------------------------------------------------
      ! Point sources, with inter-hour interpolation 
      ! --------------------------------------------------------------------------
      ! Check if you need to handle point sources
      
  if (l_points .and. emistime_e < emistime_s) then
    if (datex(4) == 23) then
      if (npoints == 0) then
        call inquirepoints(nextday(1), nextday(2), nextday(3), 0)
      elseif (npoints > 0) then
        ! Transfer data from 'ahead-of-modeltime' field to 'past-modeltime' field
        point_source_data(:,:,1) = point_source_data(:,:,2)
        call readpoints(nextday(1), nextday(2), nextday(3), 0, point_source_data(:,:,2))
      endif
    else
      if (npoints == 0) then
        call inquirepoints(datex(1), datex(2), datex(3), datex(4) + 1)
      elseif (npoints > 0) then
        ! Transfer data from 'ahead-of-modeltime' field to 'past-modeltime' field
        point_source_data(:,:,1) = point_source_data(:,:,2)
        call readpoints(datex(1), datex(2), datex(3), datex(4) + 1, point_source_data(:,:,2))
      endif
    endif
  endif
endif
    

  end subroutine emission

  ! --------------------------------------------------------------------------
  ! Cleanup after run.
  ! --------------------------------------------------------------------------
  subroutine exitemission
  
    implicit none

    if (.not. (l_emission)) return

    deallocate(emisfield)
    deallocate(co2fields)
    
    if (l_points .and. (npoints > 0.)) deallocate(point_source_data)
    

    !if (l_points .and. (npoints > 0.)) deallocate(point_source_data1)
    !if (l_points .and. (npoints > 0.)) deallocate(point_source_data2)
	
  end subroutine exitemission
 
  subroutine inquirepoints(iyear, imonth, iday, ihour)
    use netcdf
    use modmpi,      only : myidx, myidy
    use modglobal,   only : imax, jmax, dx, dy, nsv
    
    implicit none
   
    integer, intent(in)  :: iyear, imonth, iday, ihour
    integer :: ncid, ndimid

    logical :: points_exist

    character(80) :: pointsource_file = 'pointsources.____________.x___.y___.nc'

    write(pointsource_file(14:17), '(i4.4)') iyear
    write(pointsource_file(18:19), '(i2.2)') imonth
    write(pointsource_file(20:21), '(i2.2)') iday
    write(pointsource_file(22:23), '(i2.2)') ihour
    write(pointsource_file(24:25), '(i2.2)') 0 !minutes
    write(pointsource_file(28:30), '(i3.3)') myidx
    write(pointsource_file(33:35), '(i3.3)') myidy

    inquire( file=pointsource_file, exist=points_exist )

    if ( points_exist ) then
      ! Read point source emission file
      call check( nf90_open(pointsource_file, NF90_NOWRITE, ncid))

      ! Determine total amount of sources for initialisation
      call check( nf90_inq_dimid(ncid, "n", ndimid ))
      call check( nf90_inquire_dimension(ncid, ndimid, len = npoints) )
      
      call check( nf90_close( ncid ) )

      write(6,*) pointsource_file, npoints
    else
      write(6,*) pointsource_file, 'NONE'
    end if

   contains

    subroutine check(status)
      integer, intent(in) :: status

      if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop 'NetCDF error in modemission inquirepoints. See outputfile for more information.'
      end if
    end subroutine check

  end subroutine inquirepoints

  
  subroutine readpoints(iyear, imonth, iday, ihour, point_source_data)
    use netcdf
    use modmpi, only: myidx, myidy
    use modglobal, only: imax, jmax, dx, dy, nsv

    implicit none

    integer, intent(in) :: iyear, imonth, iday, ihour
    real, intent(out) :: point_source_data(npoints, 7) 

    integer :: ncid, varid, tstep
    character(80) :: pointsource_file

   
        pointsource_file = 'pointsources.____________.x___.y___.nc'

        write(pointsource_file(14:17), '(i4.4)') iyear
        write(pointsource_file(18:19), '(i2.2)') imonth
        write(pointsource_file(20:21), '(i2.2)') iday
        write(pointsource_file(22:23), '(i2.2)') ihour
        write(pointsource_file(24:25), '(i2.2)') 0 !minutes
        write(pointsource_file(28:30), '(i3.3)') myidx
        write(pointsource_file(33:35), '(i3.3)') myidy

        ! Read point source emission file
        call check(nf90_open(pointsource_file, NF90_NOWRITE, ncid))

        ! Load data
        call check(nf90_inq_varid(ncid, "x_idx", varid))
        call check(nf90_get_var(ncid, varid, point_source_data(1:npoints, 1), count = (/npoints, 1/)))
        call check(nf90_inq_varid(ncid, "y_idx", varid))
        call check(nf90_get_var(ncid, varid, point_source_data(1:npoints, 2), count = (/npoints, 1/)))
        call check( nf90_inq_varid(ncid, "height", varid))       
    	call check( nf90_get_var  (ncid, varid, point_source_data(1:npoints, 3), count = (/npoints,1/) ) )
    	call check( nf90_inq_varid(ncid, "temperature", varid))       
    	call check( nf90_get_var  (ncid, varid, point_source_data(1:npoints, 4), count = (/npoints,1/) ) )
    	call check( nf90_inq_varid(ncid, "volume", varid))       
    	call check( nf90_get_var  (ncid, varid,	point_source_data(1:npoints, 5), count = (/npoints,1/) ) )
    	call check( nf90_inq_varid(ncid, "emission", varid))       
    	call check( nf90_get_var  (ncid, varid, point_source_data(1:npoints, 6), count = (/npoints,1/) ) )
    	call check( nf90_inq_varid(ncid, "tracer_idx", varid))       
    	call check( nf90_get_var  (ncid, varid,point_source_data(1:npoints, 7), count = (/npoints,1/) ) )


        call check(nf90_close(ncid))

contains

    subroutine check(status)
        integer, intent(in) :: status

        if (status /= nf90_noerr) then
            print *, trim(nf90_strerror(status))
            stop 'NetCDF error in modemission point sources. See outputfile for more information.'
        end if
    end subroutine check

end subroutine readpoints
  


  subroutine applypoints

    use modfields, only : svp, u0, v0, tmp0, rhof
    use modglobal, only : kmax, dx, dy, dzf,rdt, rtimee, nsv

    implicit none
	
    integer :: ipoint, ix, iy, iz, isv, izt, izb, l, iheight, i, count
    real    :: emis_b,emis_a, emis_int
    real    :: plume_top_fraction, plume_bottom_fraction, plumefactor
	
	real            :: emistime_s, emistime_e ! Emission timers
    real, parameter :: div3600 = 1./3600.     ! Quick division
    
	! Create a logical mask for elements not in nullemiss
  	logical :: is_in_nullemiss(nsv)
  	
  	integer, dimension(:), allocatable :: new_nullemiss2
    
    ! Create a logical array indicating if each 'l' is in nullemiss
	is_in_nullemiss = .false.
	
	
	! Count the number of elements greater than or equal to 1
  count = 0
  do i = 1, size(nullemiss)
    if (nullemiss(i) >= 1) then
      count = count + 1
    end if
  end do

  ! Allocate the new array
  allocate(new_nullemiss2(count))

  ! Copy elements greater than or equal to 1 to the new array
  count = 0
  do i = 1, size(nullemiss)
    if (nullemiss(i) >= 1) then
      count = count + 1
      new_nullemiss2(count) = nullemiss(i)
    end if
  end do
  
	is_in_nullemiss(new_nullemiss2) = .true.
	
    do ipoint = 1, npoints
        
       
	   
	   ix  =  int(point_source_data(ipoint,1,1)+0.1) !TODO change this! +0.1 is to prevent the conversion to result in an integer that's one too low
       iy  =  int(point_source_data(ipoint,2,1)+0.1) 
       !isv =  int(point_source_data(ipoint,7,1)+0.1)
       emis_b =     point_source_data(ipoint,6,1)  !point source for the 'past-modeltime' step
	   emis_a =     point_source_data(ipoint,6,2)  !point source for the 'ahead-of-modeltime' step
	   

       call briggs( tmp0(ix,iy,1:kmax), &                              !temperature profile at point source 
                    sqrt(v0(ix,iy,1:kmax)**2 + u0(ix,iy,1:kmax)**2), & !windspeed profile at point source

                    point_source_data(ipoint,4,1), &       ! Source temperature
                    point_source_data(ipoint,5,1), &       ! Source volume flux
                    point_source_data(ipoint,3,1), &       ! Source stack height

                    izt, plume_top_fraction, &  ! This index should be in full levels
                    izb, plume_bottom_fraction) ! plume top/bottom include the partially filled levels
       
       
       	
       ! --------------------------------------------------------------------------
	   ! Interpolate emission (now, the temporal interpolation is in the same way as for area emissions)
	   ! --------------------------------------------------------------------------
		
		
		emistime_s = mod(rtimee +       1800., 3600.)*div3600
	   
		emis_int=((1. - emistime_s)*emis_b + emistime_s *emis_a)
       
       !-------------------------------------------------------------------------------------------------
       
       ! Emissions are per source, so refactor to emission per gridbox
       ! ALSO: Emissions are per source, per hour so refactor to account for pressure, gridboxsize and seconds below:
       
       
       if (izt - izb > 1) then
       		 plumefactor =  1/(izt - izb + plume_top_fraction + plume_bottom_fraction - 1)
       end if
       
	    
	   ! Apply interpolated-in-time point source emissions:
		!-------------------------------------------------------------------------------------------------
       
       !emis_int is divided by number of vertical layers from plume_bottom to plume_top
       ! to put fraction of emission at each layer!
	   
	  do l = svskip+1, nsv
    	
    	if (.not. is_in_nullemiss(l)) then
  			if (izb == izt) then
    			svp(ix, iy, izb, l) = svp(ix, iy, izb, l) + emis_int / (3600. * rhof(izb) * dzf(izb) * dx * dy * 1e-6)
  			else if (izt - izb == 1) then
    			svp(ix, iy, izb, l) = svp(ix, iy, izb, l) + (emis_int/2) / (3600. * rhof(izb) * dzf(izb) * dx * dy * 1e-6)
    			svp(ix, iy, izt, l) = svp(ix, iy, izt, l) + (emis_int/2) / (3600. * rhof(izt) * dzf(izt) * dx * dy * 1e-6)
  			else if (izt - izb > 1) then
    			svp(ix, iy, izb, l) = svp(ix, iy, izb, l) + emis_int * plumefactor * plume_bottom_fraction / (3600. * rhof(izb) * dzf(izb) * dx * dy * 1e-6)
    			svp(ix, iy, izt, l) = svp(ix, iy, izt, l) + emis_int * plumefactor * plume_top_fraction / (3600. * rhof(izt) * dzf(izt) * dx * dy * 1e-6)
    			svp(ix, iy, izb+1:izt-1, l) = svp(ix, iy, izb+1:izt-1, l) + emis_int * plumefactor / (3600. * rhof(izb+1:izt-1) * dzf(izb+1:izt-1) * dx * dy * 1e-6)
  			end if
		end if
	
	enddo
  enddo

  end subroutine applypoints

  subroutine briggs(Ta, u, Ts, Vs, hs, iztop, ztop_frac, izbottom, zbottom_frac)
  
	!Algorithm to calculate the vertical plume rise above the stack height
	!The detail description can be found in Gordon et al., (2017) and Akingunola et al., (2018)

    use modglobal,   only : zh, zf, dzf, kmax
    !----------
    ! Ta Atmospheric temperature
    ! u  Wind speed
    ! Ts Emission temperature
    ! Vs Emission volume flux
    ! hs Emission stack height
    ! tzh Atmospheric temperature at half-level grid
    ! uzh Wind speed at half-level grid
    ! ths Atmospheric temperature at stack height
    ! uhs Wind speed at stack height
    ! hmax plume rise height (calculated relative to hs)
    ! S the stability parameter
    ! Fb buoyancy flux at stack height
    ! F1 Residual buoyancy flux ( iz +1)
    ! F0 Residual buoyancy flux ( iz)
    ! F0_old Residual buoyancy flux ( iz-1)
    ! iz half-grid index starts from the top of "stack" layer
    ! ---------

    implicit none

    real, intent(in)  :: Ts, hs, vs
    real, dimension(kmax), intent(in) :: Ta, u
    integer, intent(out) :: iztop, izbottom
    real,    intent(out) :: ztop_frac, zbottom_frac

    integer :: iz, i, imax, imin, ieq, iz0
    real    :: F0, F0_old, F1, Fb, dT, dU, ths, uhs
    real    :: ztop, zbottom
	real             :: gradT, S, hmax, upperh, lowerh, lowerw
    real, parameter  :: g = 9.81, &
                        cp = 1005., &
						pi = 3.1415926535897932
						
	real, dimension(kmax+1) :: tzh, uzh

    
	!Ta and U at half-level grid using simple avaraging:
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	do i=1,kmax
	
		if(i.eq.kmax) then
		
			dT= Ta(i) - tzh(i)
			dU= u(i) - uzh(i)
		else 
			dT= (Ta(i + 1) - Ta(i))/2.
			dU= (u(i + 1) - u(i))/2.
		
		endif 
		
		if((i==1).OR.(i==kmax)) then
		
			tzh(i) = Ta(i)-dT
			tzh(i+1) = Ta(i)+dT
			
			uzh(i) = u(i)-dU
			uzh(i+1) = u(i)+dU
				
		else
		
			tzh(i+1) = Ta(i)+dT
			uzh(i+1) = u(i)+dU
			
		endif	
	
	enddo
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!The atmospheric temperature (ths) and windspeed (uhs) at the stack height:
	
	imax = minloc(zh, dim =1, mask = zh.gt.hs)
	imin = maxloc(zh, dim =1, mask = zh.lt.hs)
	
	ieq = 0
	ieq = findloc(zh, hs, dim = 1)
	
	if(ieq.eq.0) then
	
	ths = tzh(imin)+(((tzh(imax)-tzh(imin))/(zh(imax)-zh(imin)))*(hs-zh(imin)))
	uhs = uzh(imin)+(((uzh(imax)-uzh(imin))/(zh(imax)-zh(imin)))*(hs-zh(imin)))
	
	else
	
	ths = tzh(ieq)
	uhs = uzh(ieq) 
	
	endif
	
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!The buoyancy flux at stack height:
	
	
	Fb = 0.
	F1 = 0.
	
	if (Ts.gt.ths)  then
	
		Fb = (g/pi)*Vs*(Ts-ths)/Ts
    
	endif
	
	F1 = Fb
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!Residual buoyancy flux(F1) and the exact plume rise height(hmax):
	
	iz = minloc(zh, dim=1, mask=zh.gt.hs)
	iz0 = iz
	
	
	S = 0.0

    do while ((F1 > 0.).AND.(iz<=kmax))
	

		
		if((iz.eq.iz0)) then
		
			gradT = (tzh(iz)-ths)/(zh(iz)-hs)
			S = g/tzh(iz)*(gradT + g/cp)
			
			F0_old = Fb  !Should F0_old and F0  be the same for the first step
			F0 = Fb      !since iz0 and hs are in the same layer?
			
			upperh = zh(iz)-hs
			lowerh = 0   !i.e., vertical distances are relative to the top of the stack
			lowerw = (uhs+uzh(iz))/2. !the mean windspeed (as in Gordon et al., 2017)
			
		
		elseif(iz.gt.iz0) then
			
			gradT = (tzh(iz)-tzh(iz-1))/(zh(iz)-zh(iz-1))
			S = g/tzh(iz)*(gradT + g/cp)
			F0_old = F0
			F0 = F1
			upperh = zh(iz)-hs  
			lowerh = zh(iz-1)-hs  
			lowerw = (uzh(iz-1)+uzh(iz))/2. !the mean windspeed (as in Gordon et al., 2017)
		
		endif
		
		
		
			
		if(S>=0.) then	
				
			F1 = min( F0 - 0.015*S*F0_old**(1./3.)*(upperh**(8./3.)-lowerh**(8./3.)), &
                         F0 - 0.053*S*lowerw*        (upperh**3.     -lowerh**3.) ) 
				
		else
			
			F1 = F0
					
		endif
			
		if(F1==0) then
			
			hmax = upperh
			
		endif
		
		if(F1<0) then
			
			hmax =  min(((F0/(0.015*S*F0_old**(1./3.))**(3./8.))+lowerh), &
						(((F0/(0.053*S*lowerw))**(1./3.))+lowerh)) 
		
		endif		  
		
		iz = iz + 1
		
		
    end do
	
  
    !Now, the parameterization can handle the exact plume rise height:
    
    ! Initialize zbottom and S to valid values
	zbottom = 0.0000001
	ztop = 0.0000001
 
    zbottom = hs - hmax*0.5            	


  ! Attempt to handle the problematic conditions
  if (isnan(zbottom) .or. isnan(S)) then
    ! Handle NaN values
    zbottom = 0.0000001
   
  else if (zbottom <= 0 .or. S < 0) then
    ! Handle non-positive values
    zbottom = 0.0000001
  endif

  
	ztop    = hs + hmax*1.5 
	
	if (isnan(ztop)) ztop = 0.0000001
	
	if(ztop<zbottom) ztop = zbottom
	
              
	
	izbottom = minloc(zh, dim=1, mask=zh.ge.zbottom)-1 !TODO makes sure this index is consistent with what is input in svp later
    iztop    = minloc(zh, dim=1, mask=zh.ge.ztop   )-1 !TODO makes sure this index is consistent with what is input in svp later

	
	
	ztop_frac    = (ztop           - zh(iztop))/dzf(iztop)    !(zh(iztop   +1) - zh(iztop   ))
    zbottom_frac = (zh(izbottom+1) - zbottom  )/dzf(izbottom) !(zh(izbottom+1) - zh(izbottom))
    

  end subroutine

  ! pure real function briggs_iter(F0, F1, z0, z1, T0, T1, U0)
    ! real, intent(in) :: F0, F1, z0, z1, T0, T1, U0
    ! real             :: S
    ! real, parameter  :: g = 9.81, &
                        ! cp = 1005.
    
    ! S = g/T1*(((T1-T0)/(z1-z0)) + g/cp)

    ! if (S > 0.) then
     
      ! briggs_iter = min( F0 - 0.015*S*F1**(1/3.)*(z1**(8./3.)-z0**(8./3.)), &
                         ! F0 - 0.053*S*U0*        (z1**3.     -z0**3.) ) 
    ! else
      ! briggs_iter = F0
    ! end if

  ! end function

end module modemission
