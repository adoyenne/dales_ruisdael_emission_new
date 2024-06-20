!
! Copyright (c) 2020-2020 Wageningen University and Research (WUR)
!
! This file is part of DALES
!
! DALES is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DALES.  If not, see <http://www.gnu.org/licenses/>.
!

module modlsm
    use netcdf
    implicit none

    public :: initlsm, lsm, exitlsm, init_lsm_tiles, lags

    logical :: llsm            ! On/off switch LSM
    logical :: lfreedrainage   ! Free drainage bottom BC for soil moisture
    logical :: lags            ! Switch for A-Gs scheme

    ! Interpolation types soil from full to half level
    integer :: iinterp_t, iinterp_theta
    integer, parameter :: iinterp_amean = 1  ! val = 0.5*(v1+v2)
    integer, parameter :: iinterp_gmean = 2  ! val = sqrt(v1*v2)
    integer, parameter :: iinterp_hmean = 3  ! val = ((dz1+dz2)*v1*v2)/(dz1*v2+dz2*v1)
    integer, parameter :: iinterp_max   = 4  ! val = max(a,b)

    ! Soil grid
    integer :: kmax_soil
    real :: z_size_soil
    real, allocatable :: z_soil(:), zh_soil(:)
    real, allocatable :: dz_soil(:), dzh_soil(:)
    real, allocatable :: dzi_soil(:), dzhi_soil(:)

    ! Soil index in `van_genuchten_parameters.nc` lookup table:
    integer, allocatable :: soil_index(:,:,:)

    ! Source term in soil water budget
    real, allocatable :: phiw_source(:,:,:)

    ! Precipitation, interception et al.
    real, allocatable :: throughfall(:,:)
    real, allocatable :: interception(:,:)
    real, allocatable :: wl_max(:,:)

    ! Dependency factor canopy resistance on VPD (high veg only)
    real, allocatable :: gD(:,:)
    

    ! Reduction functions canopy resistance
    real, allocatable :: f1(:,:), f2_lv(:,:), f2_hv(:,:), f2b(:,:), f3(:,:)

    ! Random
    real, allocatable :: du_tot(:,:), thv_1(:,:), land_frac(:,:)

    ! A-gs:
    
    ! A-gs vegetation type/index
    integer, allocatable :: index_ags(:,:)
    
    !net assimilation and soil respiration:
    real, allocatable :: an_co2(:,:), resp_co2(:,:)
    integer :: co2_index = -1

    ! Data structure for sub-grid tiles
    type lsm_tile
        ! Static properties:
        real, allocatable :: z0m(:,:), z0h(:,:)
        ! Base tile fraction (i.e. without liquid water)
        real, allocatable :: base_frac(:,:)
        ! Dynamic tile fraction:
        real, allocatable :: frac(:,:)
        ! Monin-obukhov / surface layer:
        real, allocatable :: obuk(:,:), ustar(:,:), ra(:,:)
        ! Conductivity skin layer:
        real, allocatable :: lambda_stable(:,:), lambda_unstable(:,:)
        ! Surface fluxes:
        real, allocatable :: H(:,:), LE(:,:), G(:,:)
        real, allocatable :: wthl(:,:), wqt(:,:)
        ! Surface (potential) temperature and humidity:
        real, allocatable :: tskin(:,:), thlskin(:,:), qtskin(:,:)
        ! Buoyancy difference surface - atmosphere
        real, allocatable :: db(:,:)
        ! Vegetation properties:
        real, allocatable :: lai(:,:), rs_min(:,:), rs(:,:)
        ! Root fraction in soil
        real, allocatable :: a_r(:,:), b_r(:,:)
        real, allocatable :: root_frac(:,:,:)
        ! Root fraction weighted mean soil water content
        real, allocatable :: phiw_mean(:,:)
    end type lsm_tile

    ! Tiles for low veg (lv), high veg (hv), bare soil (bs), wet skin (ws), water (aq):
    type(lsm_tile) tile_lv, tile_hv, tile_bs, tile_ws, tile_aq

    ! Land-surface / van Genuchten parameters from NetCDF input table.
    real, allocatable :: &
        theta_res(:), theta_wp(:), theta_fc(:), theta_sat(:), &
        gamma_theta_sat(:), vg_a(:), vg_l(:), vg_n(:)

    ! Derived soil parameters
    real, allocatable :: &
        vg_m(:), &
        lambda_theta_min(:), lambda_theta_max(:), &
        gamma_theta_min(:), gamma_theta_max(:), &
        gamma_t_dry(:), rho_C(:)

contains

subroutine lsm
    implicit none

    if (.not. llsm) return

    ! Calculate dynamic tile fractions,
    ! based on the amount of liquid water on vegetation.
    call calc_tile_fractions

    ! Calculate root fraction weighted mean soil water content.
    call calc_theta_mean(tile_lv)
    call calc_theta_mean(tile_hv)

    ! Calculate canopy/soil resistances.
    if (lags) then
        call calc_canopy_resistance_ags
    else
        call calc_canopy_resistance_js
    endif

    ! Calculate aerodynamic resistance (and u*, obuk).
    call calc_stability

    ! Set grid point averaged boundary conditions (thls, qts, gradients, ..)
    call calc_bulk_bcs

    ! Calculate soil tendencies
    ! Calc diffusivity heat:
    call calc_thermal_properties
    ! Solve diffusion equation:
    call integrate_t_soil

    ! Calc diffusivity and conductivity soil moisture:
    call calc_hydraulic_properties
    ! Calculate tendency due to root water extraction
    call calc_root_water_extraction
    ! Solve diffusion equation:
    call integrate_theta_soil

    ! Update liquid water reservoir
    call calc_liquid_reservoir

end subroutine lsm

!
! Calculate dynamic tile fractions, based on liquid water on vegetation
!
subroutine calc_tile_fractions
    use modglobal, only : i1, j1
    use modsurfdata, only : wl, wmax
    implicit none

    integer :: i, j
    real :: c_liq

    ! Water tile fraction is not variable
    tile_aq%frac(:,:) = tile_aq%base_frac(:,:)

    do j=2, j1
        do i=2, i1
            c_liq = min(1., wl(i,j)/wl_max(i,j))
            tile_ws%frac(i,j) = land_frac(i,j) * c_liq
            tile_lv%frac(i,j) = (1.-c_liq) * tile_lv%base_frac(i,j)
            tile_hv%frac(i,j) = (1.-c_liq) * tile_hv%base_frac(i,j)
            tile_bs%frac(i,j) = (1.-c_liq) * (1.-tile_lv%base_frac(i,j)-tile_hv%base_frac(i,j)-tile_aq%base_frac(i,j))
        end do
    end do

end subroutine calc_tile_fractions

!
! Calculate changes in the liquid water reservoir
!
subroutine calc_liquid_reservoir
    use modglobal,    only : rk3step, rdt, i1, j1, rhow, rlv
    use modsurfdata,  only : wl, wlm
    use modmicrodata, only : imicro, sed_qr
    implicit none

    integer :: i, j
    real :: tend, rk3coef, rainrate, c_veg, wl_tend_max, wl_tend_min
    real :: wl_tend_liq, wl_tend_dew, wl_tend_precip, wl_tend_sum, wl_tend_lim

    real, parameter :: intercept_eff = 0.5
    real, parameter :: to_ms  = 1./(rhow*rlv)

    rk3coef = rdt / (4. - dble(rk3step))
    if(rk3step == 1) wlm(:,:) = wl(:,:)

    do j=2, j1
        do i=2, i1
            c_veg = tile_lv%base_frac(i,j)+tile_hv%base_frac(i,j)

            ! Max and min possible tendencies
            wl_tend_min = -wlm(i,j) / rk3coef
            wl_tend_max = (wl_max(i,j) - wlm(i,j)) / rk3coef

            ! Tendency due to evaporation from liquid water reservoir/tile.
            wl_tend_liq = -max(0., tile_ws%frac(i,j) * tile_ws%LE(i,j) * to_ms)

            ! Tendency due to dewfall into vegetation/soil/liquid water tiles
            wl_tend_dew = &
                -( min(0., tile_lv%frac(i,j) * tile_lv%LE(i,j) * to_ms) + &
                   min(0., tile_hv%frac(i,j) * tile_hv%LE(i,j) * to_ms) + &
                   min(0., tile_bs%frac(i,j) * tile_bs%LE(i,j) * to_ms) + &
                   min(0., tile_ws%frac(i,j) * tile_ws%LE(i,j) * to_ms) )

            ! Tendency due to interception of precipitation by vegetation
            if (imicro == 0) then
                rainrate = 0.
            else
                rainrate = -sed_qr(i,j,1)/rhow
            end if

            wl_tend_precip = intercept_eff * c_veg * rainrate

            ! Total and limited tendencies
            wl_tend_sum = wl_tend_liq + wl_tend_dew + wl_tend_precip;
            wl_tend_lim = min(wl_tend_max, max(wl_tend_min,  wl_tend_sum));

            ! Diagnose interception and throughfall
            throughfall(i,j) = &
                -(1.-c_veg) * rainrate &
                -(1.-intercept_eff) * c_veg * rainrate &
                + min(0., wl_tend_lim - wl_tend_sum)
            interception(i,j) = max(0., wl_tend_lim)

            ! Integrate
            wl(i,j) = wlm(i,j) + rk3coef * wl_tend_lim

        end do
    end do

end subroutine calc_liquid_reservoir


!
! Calculate root fraction weighted mean soil water content
!
subroutine calc_theta_mean(tile)
    use modglobal,   only : i1, j1
    use modsurfdata, only : phiw
    implicit none

    type(lsm_tile), intent(inout) :: tile
    integer :: i, j, k, si
    real :: theta_lim

    tile%phiw_mean(:,:) = 0.

    do k=1, kmax_soil
        do j=2,j1
            do i=2,i1
                
                
                
                if((phiw(i,j,k) <= 1e-5)) phiw(i,j,k)=0.  !if warmstart, 
    											!prognostic phiw(i,j,k) might turns into negative values ( from the previous run, so bag somewhere in the code!!!)
    											!causing overflow and seg. fault, so  
    											!I set soil water content equals 0 if it is negative  !Arseni
    			    
                
                si = soil_index(i,j,k)

                theta_lim = max(phiw(i,j,k), theta_wp(si))
                tile%phiw_mean(i,j) = tile%phiw_mean(i,j) + tile%root_frac(i,j,k) * &
                    (theta_lim - theta_wp(si)) / (theta_fc(si) - theta_wp(si))
            end do
        end do
    end do

end subroutine calc_theta_mean

!
! Calculate canopy and soil resistances using Jarvis-Stewart method.
!
subroutine calc_canopy_resistance_js
    use modglobal,   only : i1, j1
    use modfields,   only : thl0, qt0, exnf, presf
    use modsurface,  only : ps
    use modraddata,  only : swd
    use modsurfdata, only : phiw
    implicit none

    integer :: i, j, k, si
    real :: swd_pos, T, esat, e, theta_min, theta_rel, c_veg

    ! Constants f1 calculation:
    real, parameter :: a_f1 = 0.81
    real, parameter :: b_f1 = 0.004
    real, parameter :: c_f1 = 0.05

    k = kmax_soil
    do j=2,j1
        do i=2,i1
            si = soil_index(i,j,k)

            ! f1: reduction vegetation resistance as f(sw_in):
            swd_pos = max(0., -swd(i,j,1))
            f1(i,j) = 1./min(1., (b_f1*swd_pos + c_f1) / (a_f1 * (b_f1*swd_pos + 1.)))

            ! f2: reduction vegetation resistance as f(theta):
            f2_lv(i,j) = 1./min(1., max(1.e-9, tile_lv%phiw_mean(i,j)))
            f2_hv(i,j) = 1./min(1., max(1.e-9, tile_hv%phiw_mean(i,j)))

            ! f3: reduction vegetation resistance as f(VPD) (high veg only):
            T    = thl0(i,j,1) * exnf(1)
            esat = 0.611e3 * exp(17.2694 * (T - 273.16) / (T - 35.86))
            e    = qt0(i,j,1) * presf(1) / 0.622

            f3(i,j) = 1./exp(-gD(i,j) * (esat-e))

            ! f2b: reduction soil resistance as f(theta)
            c_veg     = tile_lv%base_frac(i,j) + tile_hv%base_frac(i,j)
            theta_min = c_veg * theta_wp(si) + (1.-c_veg) * theta_res(si);
            theta_rel = (phiw(i,j,k) - theta_min) / (theta_fc(si) - theta_min);
            f2b(i,j)  = 1./min(1., max(1.e-9, theta_rel))

            ! Calculate canopy and soil resistance
            tile_lv%rs(i,j) = tile_lv%rs_min(i,j) / tile_lv%lai(i,j) * f1(i,j) * f2_lv(i,j)
            tile_hv%rs(i,j) = tile_hv%rs_min(i,j) / tile_hv%lai(i,j) * f1(i,j) * f2_hv(i,j) * f3(i,j)
            tile_bs%rs(i,j) = tile_bs%rs_min(i,j) / f2b(i,j)
            tile_ws%rs(i,j) = 0.

        end do
    end do

end subroutine calc_canopy_resistance_js


subroutine calc_canopy_resistance_ags
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate canopy and soil resistances using A-Gs (plant physiology).
! In addition, this calculates/sets the net CO2 fluxes...
! https://studenttheses.uu.nl/bitstream/handle/20.500.12932/31851/5764149,%20E.%20Merritt%20Master%27s%20Thesis%20SD.pdf?sequence=3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use modglobal,   only : i1, j1
    use modfields,   only : rhof, exnh, qt0, presf, svm
    use modsurface,  only : E1, ps
    use modsurfdata, only : svflux, tskin, phiw, tsoil
    use modraddata,  only : swd
    use modemisdata, only : l_emission, svco2ags, svco2sum, svco2veg

    implicit none

    ! NOTE: these should become a namelist switches...
    logical, parameter :: lsplitleaf = .false.
    logical, parameter :: lrelaxgc = .false.
    

    real :: Ts, co2_comp, gm, fmin0, fmin, esatsurf, e, Ds, Dmax, cfrac, co2_abs, ci, to_ppb, from_ppb  
    real :: Ammax, fstr, Am, Rdark, PAR, alphac, AGSa1, Dstar, tempy, An, gc_inf, gcco2, fw, t_mean, th_mean, rs_co2
    real :: cveg, cland, theta_min, theta_rel
    integer :: i, j, k, si
    
    real, dimension(i1, j1) :: tskin_copy
    
    ! Define the vegetation type variable
    integer :: veg_type

    ! Fixed constants (** = same in DALES and IFS, !! = different in DALES and IFS)
    ! Now, some constants are put in array with two separate vegetation types: grassland (1) and forest (2): Arseni
    
    ! Define constants: 
    real :: Q10gm, Q10am, Q10lambda, T1gm, T2gm, T1Am, T2Am, gmin, ad, Kx, alpha0, R10, gm298, Ammax298, f0, co2_comp298
    
    
    real, dimension(2) :: Q10gm_array = [2.0, 2.0]       ! (**) Parameter to calculate the mesophyll conductance [-]
    real, dimension(2) :: Q10am_array = [2.0, 2.0]        ! (**) Parameter to calculate max primary productivity [-]
    real, dimension(2) :: Q10lambda_array = [1.5, 1.5]    ! (!!) Parameter to calculate the CO2 compensation concentration. (2 in IFS, 1.5 in DALES) [-]

	! Reference temperatures calculation mesophyll conductance:
    real, dimension(2) :: T1gm_array = [278, 278]         ! (**)
    real, dimension(2) :: T2gm_array = [301, 305]         ! (!!: IFS=309, DALES=301 (default), 309 (C4))

    ! Reference temperatues calculation max primary productivity:
    real, dimension(2) :: T1Am_array = [286, 281]        ! (!!: IFS=281, DALES=286 (C4))
    real, dimension(2) :: T2Am_array = [311, 311]         ! (**)


    real, parameter :: nuco2q = 1.6       ! Ratio molecular viscosity water to carbon dioxide
	real, dimension(2) :: gmin_array = [2.5e-4, 2.5e-4]      ! Cuticular (minimum) conductance. NOTE: = g_cu in IFS?
    real, dimension(2) :: ad_array =  [0.07, 0.07]         ! Regression coefficient to calculate Cfrac
    real, dimension(2) :: Kx_array = [0.7,  0.7]          ! Extinction coefficient PAR [mground m−1leaf]


    !????
    real, dimension(2) :: alpha0_array =  [0.014, 0.017]    ! Light use efficiency at low light conditions (mg J−1)


    real, parameter :: Mair = 28.97
    real, parameter :: Mco2 = 44.01

    ! Parameters respiration Jacobs (2006)
    real, parameter :: Cw = 1.6e-3        ! Constant water stress correction
    real, parameter :: wsmax = 0.55       ! Upper reference value soil water   # NOTE: I guess these are soil dependant?
    real, parameter :: wsmin = 0.005      ! Lower reference value soil water   # NOTE: see line above :-)
    real, dimension(2) :: R10_array =   [0.23, 0.1]       ! Respiration at 10oC (Jacobs 2007)
    real, parameter :: Eact0 = 53.3e3     ! Activation energy

    ! Vegetation specific constants
    ! TODO: link to lookup table IFS, for now hardcoded for single veg type.....
    ! Vegetation specific constants  
    real, dimension(2) :: gm298_array = [7.0, 3.0] ! Mesophyll conductance at 298 K [mm s-1] NOTE: Much lower than DALES...
    real, dimension(2) :: Ammax298_array = [1.7, 2.2] ! CO2 maximal primary productivity [mg CO2 m-2 s-1]
    real, dimension(2) :: f0_array = [0.85, 0.89]          ! Maximum value Cfrac (constant or equation in IFS)
    real, dimension(2) :: co2_comp298_array = [68.5, 68.5] ! CO2 compensation concentration = Lambda(25) in IFS [ppm]

    ! For sw_splitleaf
    !nr_gauss = 3                                 ! Amount of bins to use for Gaussian integrations
    !weight_g = np.array([0.2778,0.4444,0.2778])  ! Weights of the Gaussian bins (must add up to 1)
    !angle_g  = np.array([0.1127,   0.5,0.8873])  ! Sines of the leaf angles compared to the sun in the first Gaussian integration
    !LAI_g    = np.array([0.1127,   0.5,0.8873])  ! Ratio of integrated LAI at locations where shaded leaves are evaluated in the second Gaussian integration
    !sigma    = 0.2                               ! Scattering coefficient
    !kdfbl    = 0.8                               ! Diffuse radiation extinction coefficient for black leaves

    ! Conversion factors 
    ! NOTE: this should use the surface density, not `rhof(1)`.
    to_ppb   = Mair/Mco2/rhof(1)*1000
    from_ppb = 1./to_ppb
    

    do j=2,j1
        do i=2,i1
       	 	
       	 	tskin_copy(i,j)=tskin(i, j)
       	 	
       	 	!Get the vegetation type for current grid point:
            veg_type = index_ags(i, j)

            !Check the veg_type and select the appropriate constants:
            if (veg_type == 1) then
            
                ! Grassland constants
                Q10gm = Q10gm_array(1)
                Q10am = Q10am_array(1)
                Q10lambda = Q10lambda_array(1)
                T1gm = T1gm_array(1)
                T2gm = T2gm_array(1)
                T1Am = T1Am_array(1)
                T2Am = T2Am_array(1)
                gmin = gmin_array(1)
                ad = ad_array(1)
                Kx = Kx_array(1)
                alpha0 = alpha0_array(1)
                R10 = R10_array(1)
                gm298 = gm298_array(1)
                Ammax298 = Ammax298_array(1)
                f0 = f0_array(1)
                co2_comp298 = co2_comp298_array(1)

            else if (veg_type == 2) then
            
                ! Forest constants
                Q10gm = Q10gm_array(2)
                Q10am = Q10am_array(2)
                Q10lambda = Q10lambda_array(2)
                T1gm = T1gm_array(2)
                T2gm = T2gm_array(2)
                T1Am = T1Am_array(2)
                T2Am = T2Am_array(2)
                gmin = gmin_array(2)
                ad = ad_array(2)
                Kx = Kx_array(2)
                alpha0 = alpha0_array(2)
                R10 = R10_array(2)
                gm298 = gm298_array(2)
                Ammax298 = Ammax298_array(2)
                f0 = f0_array(2)
                co2_comp298 = co2_comp298_array(2)
            end if


       	 	!Check if tskin(i, j) is NaN
        	if (isnan(tskin_copy(i, j))) then
        	
            	
            	! Surface resistances for moisture and carbon dioxide:
           		tile_lv%rs(i,j) = 0.
            	tile_hv%rs(i,j) = 0.
            	tile_bs%rs(i,j) = 0.
            	tile_ws%rs(i,j) = 0.
            	
            	!respiration and vegetation:
            	resp_co2(i,j) = 0.
            	an_co2(i,j) = 0.        	
            	        	
            	
       	 	else
       	 
  
    			
    			!Check if tskin_copy(i, j) less than 200 or more than 350, and replace with 287.
        		if ((tskin_copy(i, j).lt.200.).OR.(tskin(i, j).gt.350.)) then
            		tskin_copy(i, j) = 287.
       	 		endif
            
            	si = soil_index(i, j, kmax_soil)
			
				Ts = tskin_copy(i,j) * exnh(1)
			

            	! Calculate the CO2 compensation concentration (IFS eq. 8.92)
            	! "The compensation point Γ is defined as the CO2 concentration at which the net CO2 assimilation of a fully lit leaf becomes zero."
            	! NOTE: The old DALES LSM used the atmospheric `thl`, IFS uses the surface temperature.
            	co2_comp = rhof(1) * co2_comp298 * Q10lambda**(0.1 * (Ts - 298.0)) 

           	 	! Calculate the mesophyll conductance (IFS eq. 8.93)
            	! "The mesophyll conductance gm describes the transport of CO2 from the substomatal cavities to the mesophyll cells where the carbon is fixed."
            	! NOTE: The old DALES LSM used the atmospheric `thl`, IFS uses the surface temperature.
            	gm = gm298 * Q10gm**(0.1 * (Ts - 298.0)) / ((1. + exp(0.3 * (T1gm - Ts))) * (1. + exp(0.3 * (Ts - T2gm)))) / 1000.
	
            	! Calculate CO2 concentration inside the leaf (ci)
            	! NOTE: Differs from IFS
            	fmin0 = gmin/nuco2q - (1./9.)*gm
            	fmin  = (-fmin0 + (fmin0**2 + 4*gmin/nuco2q*gm)**0.5) / (2.*gm)

            	! Calculate atmospheric moisture deficit
           	 	! NOTE: "Therefore Ci/Cs is specified as a function of atmospheric moisture deficit Ds at the leaf surface"
            	!       In DALES, this uses a mix between esat(surface) and e(atmosphere)
            	!       In IFS, this uses (in kg kg-1): qsat(Ts)-qs instead of qsat(Ts)-qa!
            	! NOTE: Old DALES LSM used the surface pressure in the calculation of `e`, not sure why...
            	esatsurf = 0.611e3 * exp(17.2694 * (Ts - 273.16) / (Ts - 35.86))
            	e = qt0(i,j,1) * presf(1) / 0.622
            	Ds = (esatsurf - e) / 1000.
             
            
            	! This seems to differ from IFS?
            	Dmax = (f0 - fmin) / ad

            	! Coupling factor (IFS eq. 8.101)
            	!cfrac = f0 * (1.0 - Ds/Dmax) + fmin * (Ds/Dmax)
            	cfrac = max(0.01, f0 * (1.0 - Ds/Dmax) + fmin * (Ds/Dmax))

            	! Absolute CO2 concentration:
            	if (l_emission) then
            		co2_abs =  (svm(i,j,1,svco2sum)*1000*(Mair/Mco2))*from_ppb !(svm scalar tacers for gases are in ug/g now, so converted first to ppb and after to abs)
			    	
            	else
                	co2_abs =  (svm(i,j,1,co2_index)*1000*(Mair/Mco2))*from_ppb !(svm scalar tacers for gases are in ug/g now, so converted first to ppb and after to abs)

				endif
				

            !if (lrelaxci) then
            !  if (ci_old_set) then
            !    ci_inf        = cfrac * (co2_abs - co2_comp) + co2_comp
            !    ci            = ci_old(i,j) + min(kci*rk3coef, 1.0) * (ci_inf - ci_old(i,j))
            !    if (rk3step  == 3) then
            !      ci_old(i,j) = ci
            !    endif
            !  else
            !    ci            = cfrac * (co2_abs - co2_comp) + co2_comp
            !    ci_old(i,j)   = ci
            !  endif
            !else
            !  ci              = cfrac * (co2_abs - co2_comp) + co2_comp
            !endif

            	! CO2 concentration in leaf (IFS eq. ???): 
            	ci = cfrac * (co2_abs - co2_comp) + co2_comp

            	! Max gross primary production in high light conditions (Ag) (IFS eq. 8.94)
            	! NOTE: The old DALES LSM used the atmospheric `thl`, IFS uses the surface temperature.
            	Ammax = Ammax298 * Q10am**(0.1 * (Ts - 298.0)) / ((1.0 + exp(0.3 * (T1Am - Ts))) * (1. + exp(0.3 * (Ts - T2Am))))

            	! Effect of soil moisture stress on gross assimilation rate.
           	 	! NOTE: this seems to be different in IFS...
            	! NOTE: for now, this uses the relative soil moisture content from the low vegetation only!
            	!fstr = max(1.0e-3, min(1.0, tile_lv%phiw_mean(i,j)))
            	
            	!NOTE: this uses the relative soil moisture content for both low and high vegetation
            	! We check the veg_type and select the appropriate tile (for grass -> lv for forest -> hv):
            	if (veg_type == 1) then
            		fstr = max(1.0e-3, min(1.0, tile_lv%phiw_mean(i,j)))
            	else if (veg_type == 2) then
            		fstr = max(1.0e-3, min(1.0, tile_hv%phiw_mean(i,j)))
            	end if

            	! Gross assimilation rate (Am, IFS eq. 8.97)
            	Am = Ammax * (1 - exp(-(gm * (ci - co2_comp) / Ammax)))
            	
            	if((Am.lt.0.01).OR.(isnan(Am))) then	!I fix Am.lt.0.01 as threshold 
            											!since Am below could cause problems later in the An calculation !Arseni
            
            		! Surface resistances for moisture and carbon dioxide:
           			tile_lv%rs(i,j) = 0.
            		tile_hv%rs(i,j) = 0.
            		tile_bs%rs(i,j) = 0.
            		tile_ws%rs(i,j) = 0.
            	
            		!respiration and vegetation:
            		resp_co2(i,j) = 0.
            		an_co2(i,j) = 0.
            
           	 	else
            

            		! Autotrophic dark respiration (IFS eq. 8.99)
            		Rdark = Am/9.
            

            		!PAR = 0.40 * max(0.1,-swdav * cveg(i,j))
            		PAR = 0.5 * max(0.1, -swd(i,j,1))

            		! Light use efficiency
            		alphac = alpha0 * (co2_abs - co2_comp) / (co2_abs + 2*co2_comp)
                      
            
            		if (lsplitleaf) then
                		print*,'Splitleaf A-Gs not (yet) implemented!'
                		stop
            	!    PARdir   = 0.5 * max(0.1, sw_in[j,i])
            	!    PARdif   = 0.5 * max(0.1, sw_in[j,i])
            	!
            	!    cos_sza = max(1.e-10, calc_zenith_angle(lat, lon, doy, time))
            	!
            	!    kdrbl = 0.5 / cos_sza   ! Direct radiation extinction coefficient for black leaves
            	!    kdf = kdfbl * np.sqrt(1.0-sigma)
            	!    kdr = kdrbl * np.sqrt(1.0-sigma)
            	!    ref = (1.0 - np.sqrt(1.0-sigma)) / (1.0 + np.sqrt(1.0-sigma)) ! Reflection coefficient
            	!    ref_dir = 2 * ref / (1.0 + 1.6 * cos_sza)
            	!
            	!    Hleaf = np.zeros(nr_gauss+1)
            	!    Fleaf = np.zeros(nr_gauss+1)
            	!    Agl   = np.zeros(nr_gauss+1)
            	!
            	!    Fnet  = np.zeros(nr_gauss)
            	!    gnet  = np.zeros(nr_gauss)
            	!
            	!    ! Loop over different LAI locations
            	!    for it in range(nr_gauss):
            	!
            	!        iLAI   = LAI * LAI_g[it]        ! Integrated LAI between here and canopy top; Gaussian distributed
            	!        fSL    = np.exp(-kdrbl * iLAI)  ! Fraction of sun-lit leaves
            	!
            	!        PARdfD = PARdif * (1.0-ref)     * np.exp(-kdf * iLAI)     ! Total downward PAR due to diffuse radiation at canopy top
            	!        PARdrD = PARdir * (1.0-ref_dir) * np.exp(-kdr * iLAI)     ! Total downward PAR due to direct radiation at canopy top
            	!        PARdfU = PARdif * (1.0-ref)     * np.exp(-kdf * LAI) * \
            	!                albedo * (1.0-ref) * np.exp(-kdf * (LAI-iLAI))    ! Total upward (reflected) PAR that originates as diffuse radiation
            	!        PARdrU = PARdir * (1.0-ref_dir) * np.exp(-kdr * LAI) * \
            	!                albedo * (1.0-ref) * np.exp(-kdf * (LAI-iLAI))    ! Total upward (reflected) PAR that originates as direct radiation
            	!        PARdfT = PARdfD + PARdfU                                  ! Total PAR due to diffuse radiation at canopy top
            	!        PARdrT = PARdrD + PARdrU                                  ! Total PAR due to direct radiation at canopy top
            	!
            	!        dirPAR = (1.0-sigma) * PARdir * fSL                       ! Purely direct PAR (can only be downward)
            	!        difPAR = PARdfT + PARdrT - dirPAR                         ! Total diffuse radiation
            	!
            	!        HdfT   = kdf * PARdfD + kdf * PARdfU
            	!        HdrT   = kdr * PARdrD + kdf * PARdrU
            	!        dirH   = kdrbl * dirPAR
            	!        Hshad  = HdfT + HdrT - dirH
            	!
            	!        Hsun   = Hshad + angle_g * (1.0-sigma) * kdrbl * PARdir / np.sum(angle_g * weight_g)
            	!
            	!        Hleaf[0]  = Hshad
            	!        Hleaf[1:] = Hsun
            	!
            	!        Agl    = fstr * (Am + Rdark) * (1 - np.exp(-alphac*Hleaf/(Am + Rdark)))
            	!        gleaf  = gmin/nuco2q +  Agl/(co2_abs-ci)
            	!        Fleaf  = Agl - Rdark
            	!
            	!        Fshad  = Fleaf[0]
            	!        Fsun   = np.sum(weight_g * Fleaf[1:])
            	!        gshad  = gleaf[0]
            	!        gsun   = np.sum(weight_g * gleaf[1:])
            	!
            	!        Fnet[it] = Fsun * fSL + Fshad * (1 - fSL)
            	!        gnet[it] = gsun * fSL + gshad * (1 - fSL)
            	!
            	!    An     = LAI * np.sum(weight_g * Fnet)
            	!    gc_inf = LAI * np.sum(weight_g * gnet)
            	!
            		else
                		! Calculate upscaling from leaf to canopy: net flow CO2 into the plant (An)
                		! NOTE: this only uses LAI from low vegetation (for now..)
                		AGSa1  = 1.0 / (1 - f0)
                		Dstar  = Dmax / (AGSa1 * (f0 - fmin))
                		tempy  = alphac * Kx * PAR / (Am + Rdark)

    					!An = (Am + Rdark) * (1 - 1.0 / (Kx * tile_lv%lai(i, j)) * (E1(tempy * exp(-Kx * tile_lv%lai(i, j))) - E1(tempy)))
						
                		!!An = (Am + Rdark) * (1.0 - 1.0 / denominator * (E1(tempy * exp_term) - E1(tempy)))
                		!gc_inf = tile_lv%lai(i,j) * (gmin/nuco2q + AGSa1 * fstr * An / ((co2_abs - co2_comp) * (1 + Ds / Dstar)))
            		
            			
            			! NOTE: this uses LAI from both low and high vegetation:
            			! We check the veg_type and select the appropriate tile (for grass -> lv for forest -> hv)
            			if (veg_type == 1) then
            				
            				An = (Am + Rdark) * (1 - 1.0 / (Kx * tile_lv%lai(i, j)) * (E1(tempy * exp(-Kx * tile_lv%lai(i, j))) - E1(tempy)))
                			gc_inf = tile_lv%lai(i,j) * (gmin/nuco2q + AGSa1 * fstr * An / ((co2_abs - co2_comp) * (1 + Ds / Dstar)))
            		
            			else if (veg_type == 2) then
            				
            				An = (Am + Rdark) * (1 - 1.0 / (Kx * tile_hv%lai(i, j)) * (E1(tempy * exp(-Kx * tile_hv%lai(i, j))) - E1(tempy)))
                			gc_inf = tile_hv%lai(i,j) * (gmin/nuco2q + AGSa1 * fstr * An / ((co2_abs - co2_comp) * (1 + Ds / Dstar)))
            		
            			end if
            		
            		
            		
            		endif

            		if (lrelaxgc) then
               	 		print*,'Relax GC A-Gs not (yet) implemented!'
               	 		stop
            	!  if (gc_old_set) then
            	!    gcco2       = gc_old(i,j) + min(kgc*rk3coef, 1.0) * (gc_inf - gc_old(i,j))
            	!    if (rk3step ==3) then
            	!      gc_old(i,j) = gcco2
            	!    endif
            	!  else
            	!    gcco2 = gc_inf
            	!    gc_old(i,j) = gcco2
            	!  endif
            		else
                		gcco2 = gc_inf
            		endif

            		! Calculate mean t_soil and theta_soil for calculation soil respiration
            		t_mean = 0
            		th_mean = 0
            		k = kmax_soil
            		do k=1, kmax_soil
                		t_mean  = t_mean  + tsoil(i,j,k) * dz_soil(k)
                		th_mean = th_mean + phiw(i,j,k)  * dz_soil(k)
            		enddo
            		t_mean  = t_mean  / (-zh_soil(1))
            		th_mean = th_mean / (-zh_soil(1))

            		! Sub-grid fractions of low+high veg, and low+high veg + bare soil
            		!cveg = tile_lv%base_frac(i,j) + tile_hv%base_frac(i,j)
            		
            		!veg is selected based on veg_type:
            		!We check the veg_type and select the appropriate tile (for grass -> lv for forest -> hv)
            		if (veg_type == 1) then
            			cveg = tile_lv%base_frac(i,j)
            		else if (veg_type == 2) then
            			cveg = tile_hv%base_frac(i,j)
            		end if
            
            		
            		cland = cveg + tile_bs%base_frac(i,j)

            		! Water stress function:
            		fw = Cw * wsmax / (th_mean + wsmin)

            		! Soil resistance, using old/JS approach
            		theta_min = cveg * theta_wp(si) + (1.-cveg) * theta_res(si)
            		theta_rel = (phiw(i,j,kmax_soil) - theta_min) / (theta_fc(si) - theta_min);
            		f2b(i,j)  = 1./min(1., max(1.e-9, theta_rel))

            		! Surface resistances for moisture and carbon dioxide
            		tile_lv%rs(i,j) = 1.0 / (1.6 * gcco2)
            		tile_hv%rs(i,j) = 1.0 / (1.6 * gcco2)
            		tile_bs%rs(i,j) = tile_bs%rs_min(i,j) / f2b(i,j)
            		tile_ws%rs(i,j) = 0.
            		
            		!Fix to replace possible NaNs in canopy resistances:
            		if(isnan(tile_lv%rs(i,j))) tile_lv%rs(i,j)=100.  ! s m-1
            		if(isnan(tile_hv%rs(i,j))) tile_hv%rs(i,j)=100.  ! s m-1
            		if(isnan(tile_bs%rs(i,j))) tile_bs%rs(i,j)=100.  ! s m-1

            		rs_co2 = 1. / gcco2

            		! NOTE:Combined assimilation for low and high veg, using only low veg as vegetation type:
            		resp_co2(i,j) = R10 * (1.-fw) * exp(Eact0 / (283.15 * 8.314) * (1.0 - 283.15 / t_mean))
            		! Scale with vegetation fraction, and translate to `ug/g m s-1`.
			
					resp_co2(i,j) = cveg * resp_co2(i,j)*((Mco2/Mair)/1000) * to_ppb
					

            		! Calculate net assimilation
            		! NOTE: this uses the aerodynamic resistance from low vegetation...
            		!an_co2(i,j) = -(co2_abs - ci) / (tile_lv%ra(i,j) + rs_co2)
            		
            		!NOTE: this uses the aerodynamic resistance from both low and high vegetations:
            		! We check the veg_type and select the appropriate tile (for grass -> lv for forest -> hv)
            		if (veg_type == 1) then
            			an_co2(i,j) = -(co2_abs - ci) / (tile_lv%ra(i,j) + rs_co2)

            		else if (veg_type == 2) then
            			an_co2(i,j) = -(co2_abs - ci) / (tile_hv%ra(i,j) + rs_co2)
            		end if
            		
            		
            		! Scale with vegetation + bare soil fraction, and translate to `ug/g m s-1`.
            		an_co2(i,j) = cland * an_co2(i,j)*((Mco2/Mair)/1000) * to_ppb
            		
            		!Replace unreasonable values with 0:
            		if ((an_co2(i,j)>0.) .OR. (an_co2(i,j)<-150.)) an_co2(i,j)=0.
					if ((resp_co2(i,j)<0.) .OR. (resp_co2(i,j)>150.)) resp_co2(i,j)=0.
					
					
					!Only for the sanity check of forest vegetation, here, we nullify all others grid points:
					!if (veg_type .ne. 2) then
						!an_co2(i,j)=0
					!endif
					
						
				endif
			endif		
			
			! Set CO2 fluxes:
			
            if (l_emission) then
                ! Total flux into the CO2 field holding the sum of all CO2 concentrations:
                svflux(i,j,svco2sum) = resp_co2(i,j) + an_co2(i,j) !Kinematic scalar flux [ug/g m/s]
				
                ! Respiration flux into the respiration specific field:
                svflux(i,j,svco2ags) = resp_co2(i,j)  !Kinematic scalar flux [ug/g m/s]

                ! Net update into the vegetation specific field:
                !svflux(i,j,svco2veg) = an_co2(i,j) !Kinematic scalar flux [ug/g m/s]
            else
            
            	! Respiration flux into the respiration specific field:
                svflux(i,j,svco2ags) = resp_co2(i,j)  !Kinematic scalar flux [ug/g m/s]

                svflux(i,j,co2_index) = resp_co2(i,j) + an_co2(i,j)  !Kinematic scalar flux [ug/g m/s] 
            endif

        end do
    end do

end subroutine calc_canopy_resistance_ags

!
! Calculate Obukhov length, ustar, and aerodynamic resistance, for all tiles
!
subroutine calc_stability
    use modglobal, only : i1, j1, cu, cv, rv, rd
    use modfields, only : u0, v0, thl0, qt0
    implicit none

    real, parameter :: du_min = 0.1
    real :: du, dv
    integer :: i, j

    ! Calculate properties shared by all tiles:
    ! Absolute wind speed difference, and virtual potential temperature atmosphere
    do j=2,j1
        do i=2,i1
            du = 0.5*(u0(i,j,1) + u0(i+1,j,1)) + cu
            dv = 0.5*(v0(i,j,1) + v0(i,j+1,1)) + cv
            du_tot(i,j) = max(0.1, sqrt(du**2 + dv**2))

            thv_1(i,j) = thl0(i,j,1)  * (1.+(rv/rd-1.)*qt0(i,j,1))
        end do
    end do

    call calc_obuk_ustar_ra(tile_lv)
    call calc_obuk_ustar_ra(tile_hv)
    call calc_obuk_ustar_ra(tile_bs)
    call calc_obuk_ustar_ra(tile_ws)
    call calc_obuk_ustar_ra(tile_aq)

end subroutine calc_stability

!
! Calculate Obukhov length and ustar, for single tile
!
subroutine calc_obuk_ustar_ra(tile)
    use modglobal, only : i1, j1, rd, rv, grav, zf
    use modfields, only : u0, v0
    implicit none

    type(lsm_tile), intent(inout) :: tile
    integer :: i, j
    real :: thvs

    do j=2,j1
        do i=2,i1
            if (tile%frac(i,j) > 0) then
                ! Buoyancy difference surface - atmosphere
                thvs = tile%thlskin(i,j) * (1.+(rv/rd-1.)*tile%qtskin(i,j))
                !tile%db(i,j) = grav/thvs * (thvs - thv_1(i,j))
                tile%db(i,j) = grav/thvs * (thv_1(i,j) - thvs)

                ! Iteratively find Obukhov length
                tile%obuk(i,j) = calc_obuk_dirichlet( &
                    tile%obuk(i,j), du_tot(i,j), tile%db(i,j), zf(1), tile%z0m(i,j), tile%z0h(i,j))
            end if
        end do
    end do

    do j=2,j1
        do i=2,i1
            if (tile%frac(i,j) > 0) then
                ! Calculate friction velocity and aerodynamic resistance
                tile%ustar(i,j) = du_tot(i,j) * fm(zf(1), tile%z0m(i,j), tile%obuk(i,j))
                tile%ra(i,j)    = 1./(tile%ustar(i,j) * fh(zf(1), tile%z0h(i,j), tile%obuk(i,j)))
            end if
        end do
    end do

end subroutine calc_obuk_ustar_ra

!
! Calculate surface temperature for each tile, and calculate
! surface fluxes (H, LE, G0, wthl, wqt) and values (thlskin, qtskin)
!
subroutine calc_tile_bcs(tile)
    use modglobal,   only : i1, j1, cp, rlv, boltz
    use modfields,   only : exnh, exnf, presh, thl0, qt0, rhof
    use modraddata,  only : swd, swu, lwd, lwu
    use modsurfdata, only : ps, tsoil
    implicit none

    type(lsm_tile), intent(inout) :: tile
    integer :: i, j
    real :: Ts, thvs, esats, qsats, desatdTs, dqsatdTs, &
        rs_lim, fH, fLE, fG, num, denom, Ta, qsat_new, &
        rhocp_i, rholv_i, Qnet

    rhocp_i = 1. / (rhof(1) * cp)
    rholv_i = 1. / (rhof(1) * rlv)

    ! tmp
    !real :: HLEG, lwu_new, Qnet_new


    do j=2, j1
        do i=2, i1
            if (tile%frac(i,j) > 0) then

                ! Disable canopy resistance in case of dew fall
                

				Ts  = tile%thlskin(i,j) * exnh(1) 

                esats = 0.611e3 * exp(17.2694 * (Ts - 273.16) / (Ts - 35.86))
                qsats = 0.622 * esats / ps

                if (qsats < qt0(i,j,1)) then
                    rs_lim = 0.
                else
                    rs_lim = tile%rs(i,j)
                end if

                ! Calculate recuring terms
                ! NOTE: this should use the surface density, not rhof(1).
                !       Not sure if rhoh is available somewhere...
                fH  = rhof(1) * cp  / tile%ra(i,j)
                fLE = rhof(1) * rlv / (tile%ra(i,j) + rs_lim)

                if (tile%db(i,j) > 0) then
                    fG = tile%lambda_stable(i,j)
                else
                    fG = tile%lambda_unstable(i,j)
                end if

                ! Net radiation; negative sign = net input of radiation at surface
                Qnet = swd(i,j,1) + swu(i,j,1) + lwd(i,j,1) + lwu(i,j,1)

                ! Solve new skin temperature from SEB
                desatdTs = esats * (17.2694 / (Ts - 35.86) - 17.2694 * (Ts - 273.16) / (Ts - 35.86)**2.)
                dqsatdTs = 0.622 * desatdTs / ps
                Ta = thl0(i,j,1) * exnf(1)  

                num = -(Qnet - lwu(i,j,1) &
                      - fH * Ta + (qsats - dqsatdTs * Ts - qt0(i,j,1)) * fLE &
                      - fG * tsoil(i,j,kmax_soil) - 3.*boltz * Ts**4)
                denom = (fH + fLE * dqsatdTs + fG + 4.*boltz * Ts**3)
                tile%tskin(i,j) = num / denom

                ! Update qsat with linearised relation, to make sure that the SEB closes
                qsat_new = qsats + dqsatdTs * (tile%tskin(i,j) - Ts)

                ! Calculate energetic surface fluxes
                tile%H (i,j) = fH  * (tile%tskin(i,j) - Ta)
                tile%LE(i,j) = fLE * (qsat_new - qt0(i,j,1))
                tile%G (i,j) = fG  * (tsoil(i,j,kmax_soil) - tile%tskin(i,j))

                ! Calculate kinematic surface fluxes
                tile%wthl(i,j) = tile%H (i,j) * rhocp_i
                tile%wqt (i,j) = tile%LE(i,j) * rholv_i

                ! Calculate surface values
                tile%thlskin(i,j) = thl0(i,j,1) + tile%wthl(i,j) * tile%ra(i,j)
                tile%qtskin (i,j) = qt0 (i,j,1) + tile%wqt (i,j) * tile%ra(i,j)
            end if
        end do
    end do

end subroutine calc_tile_bcs

!
! Calculate BCs and fluxes for the water tile
!
subroutine calc_water_bcs
    use modglobal,   only : i1, j1, cp, rlv
    use modfields,   only : exnh, thl0, qt0, rhof
    use modsurfdata, only : ps
    implicit none

    integer :: i, j
    real :: esats

    do j=2, j1
        do i=2, i1
            if (tile_aq%frac(i,j) > 0) then
                ! Calculate BCs. `tile_aq%tskin` is fixed in time (for now...)
                tile_aq%thlskin(i,j) = tile_aq%tskin(i,j) / exnh(1)
                esats = 0.611e3 * exp(17.2694 * (tile_aq%tskin(i,j) - 273.16) / (tile_aq%tskin(i,j) - 35.86))
                tile_aq%qtskin(i,j) = 0.622 * esats / ps

                ! Calculate kinematic fluxes
                tile_aq%wthl(i,j) = 1./tile_aq%ra(i,j) * (tile_aq%thlskin(i,j) - thl0(i,j,1))
                tile_aq%wqt (i,j) = 1./tile_aq%ra(i,j) * (tile_aq%qtskin (i,j) - qt0 (i,j,1))

                ! Calculate energetic fluxes
                tile_aq%H (i,j) = tile_aq%wthl(i,j) * rhof(1) * cp
                tile_aq%LE(i,j) = tile_aq%wqt (i,j) * rhof(1) * rlv
                tile_aq%G (i,j) = 0.
            end if
        end do
    end do

end subroutine calc_water_bcs

!
! Set grid point mean boundary conditions, as used by
! the diffusion scheme, thermodynamics, ...
!
subroutine calc_bulk_bcs
    use modglobal,   only : i1, j1, cp, rlv, fkar, zf, cu, cv, grav, rv, rd, lopenbc,lboundary,lperiodic
    use modfields,   only : rhof, thl0, qt0, u0, v0, thvh
    use modsurface,  only : phim, phih
    use modmpi,      only : excjs
    use modopenboundary, only : openboundary_excjs
    use modsurfdata, only : &
        H, LE, G0, tskin, qskin, thlflux, qtflux, dthldz, dqtdz, &
        dudz, dvdz, ustar, obl, cliq, ra, rsveg, rssoil
    implicit none

    integer :: i, j
    real :: rhocp_i, rholv_i, ucu, vcv, cveg, bflux

    rhocp_i = 1. / (rhof(1) * cp)
    rholv_i = 1. / (rhof(1) * rlv)

    ! Calculate surface temperature for each tile, and calculate
    ! surface fluxes (H, LE, G0, wthl, wqt) and values (thlskin, qtskin)
    call calc_tile_bcs(tile_lv)
    call calc_tile_bcs(tile_hv)
    call calc_tile_bcs(tile_bs)
    call calc_tile_bcs(tile_ws)
    call calc_water_bcs

    do j=2,j1
        do i=2,i1
            ! Calc grid point averaged quantities
            H(i,j) = tile_lv%frac(i,j) * tile_lv%H(i,j) + &
                     tile_hv%frac(i,j) * tile_hv%H(i,j) + &
                     tile_bs%frac(i,j) * tile_bs%H(i,j) + &
                     tile_ws%frac(i,j) * tile_ws%H(i,j) + &
                     tile_aq%frac(i,j) * tile_aq%H(i,j)

            LE(i,j) = tile_lv%frac(i,j) * tile_lv%LE(i,j) + &
                      tile_hv%frac(i,j) * tile_hv%LE(i,j) + &
                      tile_bs%frac(i,j) * tile_bs%LE(i,j) + &
                      tile_ws%frac(i,j) * tile_ws%LE(i,j) + &
                      tile_aq%frac(i,j) * tile_aq%LE(i,j)

            G0(i,j) = tile_lv%frac(i,j) * tile_lv%G(i,j) + &
                      tile_hv%frac(i,j) * tile_hv%G(i,j) + &
                      tile_bs%frac(i,j) * tile_bs%G(i,j) + &
                      tile_ws%frac(i,j) * tile_ws%G(i,j) + &
                      tile_aq%frac(i,j) * tile_aq%G(i,j)

            ustar(i,j) = tile_lv%frac(i,j) * tile_lv%ustar(i,j) + &
                         tile_hv%frac(i,j) * tile_hv%ustar(i,j) + &
                         tile_bs%frac(i,j) * tile_bs%ustar(i,j) + &
                         tile_ws%frac(i,j) * tile_ws%ustar(i,j) + &
                         tile_aq%frac(i,j) * tile_aq%ustar(i,j)

            tskin(i,j) = tile_lv%frac(i,j) * tile_lv%thlskin(i,j) + &
                         tile_hv%frac(i,j) * tile_hv%thlskin(i,j) + &
                         tile_bs%frac(i,j) * tile_bs%thlskin(i,j) + &
                         tile_ws%frac(i,j) * tile_ws%thlskin(i,j) + &
                         tile_aq%frac(i,j) * tile_aq%thlskin(i,j)

            qskin(i,j) = tile_lv%frac(i,j) * tile_lv%qtskin(i,j) + &
                         tile_hv%frac(i,j) * tile_hv%qtskin(i,j) + &
                         tile_bs%frac(i,j) * tile_bs%qtskin(i,j) + &
                         tile_ws%frac(i,j) * tile_ws%qtskin(i,j) + &
                         tile_aq%frac(i,j) * tile_aq%qtskin(i,j)

            ! Kinematic surface fluxes
            thlflux(i,j) =  H(i,j)  * rhocp_i
            qtflux (i,j) =  LE(i,j) * rholv_i

            ! Calculate mean Obukhov length from mean fluxes
            bflux = grav/thvh(1) * (thlflux(i,j) * (1.-(1.-rv/rd)*qskin(i,j)) - &
                (1.-rv/rd)*tskin(i,j)*qtflux(i,j))
            obl(i,j) = -ustar(i,j)**3 / (fkar * bflux)

            ! MO gradients
            dthldz(i,j) = -thlflux(i,j) / (fkar * zf(1) * ustar(i,j)) * phih(zf(1)/obl(i,j))
            dqtdz (i,j) = -qtflux (i,j) / (fkar * zf(1) * ustar(i,j)) * phih(zf(1)/obl(i,j))

            ! NOTE: dudz, dvdz are at the grid center (full level), not the velocity locations.
            ucu = 0.5*(u0(i,j,1) + u0(i+1,j,1))+cu
            vcv = 0.5*(v0(i,j,1) + v0(i,j+1,1))+cv

            dudz(i,j) = ustar(i,j) / (fkar * zf(1)) * phim(zf(1)/obl(i,j)) * (ucu/du_tot(i,j))
            dvdz(i,j) = ustar(i,j) / (fkar * zf(1)) * phim(zf(1)/obl(i,j)) * (vcv/du_tot(i,j))

            ! Cyclic BCs where needed.
            if(lopenbc) then ! Only use periodicity for non-domain boundaries when openboundaries are used
              call openboundary_excjs(ustar, 2,i1,2,j1,1,1,1,1, &
                & (.not.lboundary(1:4)).or.lperiodic(1:4))
            else
              call excjs(ustar, 2,i1,2,j1,1,1,1,1)
            endif

            ! Just for diagnostics (modlsmcrosssection)
            cliq(i,j) = tile_ws%frac(i,j) / land_frac(i,j)

            ! Calculate ra consistent with mean flux and temperature difference:
            ra(i,j) = (tskin(i,j) - thl0(i,j,1)) / thlflux(i,j)

            cveg = tile_lv%frac(i,j) + tile_hv%frac(i,j)
            if (cveg > 0) then
                rsveg(i,j) = (tile_lv%frac(i,j) * tile_lv%rs(i,j) + &
                              tile_hv%frac(i,j) * tile_hv%rs(i,j)) / cveg
            else
                rsveg(i,j) = 0.
            end if

            if (tile_bs%frac(i,j) > 0) then
                rssoil(i,j) = tile_bs%rs(i,j)
            else
                rssoil(i,j) = 0.
            end if

        end do
    end do

end subroutine calc_bulk_bcs

!
! Interpolation soil from full to half levels,
! using various interpolation methods
!
subroutine interpolate_soil(fieldh, field, iinterp)
    use modglobal, only : i1, j1
    implicit none
    real, intent(inout) :: fieldh(:,:,:)
    real, intent(in)    :: field(:,:,:)
    integer, intent(in) :: iinterp
    integer i, j, k

    if (iinterp == iinterp_amean) then
        do k=2,kmax_soil
            do j=2,j1
                do i=2,i1
                    fieldh(i,j,k) = 0.5*(field(i,j,k-1) + field(i,j,k))
                end do
            end do
        end do
    else if (iinterp == iinterp_gmean) then
        do k=2,kmax_soil
            do j=2,j1
                do i=2,i1
                    fieldh(i,j,k) = sqrt(field(i,j,k-1) * field(i,j,k))
                end do
            end do
        end do
    else if (iinterp == iinterp_hmean) then
        do k=2,kmax_soil
            do j=2,j1
                do i=2,i1
                    fieldh(i,j,k) = ((dz_soil(k-1)+dz_soil(k))*field(i,j,k-1)*field(i,j,k)) / &
                            (dz_soil(k-1)*field(i,j,k) + dz_soil(k)*field(i,j,k-1))
                end do
            end do
        end do
    else if (iinterp == iinterp_max) then
        do k=2,kmax_soil
            do j=2,j1
                do i=2,i1
                    fieldh(i,j,k) = max(field(i,j,k-1), field(i,j,k))
                end do
            end do
        end do
    else
        print*,'ERROR: invalid soil interpolation type iinterp=',iinterp
        stop
    end if

end subroutine interpolate_soil

!
! Calculate temperature diffusivity soil at full and half levels
!
subroutine calc_thermal_properties
    use modglobal, only : i1, j1, gamma_t_matrix, gamma_t_water
    use modsurfdata, only : lambda, lambdah, phiw
    implicit none

    integer :: i, j, k, si
    real :: lambda_t_sat, kersten, gamma_t, part1, part2, part3

    ! Calculate diffusivity heat
    do k=1,kmax_soil
        do j=2,j1
            do i=2,i1
                si = soil_index(i,j,k)
				
                ! Heat conductivity at saturation (from IFS code..)    
        
    			lambda_T_sat =   gamma_T_matrix**(1.-theta_sat(si)) &
                               * gamma_T_water**phiw(i,j,k) &
                               * 2.2**(theta_sat(si) - phiw(i,j,k))
                               
                ! Kersten number for fine soils [IFS eq 8.64] (-)
                kersten = log10(max(0.1, phiw(i,j,k) / theta_sat(si))) + 1.

                ! Heat conductivity soil [IFS eq 8.62] (W m-1 K-1)
                gamma_t = kersten * (lambda_T_sat - gamma_t_dry(si)) + gamma_t_dry(si)

                ! Heat diffusivity (m2 s-1)
                lambda(i,j,k) = gamma_t / rho_C(si)
            end do
        end do
    end do

    ! Interpolate to half levels
    call interpolate_soil(lambdah, lambda, iinterp_t)

end subroutine calc_thermal_properties

!
! Calculate soil moisture diffusivity and conductivity,
! using van Genuchten parameterisation.
!
subroutine calc_hydraulic_properties
    use modglobal, only : i1, j1
    use modsurfdata, only : phiw, lambdas, lambdash, gammas, gammash
    implicit none

    integer :: i, j, k, si
    real :: theta_lim, theta_norm

    ! Calculate diffusivity and conductivity soil moisture
    do k=1,kmax_soil
        do j=2,j1
            do i=2,i1
                si = soil_index(i,j,k)

                ! Limit soil moisture just above the residual soil moisture content
                !theta_lim = max(phiw(i,j,k), 1.001*theta_res(si))
                theta_lim = min(max(phiw(i,j,k), 1.001*theta_res(si)), 0.999*theta_sat(si))

                ! Dimensionless soil water content
                theta_norm = (theta_lim - theta_res(si)) / (theta_sat(si) - theta_res(si))
                
                if((isnan(theta_norm)).OR.(theta_norm<=0.)) theta_norm=0.1          !theta_norm could be 0, so to avoid a complex number I put this! Arseni
            

                ! Calculate & limit the diffusivity
                lambdas(i,j,k) = calc_diffusivity_vg( &
                        theta_norm, vg_a(si), vg_l(si), vg_m(si), &
                        gamma_theta_sat(si), theta_sat(si), theta_res(si))
                
                lambdas(i,j,k) = max(min(lambdas(i,j,k), lambda_theta_max(si)), lambda_theta_min(si))

                ! Calculate & limit the conductivity
                gammas(i,j,k) = calc_conductivity_vg( &
                        theta_norm, vg_l(si), vg_m(si), gamma_theta_sat(si))
                gammas(i,j,k) = max(min(gammas(i,j,k), gamma_theta_max(si)), gamma_theta_min(si))
            end do
        end do
    end do

    ! Interpolate to half levels
    call interpolate_soil(lambdash, lambdas, iinterp_theta)
    call interpolate_soil(gammash,  gammas,  iinterp_theta)

    ! Optionally, set free drainage bottom BC
    if (lfreedrainage) then
        gammash(:,:,1) = gammash(:,:,2)
    else
        gammash(:,:,1) = 0.
    end if

end subroutine calc_hydraulic_properties

!
! Calculate soil moisture tendency from root water extraction
!
subroutine calc_root_water_extraction
    use modglobal, only : i1, j1, rhow, rlv
    use modsurfdata, only : phiw
    implicit none

    integer :: i, j, k
    real :: phiw_rf_lv, phiw_rf_hv, phi_frac_lv, phi_frac_hv, LE_lv, LE_hv
    real, parameter :: fac = 1./(rhow * rlv)

    do j=2, j1
        do i=2, i1

            LE_lv = tile_lv%frac(i,j) * tile_lv%LE(i,j)
            LE_hv = tile_hv%frac(i,j) * tile_hv%LE(i,j)

            phiw_rf_lv = 0.
            phiw_rf_hv = 0.
            do k=1, kmax_soil
                phiw_rf_lv = phiw_rf_lv + tile_lv%root_frac(i,j,k) * phiw(i,j,k)
                phiw_rf_hv = phiw_rf_hv + tile_hv%root_frac(i,j,k) * phiw(i,j,k)
            end do

            do k=1, kmax_soil
                phi_frac_lv = tile_lv%root_frac(i,j,k) * phiw(i,j,k) / phiw_rf_lv
                phi_frac_hv = tile_hv%root_frac(i,j,k) * phiw(i,j,k) / phiw_rf_hv

                phiw_source(i,j,k) = -max(0., LE_lv) * fac * dzi_soil(k) * phi_frac_lv &
                                     -max(0., LE_hv) * fac * dzi_soil(k) * phi_frac_hv
            end do
        end do
    end do

end subroutine calc_root_water_extraction

!
! Solve diffusion equation (explicit) for soil temperature.
! Top flux = G0/rho_C, bottom flux = zero.
!
subroutine integrate_t_soil
    use modglobal, only : rk3step, rdt, i1, j1
    use modsurfdata, only : tsoil, tsoilm, lambdah, G0

    implicit none
    integer :: i, j, k, si
    real :: tend, rk3coef, flux_top

    rk3coef = rdt / (4. - dble(rk3step))
    if(rk3step == 1) tsoilm(:,:,:) = tsoil(:,:,:)

    ! Top soil layer
    k = kmax_soil
    do j=2,j1
        do i=2,i1
            si = soil_index(i,j,k)
            flux_top = G0(i,j) / rho_C(si)
            tend = (-flux_top - (lambdah(i,j,k) * (tsoil(i,j,k) - tsoil(i,j,k-1)) * dzhi_soil(k)))*dzi_soil(k)

            tsoil(i,j,k) = tsoilm(i,j,k) + rk3coef * tend
        end do
    end do

    ! Bottom soil layer
    k = 1
    do j=2,j1
        do i=2,i1
            tend = ((lambdah(i,j,k+1) * (tsoil(i,j,k+1) - tsoil(i,j,k)) * dzhi_soil(k+1)))*dzi_soil(k)

            tsoil(i,j,k) = tsoilm(i,j,k) + rk3coef * tend
        end do
    end do

    ! Interior
    do k=2,kmax_soil-1
        do j=2,j1
            do i=2,i1
                tend = ((lambdah(i,j,k+1) * (tsoil(i,j,k+1) - tsoil(i,j,k  )) * dzhi_soil(k+1)) &
                      - (lambdah(i,j,k  ) * (tsoil(i,j,k  ) - tsoil(i,j,k-1)) * dzhi_soil(k  ))) * dzi_soil(k)

                tsoil(i,j,k) = tsoilm(i,j,k) + rk3coef * tend
            end do
        end do
    end do

end subroutine integrate_t_soil

!
! Solve diffusion equation (explicit) for soil moisture,
! including a source term for root water extraction.
!
subroutine integrate_theta_soil
    use modglobal, only : rk3step, rdt, i1, j1, rhow, rlv
    use modsurfdata, only : phiw, phiwm, lambdash, gammash

    implicit none
    integer :: i, j, k, si
    real :: tend, rk3coef, flux_top, fac

    rk3coef = rdt / (4. - dble(rk3step))
    if(rk3step == 1) phiwm(:,:,:) = phiw(:,:,:)

    fac = 1./(rhow * rlv)

    ! Top soil layer
    k = kmax_soil
    do j=2,j1
        do i=2,i1
            flux_top = tile_bs%frac(i,j) * tile_bs%LE(i,j) * fac + throughfall(i,j)
            tend = (-flux_top - (lambdash(i,j,k) * (phiw(i,j,k) - phiw(i,j,k-1)) * dzhi_soil(k)))*dzi_soil(k) &
                   - gammash(i,j,k) * dzi_soil(k) + phiw_source(i,j,k)

            phiw(i,j,k) = phiwm(i,j,k) + rk3coef * tend
        end do
    end do

    ! Bottom soil layer
    k = 1
    do j=2,j1
        do i=2,i1
            tend = ((lambdash(i,j,k+1) * (phiw(i,j,k+1) - phiw(i,j,k)) * dzhi_soil(k+1)))*dzi_soil(k) &
                   + (gammash(i,j,k+1) - gammash(i,j,k)) * dzi_soil(k) + phiw_source(i,j,k)

            phiw(i,j,k) = phiwm(i,j,k) + rk3coef * tend
        end do
    end do

    ! Interior
    do k=2,kmax_soil-1
        do j=2,j1
            do i=2,i1
                tend = ((lambdash(i,j,k+1) * (phiw(i,j,k+1) - phiw(i,j,k  )) * dzhi_soil(k+1)) &
                      - (lambdash(i,j,k  ) * (phiw(i,j,k  ) - phiw(i,j,k-1)) * dzhi_soil(k  ))) * dzi_soil(k) &
                      + (gammash(i,j,k+1) - gammash(i,j,k)) * dzi_soil(k) + phiw_source(i,j,k)

                phiw(i,j,k) = phiwm(i,j,k) + rk3coef * tend
            end do
        end do
    end do

end subroutine integrate_theta_soil

!
! Initialise the land-surface model
!
subroutine initlsm
    use modglobal,   only : ifnamopt, fname_options, checknamelisterror, lwarmstart
    use modmpi,      only : myid, comm3d, mpierr, mpi_logical, mpi_integer, my_real
    use modsurfdata, only : isurf
    use modemisdata, only : l_emission
    implicit none

    integer :: ierr
    logical :: lheterogeneous

    ! Namelist definition
    namelist /NAMLSM/ &
        lheterogeneous, lfreedrainage, lags, dz_soil, iinterp_t, iinterp_theta, co2_index

    llsm = (isurf == 11)

    if (llsm) then

        allocate(dz_soil(kmax_soil))

        ! Read namelist
        if (myid == 0) then
            open(ifnamopt, file=fname_options, status='old', iostat=ierr)
            read(ifnamopt, NAMLSM, iostat=ierr)
            call checknamelisterror(ierr, ifnamopt, 'NAMLSM')
            write(6, NAMLSM)
            close(ifnamopt)
        end if

        ! Broadcast namelist values to all MPI tasks
        call MPI_BCAST(lheterogeneous, 1, mpi_logical, 0, comm3d, mpierr)
        call MPI_BCAST(lfreedrainage,  1, mpi_logical, 0, comm3d, mpierr)
        call MPI_BCAST(lags,           1, mpi_logical, 0, comm3d, mpierr)

        call MPI_BCAST(iinterp_t,      1, mpi_integer, 0, comm3d, mpierr)
        call MPI_BCAST(iinterp_theta,  1, mpi_integer, 0, comm3d, mpierr)

        call MPI_BCAST(co2_index,      1, mpi_integer, 0, comm3d, mpierr)

        call MPI_BCAST(dz_soil, kmax_soil, my_real, 0, comm3d, mpierr)

        ! Checks on input
        if (lags .and. .not. l_emission .and. (co2_index == -1)) then
            print*,'With A-Gs enabled and without the emission module, the `co2_index`'
            print*,'in the scalar array should be specified in the `NAMLSM` group.'
            stop
        endif

        ! Create/calculate soil grid properties
        call create_soil_grid

        ! Allocate required fields in modsurfacedata, and arrays / tiles from this module:
        call allocate_fields

        if (lheterogeneous) then
            ! Initialise heterogeneous LSM from external input
            call init_heterogeneous
        else
            ! Initialise homogeneous LSM from namelist input
            call init_homogeneous
        end if

        ! Read the soil parameter table
        call read_soil_table

        ! Calculate derived soil properties
        call calc_soil_properties

        ! Calculate root fractions soil (low and high veg)
        call calc_root_fractions
    end if

end subroutine initlsm

!
! Cleanup (deallocate) the land-surface model
!
subroutine exitlsm
    use modsurfdata, only : &
        tsoil, tsoilm, phiw, phiwm, &
        lambda, lambdah, lambdas, lambdash, gammas, gammash, &
        H, LE, G0, Qnet, wl, wlm, &
        tendskin, cliq, rsveg, rssoil
    implicit none

    if (.not. llsm) return

    ! Allocated from `read_soil_table`:
    deallocate( theta_res, theta_wp, theta_fc, theta_sat, gamma_theta_sat, vg_a, vg_l, vg_n )

    ! Allocated from `allocate_fields`:
    deallocate( soil_index, tsoil, tsoilm, phiw, phiwm, phiw_source, index_ags )
    deallocate( wl, wlm, wl_max )
    deallocate( throughfall, interception )
    deallocate( gD, f1, f2_lv, f2_hv, f2b, f3 )
    deallocate( du_tot, thv_1, land_frac )
    deallocate( lambda, lambdah, lambdas, lambdash, gammas, gammash )
    deallocate( Qnet, H, LE, G0 )
    deallocate( tendskin, cliq, rsveg, rssoil )

    ! Allocated from `create_soil_grid`:
    deallocate( z_soil, dz_soil, dzi_soil, zh_soil, dzh_soil, dzhi_soil )

    if (lags) deallocate(an_co2, resp_co2)

    ! Tiles, allocated from `allocate_tile`:
    call deallocate_tile(tile_lv)
    call deallocate_tile(tile_hv)
    call deallocate_tile(tile_bs)
    call deallocate_tile(tile_ws)
    call deallocate_tile(tile_aq)

end subroutine exitlsm

!
! Calculate soil grid properties
!
subroutine create_soil_grid
    implicit none
    integer :: k

    allocate(z_soil   (kmax_soil  ))
    allocate(dzi_soil (kmax_soil  ))
    allocate(zh_soil  (kmax_soil+1))
    allocate(dzh_soil (kmax_soil+1))
    allocate(dzhi_soil(kmax_soil+1))

    !
    ! Full level is in the middle of two half levels
    !
    z_size_soil = -sum(dz_soil)

    ! Half level height
    zh_soil(1) = z_size_soil
    zh_soil(kmax_soil+1) = 0.
    do k=2, kmax_soil
        zh_soil(k) = zh_soil(k-1) + dz_soil(k-1)
    end do

    ! Full level height
    do k=1, kmax_soil
        z_soil(k) = 0.5*(zh_soil(k) + zh_soil(k+1))
    end do

    do k=2, kmax_soil
        dzh_soil(k) = z_soil(k) - z_soil(k-1)
    end do

    dzh_soil(1) = 2*(z_soil(1) - zh_soil(1))
    dzh_soil(kmax_soil+1) = 2*(-z_soil(kmax_soil))

    ! Inverse grid spacings
    dzi_soil(:) = 1./dz_soil(:)
    dzhi_soil(:) = 1./dzh_soil(:)

    !
    ! Half level is in the middle of two full levels
    !
    !! Half level heights
    !zh_soil(1) = z_size_soil
    !zh_soil(kmax_soil+1) = 0.

    !do k=2, kmax_soil
    !    zh_soil(k) = 0.5*(z_soil(k) + z_soil(k-1))
    !end do

    !! Grid spacing full and half levels
    !do k=1, kmax_soil
    !    dz_soil(k) = zh_soil(k+1) - zh_soil(k)
    !end do

    !do k=2, kmax_soil
    !    dzh_soil(k) = z_soil(k) - z_soil(k-1)
    !end do

    !dzh_soil(1) = 2*(z_soil(1) - zh_soil(1))
    !dzh_soil(kmax_soil+1) = 2*(-z_soil(kmax_soil))

    !! Inverse grid spacings
    !dzi_soil(:) = 1./dz_soil(:)
    !dzhi_soil(:) = 1./dzh_soil(:)

end subroutine create_soil_grid

!
! Allocate all LSM fields
!
subroutine allocate_fields
    use modglobal, only : i2, j2
    use modsurfdata, only : &
        tsoil, tsoilm, phiw, phiwm, &
        lambda, lambdah, lambdas, lambdash, gammas, gammash, &
        H, LE, G0, Qnet, wl, wlm, &
        tendskin, cliq, rsveg, rssoil
    implicit none

    ! Allocate soil variables
    allocate(soil_index(i2, j2, kmax_soil))

    allocate(tsoil (i2, j2, kmax_soil))
    allocate(tsoilm(i2, j2, kmax_soil))
    allocate(phiw  (i2, j2, kmax_soil))
    allocate(phiwm (i2, j2, kmax_soil))

    allocate(wl (i2, j2))
    allocate(wlm(i2, j2))

    allocate(wl_max(i2, j2))

    allocate(phiw_source(i2, j2, kmax_soil))

    allocate(throughfall(i2, j2))
    allocate(interception(i2, j2))

    allocate(gD(i2, j2))
    allocate(f1(i2, j2))
    allocate(f2_lv(i2, j2))
    allocate(f2_hv(i2, j2))
    allocate(f2b(i2, j2))
    allocate(f3(i2, j2))

    allocate(du_tot(i2, j2))
    allocate(thv_1(i2, j2))
    allocate(land_frac(i2, j2))

    ! NOTE: names differ from what is described in modsurfdata!
    ! Diffusivity temperature:
    allocate(lambda (i2, j2, kmax_soil  ))
    allocate(lambdah(i2, j2, kmax_soil+1))

    ! Diffusivity theta:
    allocate(lambdas (i2, j2, kmax_soil))
    allocate(lambdash(i2, j2, kmax_soil+1))

    ! Conductivity theta:
    allocate(gammas (i2, j2, kmax_soil))
    allocate(gammash(i2, j2, kmax_soil+1))

    ! Tile averaged surface fluxes
    allocate(Qnet (i2, j2))
    allocate(H    (i2, j2))
    allocate(LE   (i2, j2))
    allocate(G0   (i2, j2))

    ! TMP, to prevent segfaults in modlsmcrosssection.f90
    allocate(tendskin(i2, j2))
    allocate(cliq(i2, j2))
    allocate(rsveg(i2, j2))
    allocate(rssoil(i2, j2))

    ! A-Gs
    if (lags) then
        allocate(an_co2(i2, j2))
        allocate(resp_co2(i2, j2))
    endif
    
    ! A-Gs vegetation type/index
    allocate(index_ags(i2, j2))

    ! Allocate the tiled variables
    call allocate_tile(tile_lv)
    call allocate_tile(tile_hv)
    call allocate_tile(tile_bs)
    call allocate_tile(tile_ws)
    call allocate_tile(tile_aq)

end subroutine allocate_fields

!
! Allocate all fields of a LSM tile
!
subroutine allocate_tile(tile)
    use modglobal, only : i2, j2
    implicit none
    type(lsm_tile), intent(inout) :: tile

    ! Static properties:
    allocate(tile % z0m(i2, j2))
    allocate(tile % z0h(i2, j2))

    ! Base and dynamic tile fraction:
    allocate(tile % base_frac(i2, j2))
    allocate(tile % frac(i2, j2))

    ! Monin-obukhov / surface layer:
    allocate(tile % obuk(i2, j2))
    allocate(tile % ustar(i2, j2))
    allocate(tile % ra(i2, j2))

    ! Conductivity skin layer:
    allocate(tile % lambda_stable(i2, j2))
    allocate(tile % lambda_unstable(i2, j2))

    ! Surface fluxes:
    allocate(tile % H(i2, j2))
    allocate(tile % LE(i2, j2))
    allocate(tile % G(i2, j2))
    allocate(tile % wthl(i2, j2))
    allocate(tile % wqt(i2, j2))

    ! Surface temperature and humidity:
    allocate(tile % tskin(i2, j2))
    allocate(tile % thlskin(i2, j2))
    allocate(tile % qtskin(i2, j2))

    ! Buoyancy difference surface-atmosphere
    allocate(tile % db(i2, j2))

    ! Vegetation properties:
    allocate(tile % lai(i2, j2))
    allocate(tile % rs_min(i2, j2))
    allocate(tile % rs(i2, j2))

    ! Root fraction in soil:
    allocate(tile % a_r(i2, j2))
    allocate(tile % b_r(i2, j2))
    allocate(tile % root_frac(i2, j2, kmax_soil))

    ! Root fraction weigthed mean soil water content
    allocate(tile % phiw_mean(i2, j2))

end subroutine allocate_tile

!
! Deallocate all fields of a LSM tile
!
subroutine deallocate_tile(tile)
    implicit none
    type(lsm_tile), intent(inout) :: tile

    deallocate( tile%z0m, tile%z0h, tile%base_frac, tile%frac )
    deallocate( tile%obuk, tile%ustar, tile%ra )
    deallocate( tile%lambda_stable, tile%lambda_unstable )
    deallocate( tile%H, tile%LE, tile%G, tile%wthl, tile%wqt)
    deallocate( tile%tskin, tile%thlskin, tile%qtskin)
    deallocate( tile%db, tile%lai, tile%rs_min, tile%rs )
    deallocate( tile%a_r, tile%b_r, tile%root_frac, tile%phiw_mean )

end subroutine deallocate_tile

!
! Init some of the tiled variables, in case of cold start.
! Called from modstartup -> readinitfiles()
!
subroutine init_lsm_tiles
    use modfields, only : thlprof, qtprof
    implicit none

    tile_lv % thlskin(:,:) = thlprof(1)
    tile_lv % qtskin (:,:) = qtprof(1)
    tile_lv % obuk   (:,:) = -0.1

    tile_hv % thlskin(:,:) = thlprof(1)
    tile_hv % qtskin (:,:) = qtprof(1)
    tile_hv % obuk   (:,:) = -0.1

    tile_bs % thlskin(:,:) = thlprof(1)
    tile_bs % qtskin (:,:) = qtprof(1)
    tile_bs % obuk   (:,:) = -0.1

    tile_ws % thlskin(:,:) = thlprof(1)
    tile_ws % qtskin (:,:) = qtprof(1)
    tile_ws % obuk   (:,:) = -0.1

    tile_aq % thlskin(:,:) = thlprof(1)
    tile_aq % qtskin (:,:) = qtprof(1)
    tile_aq % obuk   (:,:) = -0.1

end subroutine init_lsm_tiles

!
! Initialise the LSM homogeneous from namelist input
!
subroutine init_homogeneous
    use modglobal,   only : ifnamopt, fname_options, checknamelisterror, lwarmstart, eps1
    use modmpi,      only : myid, comm3d, mpierr, mpi_logical, my_real, mpi_integer
    use modsurfdata, only : tsoil, tsoilm, phiw, phiwm, wl, wlm, wmax
    implicit none

    integer :: ierr, k
    real :: c_low, c_high, c_bare, c_water
    real :: z0m_low, z0m_high, z0m_bare, z0m_water
    real :: z0h_low, z0h_high, z0h_bare, z0h_water
    real :: lambda_s_low, lambda_s_high, lambda_s_bare
    real :: lambda_us_low, lambda_us_high, lambda_us_bare
    real :: lai_low, lai_high
    real :: rs_min_low, rs_min_high, rs_min_bare
    real :: ar_low, br_low, ar_high, br_high
    real :: gD_high
    real :: tskin_water

    real, allocatable :: t_soil_p(:), theta_soil_p(:)
    integer, allocatable :: soil_index_p(:)

    ! Read namelist
    namelist /NAMLSM_HOMOGENEOUS/ &
        c_low, c_high, c_bare, c_water, &
        z0m_low, z0m_high, z0m_bare, z0m_water, &
        z0h_low, z0h_high, z0h_bare, z0h_water, &
        lambda_s_low, lambda_s_high, lambda_s_bare, &
        lambda_us_low, lambda_us_high, lambda_us_bare, &
        lai_low, lai_high, &
        rs_min_low, rs_min_high, rs_min_bare, &
        t_soil_p, theta_soil_p, soil_index_p, &
        ar_low, br_low, ar_high, br_high, &
        gD_high, tskin_water

    allocate(t_soil_p(kmax_soil), theta_soil_p(kmax_soil))
    allocate(soil_index_p(kmax_soil))

    if (myid == 0) then
        open(ifnamopt, file=fname_options, status='old', iostat=ierr)
        read(ifnamopt, NAMLSM_HOMOGENEOUS, iostat=ierr)
        call checknamelisterror(ierr, ifnamopt, 'NAMLSM_HOMOGENEOUS')
        write(6, NAMLSM_HOMOGENEOUS)
        close(ifnamopt)
    end if

    ! Broadcast to all MPI tasks
    call MPI_BCAST(c_low,   1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(c_high,  1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(c_bare,  1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(c_water, 1, my_real, 0, comm3d, mpierr)

    call MPI_BCAST(z0m_low,   1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(z0m_high,  1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(z0m_bare,  1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(z0m_water, 1, my_real, 0, comm3d, mpierr)

    call MPI_BCAST(z0h_low,   1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(z0h_high,  1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(z0h_bare,  1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(z0h_water, 1, my_real, 0, comm3d, mpierr)

    call MPI_BCAST(lambda_s_low,  1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(lambda_s_high, 1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(lambda_s_bare, 1, my_real, 0, comm3d, mpierr)

    call MPI_BCAST(lambda_us_low,  1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(lambda_us_high, 1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(lambda_us_bare, 1, my_real, 0, comm3d, mpierr)

    call MPI_BCAST(lai_low,  1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(lai_high, 1, my_real, 0, comm3d, mpierr)

    call MPI_BCAST(rs_min_low,  1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(rs_min_high, 1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(rs_min_bare, 1, my_real, 0, comm3d, mpierr)

    call MPI_BCAST(t_soil_p,     kmax_soil, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(theta_soil_p, kmax_soil, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(soil_index_p, kmax_soil, mpi_integer, 0, comm3d, mpierr)

    call MPI_BCAST(ar_low,  1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(br_low,  1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(ar_high, 1, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(br_high, 1, my_real, 0, comm3d, mpierr)

    call MPI_BCAST(gD_high, 1, my_real, 0, comm3d, mpierr)

    call MPI_BCAST(tskin_water, 1, my_real, 0, comm3d, mpierr)

    ! Set values
    tile_lv % base_frac(:,:) = c_low
    tile_hv % base_frac(:,:) = c_high
    tile_bs % base_frac(:,:) = c_bare
    tile_aq % base_frac(:,:) = c_water
    tile_ws % base_frac(:,:) = 0.

    tile_lv % z0m(:,:) = z0m_low
    tile_hv % z0m(:,:) = z0m_high
    tile_bs % z0m(:,:) = z0m_bare
    tile_aq % z0m(:,:) = z0m_water

    tile_lv % z0h(:,:) = z0h_low
    tile_hv % z0h(:,:) = z0h_high
    tile_bs % z0h(:,:) = z0h_bare
    tile_aq % z0h(:,:) = z0h_water

    tile_lv % lambda_stable(:,:) = lambda_s_low
    tile_hv % lambda_stable(:,:) = lambda_s_high
    tile_bs % lambda_stable(:,:) = lambda_s_bare

    tile_lv % lambda_unstable(:,:) = lambda_us_low
    tile_hv % lambda_unstable(:,:) = lambda_us_high
    tile_bs % lambda_unstable(:,:) = lambda_us_bare

    tile_lv % lai(:,:) = lai_low
    tile_hv % lai(:,:) = lai_high

    tile_lv % rs_min(:,:) = rs_min_low
    tile_hv % rs_min(:,:) = rs_min_high
    tile_bs % rs_min(:,:) = rs_min_bare

    tile_lv % a_r(:,:) = ar_low
    tile_lv % b_r(:,:) = br_low
    tile_hv % a_r(:,:) = ar_high
    tile_hv % b_r(:,:) = br_high

    gD(:,:) = gD_high

    tile_aq % tskin(:,:) = tskin_water


    if (.not. lwarmstart) then
        ! Init prognostic variables in case of cold start.
        ! For a warm start, these are read from restartfiles
        ! in modstartup.
        do k=1, kmax_soil
            tsoil(:,:,k) = t_soil_p(k)
            phiw (:,:,k) = theta_soil_p(k)
        end do

        tsoilm(:,:,:) = tsoil(:,:,:)
        phiwm (:,:,:) = phiw (:,:,:)

        wl(:,:)  = 0.
        wlm(:,:) = 0.
    end if

    do k=1, kmax_soil
        soil_index(:,:,k) = soil_index_p(k)
    end do

    ! Set properties wet skin tile
    tile_ws % z0m(:,:) = c_low*z0m_low + c_high*z0m_high + c_bare*z0m_bare
    tile_ws % z0h(:,:) = c_low*z0h_low + c_high*z0h_high + c_bare*z0h_bare

    tile_ws % lambda_stable(:,:)   = c_low*lambda_s_low  + c_high*lambda_s_high  + c_bare*lambda_s_bare
    tile_ws % lambda_unstable(:,:) = c_low*lambda_us_low + c_high*lambda_us_high + c_bare*lambda_us_bare

    ! Calculate land fraction, and limit to prevent div/0's
    land_frac(:,:) = 1.-c_water
    where (land_frac == 0) land_frac = eps1

    ! Max liquid water per grid point, accounting for LAI
    wl_max(:,:) = wmax * ( &
            tile_lv%base_frac(:,:) * tile_lv%lai(:,:) + &
            tile_hv%base_frac(:,:) * tile_hv%lai(:,:) + &
            tile_bs%base_frac(:,:)) / land_frac(:,:)
    where (wl_max == 0) wl_max = eps1

    ! Cleanup!
    deallocate(t_soil_p, theta_soil_p, soil_index_p)

end subroutine init_homogeneous

!
! Init LSM properties from external input
!
subroutine init_heterogeneous
    use modglobal,   only : i1, j1, lwarmstart, iexpnr, eps1
    use modmpi,      only : myidx, myidy
    use modsurfdata, only : tsoil, phiw, wl, wlm, wmax
    implicit none

    character(20) :: input_file = 'lsm.inp.x___y___.___'
    write(input_file(10:12), '(i3.3)') myidx
    write(input_file(14:16), '(i3.3)') myidy
    write(input_file(18:20), '(i3.3)') iexpnr
    
    ! Open the unformatted files and read the heterogeneous lsm input
	print*,'Reading heterogeneous input for lsm, files: ',input_file
    

    open(666, file=input_file, form='unformatted', status='unknown', action='read', access='stream')

    ! 2D surface fields
    read(666) tile_lv%base_frac(2:i1, 2:j1)
    read(666) tile_hv%base_frac(2:i1, 2:j1)
    read(666) tile_bs%base_frac(2:i1, 2:j1)
    read(666) tile_aq%base_frac(2:i1, 2:j1)

    read(666) tile_lv%z0m(2:i1, 2:j1)
    read(666) tile_hv%z0m(2:i1, 2:j1)
    read(666) tile_bs%z0m(2:i1, 2:j1)
    read(666) tile_aq%z0m(2:i1, 2:j1)

    read(666) tile_lv%z0h(2:i1, 2:j1)
    read(666) tile_hv%z0h(2:i1, 2:j1)
    read(666) tile_bs%z0h(2:i1, 2:j1)
    read(666) tile_aq%z0h(2:i1, 2:j1)

    read(666) tile_lv%lambda_stable(2:i1, 2:j1)
    read(666) tile_hv%lambda_stable(2:i1, 2:j1)
    read(666) tile_bs%lambda_stable(2:i1, 2:j1)

    read(666) tile_lv%lambda_unstable(2:i1, 2:j1)
    read(666) tile_hv%lambda_unstable(2:i1, 2:j1)
    read(666) tile_bs%lambda_unstable(2:i1, 2:j1)

    read(666) tile_lv%lai(2:i1, 2:j1)
    read(666) tile_hv%lai(2:i1, 2:j1)

    read(666) tile_lv%rs_min(2:i1, 2:j1)
    read(666) tile_hv%rs_min(2:i1, 2:j1)
    read(666) tile_bs%rs_min(2:i1, 2:j1)

    read(666) tile_lv%a_r(2:i1, 2:j1)
    read(666) tile_lv%b_r(2:i1, 2:j1)
    read(666) tile_hv%a_r(2:i1, 2:j1)
    read(666) tile_hv%b_r(2:i1, 2:j1)

    read(666) gD(2:i1, 2:j1)

    read(666) tile_aq%tskin(2:i1, 2:j1)

    ! 3D soil fields
    read(666) soil_index(2:i1, 2:j1, 1:kmax_soil)

    if (.not. lwarmstart) then
        read(666) tsoil(2:i1, 2:j1, 1:kmax_soil)
        read(666) phiw (2:i1, 2:j1, 1:kmax_soil)

        wl(:,:)  = 0.
        wlm(:,:) = 0.
    end if
    
    ! A-Gs vegetation type/index
    read(666) index_ags(2:i1, 2:j1)
    

    close(666)

    ! Derived quantities
    tile_ws%base_frac(:,:) = 0.

    ! Set properties wet skin tile
    tile_ws%z0m(:,:) = &
        tile_lv%base_frac(:,:)*tile_lv%z0m(:,:) + &
        tile_hv%base_frac(:,:)*tile_lv%z0m(:,:) + &
        tile_bs%base_frac(:,:)*tile_bs%z0m(:,:)

    tile_ws%z0h(:,:) = &
        tile_lv%base_frac(:,:)*tile_lv%z0h(:,:) + &
        tile_hv%base_frac(:,:)*tile_lv%z0h(:,:) + &
        tile_bs%base_frac(:,:)*tile_bs%z0h(:,:)

    tile_ws%lambda_stable(:,:) = &
        tile_lv%base_frac(:,:)*tile_lv%lambda_stable(:,:) + &
        tile_hv%base_frac(:,:)*tile_lv%lambda_stable(:,:) + &
        tile_bs%base_frac(:,:)*tile_bs%lambda_stable(:,:)

    tile_ws%lambda_unstable(:,:) = &
        tile_lv%base_frac(:,:)*tile_lv%lambda_unstable(:,:) + &
        tile_hv%base_frac(:,:)*tile_lv%lambda_unstable(:,:) + &
        tile_bs%base_frac(:,:)*tile_bs%lambda_unstable(:,:)

    ! Calculate land fraction, and limit to prevent div/0's
    land_frac(:,:) = 1.-tile_aq%base_frac(:,:)
    where (land_frac == 0) land_frac = eps1

    ! Max liquid water per grid point, accounting for LAI
    wl_max(:,:) = wmax * ( &
            tile_lv%base_frac(:,:) * tile_lv%lai(:,:) + &
            tile_hv%base_frac(:,:) * tile_hv%lai(:,:) + &
            tile_bs%base_frac(:,:)) / land_frac(:,:)
    where (wl_max == 0) wl_max = eps1

end subroutine init_heterogeneous

!
! Read the input table with the (van Genuchten) soil parameters
!
subroutine read_soil_table
    use modmpi, only : myid, comm3d, mpierr, mpi_logical, my_real, mpi_integer
    implicit none
    integer :: table_size, ncid, dimid, varid

    if (myid == 0) then
        ! Open the NetCDF file and read the table size
        print*,'Reading "van_genuchten_parameters.nc"'
        call check( nf90_open('van_genuchten_parameters.nc', nf90_nowrite, ncid) )
        call check( nf90_inq_dimid(ncid, 'index', dimid) )
        call check( nf90_inquire_dimension(ncid, dimid, len=table_size) )
    end if

    call MPI_BCAST(table_size, 1, mpi_integer, 0, comm3d, mpierr)

    ! Allocate variables
    allocate( theta_res(table_size), theta_wp(table_size), theta_fc(table_size) )
    allocate( theta_sat(table_size), gamma_theta_sat(table_size) )
    allocate( vg_a(table_size), vg_l(table_size), vg_n(table_size) )

    if (myid == 0) then
        ! Read variables
        call check( nf90_inq_varid(ncid, 'theta_res', varid) )
        call check( nf90_get_var(ncid, varid, theta_res) )

        call check( nf90_inq_varid(ncid, 'theta_wp', varid) )
        call check( nf90_get_var(ncid, varid, theta_wp) )

        call check( nf90_inq_varid(ncid, 'theta_fc', varid) )
        call check( nf90_get_var(ncid, varid, theta_fc) )

        call check( nf90_inq_varid(ncid, 'theta_sat', varid) )
        call check( nf90_get_var(ncid, varid, theta_sat) )

        call check( nf90_inq_varid(ncid, 'gamma_sat', varid) )
        call check( nf90_get_var(ncid, varid, gamma_theta_sat) )

        call check( nf90_inq_varid(ncid, 'alpha', varid) )
        call check( nf90_get_var(ncid, varid, vg_a) )

        call check( nf90_inq_varid(ncid, 'l', varid) )
        call check( nf90_get_var(ncid, varid, vg_l) )

        call check( nf90_inq_varid(ncid, 'n', varid) )
        call check( nf90_get_var(ncid, varid, vg_n) )

        call check( nf90_close(ncid) )
    end if

    ! Broadcast to other MPI tasks
    call MPI_BCAST(theta_res,       table_size, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(theta_wp,        table_size, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(theta_fc,        table_size, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(theta_sat,       table_size, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(gamma_theta_sat, table_size, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(vg_a,            table_size, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(vg_l,            table_size, my_real, 0, comm3d, mpierr)
    call MPI_BCAST(vg_n,            table_size, my_real, 0, comm3d, mpierr)

end subroutine read_soil_table

!
! Calculate derived (tabulated) soil properties
!
subroutine calc_soil_properties
    use modglobal, only : rho_solid_soil, rho_c_matrix, rho_c_water, eps1
    implicit none

    integer :: i, table_size
    real :: theta_norm_min, theta_norm_max, rho_dry

    table_size = size(vg_n)
    allocate(vg_m(table_size))
    allocate(lambda_theta_min(table_size))
    allocate(lambda_theta_max(table_size))
    allocate(gamma_theta_min(table_size))
    allocate(gamma_theta_max(table_size))
    allocate(gamma_t_dry(table_size))
    allocate(rho_C(table_size))

    do i=1, table_size
        ! van Genuchten parameter `m`
        vg_m(i) = (1. - (1. / vg_n(i)))

        ! Min/max values diffusivity soil moisture
        theta_norm_min = (1.001 * theta_res(i) - theta_res(i)) / (theta_sat(i) - theta_res(i)) + eps1
        theta_norm_max = (0.999 * theta_sat(i) - theta_res(i)) / (theta_sat(i) - theta_res(i))

        lambda_theta_min(i) = calc_diffusivity_vg( &
                theta_norm_min, vg_a(i), vg_l(i), vg_m(i), gamma_theta_sat(i), &
                theta_sat(i), theta_res(i))
        lambda_theta_max(i) = calc_diffusivity_vg( &
                theta_norm_max, vg_a(i), vg_l(i), vg_m(i), gamma_theta_sat(i), &
                theta_sat(i), theta_res(i))

        ! Min/max values conductivity soil moisture
        gamma_theta_min(i) = 0.
        gamma_theta_max(i) = gamma_theta_sat(i)

        ! Conductivity temperature
        rho_dry = (1. - theta_sat(i)) * rho_solid_soil  ! Density of soil (kg m-3)
        gamma_t_dry(i) = (0.135 * rho_dry + 64.7) / (rho_solid_soil - 0.947 * rho_dry)
        rho_C(i) = (1. - theta_sat(i)) * rho_C_matrix + theta_fc(i) * rho_C_water
    end do

end subroutine calc_soil_properties

!
! Calculate root fraction for the low and high vegetation tiles
!
subroutine calc_root_fractions
    use modglobal, only : i1, j1
    implicit none
    real :: root_sum_lv, root_sum_hv
    integer i, j, k

    do j=2,j1
        do i=2,i1
            root_sum_lv = 0
            root_sum_hv = 0
            do k=2, kmax_soil

                tile_lv%root_frac(i,j,k) = 0.5 * (&
                    exp(tile_lv%a_r(i,j) * zh_soil(k+1)) + &
                    exp(tile_lv%b_r(i,j) * zh_soil(k+1)) - &
                    exp(tile_lv%a_r(i,j) * zh_soil(k  )) - &
                    exp(tile_lv%b_r(i,j) * zh_soil(k  )))

                tile_hv%root_frac(i,j,k) = 0.5 * (&
                    exp(tile_hv%a_r(i,j) * zh_soil(k+1)) + &
                    exp(tile_hv%b_r(i,j) * zh_soil(k+1)) - &
                    exp(tile_hv%a_r(i,j) * zh_soil(k  )) - &
                    exp(tile_hv%b_r(i,j) * zh_soil(k  )))

                root_sum_lv = root_sum_lv + tile_lv%root_frac(i,j,k)
                root_sum_hv = root_sum_hv + tile_hv%root_frac(i,j,k)
            end do

            ! Make sure that the fractions sum to one.
            tile_lv%root_frac(i,j,1) = 1. - root_sum_lv
            tile_hv%root_frac(i,j,1) = 1. - root_sum_hv
        end do
    end do

end subroutine calc_root_fractions

!
! Iterative Rib -> Obukhov length solver
!
function calc_obuk_dirichlet(L_in, du, db_in, zsl, z0m, z0h) result(res)
    use modglobal, only : fkar
    implicit none
    real, intent(in) :: L_in, du, db_in, zsl, z0m, z0h

    integer :: m, n, nlim
    real :: res, L, db, Lmax, Ri, L0, Lstart, Lend, fx0, fxdif

    m = 0
    nlim = 10
    Lmax = 1e10
    L = L_in
    db = db_in

    ! The solver does not have a solution for large Ri numbers,
    ! i.e. the `fx` equation below has no zero crossing.
    ! The limit of 0.13 typically results in a minimum (positive) L of ~1.
    ! IFS caps z/L at 5, so for a typical zsl at ~L=2.
    !Ri = fkar * db * zsl / du**2
    !if (Ri > 0.13) then
    !    print*,'WARNING: Ri out of range, returning L=1'
    !    res = 1
    !    return
    !end if

    ! Avoid buoyancy difference of zero:
    if (db >= 0) then
        db = max(db, 1e-9)
    else
        db = min(db, -1e-9)
    end if

    ! Allow for one restart of iterative procedure:
    do while (m <= 1)
        ! if L and db are of different sign, or the last calculation did not converge,
        ! the stability has changed and the procedure needs to be reset
        if (L*db <= 0) then
            nlim = 200
            if (db >= 0) then
                L = 1e-9
            else
                L = -1e-9
            end if
        end if

        ! Make sure the iteration starts
        if (db >= 0) then
            L0 = 1e30
        else
            L0 = -1e30
        end if

        ! Exit on convergence or on iteration count
        n = 0
        fxdif = 1
        do while (abs((L - L0) / L0) > 0.001 .and. n < nlim .and. abs(L) < Lmax)
            L0     = L
            Lstart = L - 0.001*L
            Lend   = L + 0.001*L

            fx0    = fx(zsl, L0, du, db, z0h, z0m)
            fxdif  = (fx(zsl, Lend, du, db, z0h, z0m) - fx(zsl, Lstart, du, db, z0h, z0m)) / (Lend - Lstart)
            L      = L - fx0/fxdif
            n      = n+1
        end do

        if (n < nlim .and. abs(L) < Lmax) then
            ! Convergence has been reached
            res = L
            return
        else
            ! Convergence has not been reached, procedure restarted once
            L = 1e-9
            m = m+1
            nlim = 200
        end if
    end do

    if (m > 1) then
        print*,'WARNING: convergence has not been reached in Obukhov length iteration'
        print*,'Input: ', L_in, du, db_in, zsl, z0m, z0h
        !stop
        res = 1e-9
        return
    end if

end function calc_obuk_dirichlet

pure function fx(zsl, L, du, db, z0h, z0m) result(res)
    implicit none
    real, intent(in) :: zsl, L, du, db, z0h, z0m
    real :: res, fkar
    fkar = 0.4

    res = zsl/L - fkar * zsl * db * fh(zsl, z0h, L) / (du * fm(zsl, z0m, L))**2
end function fx

pure function fm(zsl, z0m, L) result(res)
    use modglobal, only : fkar
    use modsurface, only : psim
    implicit none
    real, intent(in) :: zsl, z0m, L
    real :: res

    res = fkar / (log(zsl/z0m) - psim(zsl/L) + psim(z0m/L))
end function fm

pure function fh(zsl, z0h, L) result(res)
    use modglobal, only : fkar
    use modsurface, only : psih
    implicit none
    real, intent(in) :: zsl, z0h, L
    real :: res

    res = fkar / (log(zsl/z0h) - psih(zsl/L) + psih(z0h/L))
end function fh

!
! Convert soil hydraulic head to soil water content, using van Genuchten parameterisation.
!
pure function psi_to_theta(theta_res, theta_sat, vg_a, vg_n, vg_m, psi) result(res)
    implicit none
    real, intent(in) :: theta_res, theta_sat, vg_a, vg_n, vg_m, psi
    real :: res

    res = theta_res + (theta_sat - theta_res) * (1. / (1.+ abs(vg_a * psi)**vg_n))**vg_m
end function psi_to_theta

!
! Convert soil water content to hydraulic head, using van Genuchten parameterisation.
!
pure function theta_to_psi(theta_res, theta_sat, vg_a, vg_n, vg_m, theta) result(res)
    implicit none
    real, intent(in) :: theta_res, theta_sat, vg_a, vg_n, vg_m, theta
    real :: res

    res = -(vg_a**(-vg_n) * (-1. + ((theta_res - theta_sat)/(theta_res - theta))**(1./vg_m)))**(1./vg_n)
end function theta_to_psi

!
! Calculate hydraulic diffusivity using van Genuchten parameterisation.
!
pure function calc_diffusivity_vg( &
        theta_norm, vg_a, vg_l, vg_m, lambda_sat, theta_sat, theta_res) result(res)
    implicit none
    real, intent(in) :: theta_norm, vg_a, vg_l, vg_m, lambda_sat, theta_sat, theta_res
    real :: res, res_old
    real :: part1, part2, part3, part4, part5

    
    
    part1 = (1. - vg_m)
    part2 = lambda_sat / (vg_a * vg_m * (theta_sat - theta_res))
    part3 = theta_norm ** (vg_l - (1. / vg_m))
    
    if((theta_norm ** (1. / vg_m))>=1.0) then  	!some values above 1 are appeared and the part4 and part5 turns into 
    		  							      	!a complex number (due to rising a negative value to a fraction power) that couses segmentation fault:
    											!here I put 0.9 to avoid negative values !Arseni
    											
    	part4 = (1. - 0.9) ** (-vg_m)
    	part5 = (1. - 0.9) ** vg_m - 2. 
    	
    else
    
    	part4 = (1. - theta_norm ** (1. / vg_m)) ** (-vg_m)
    	part5 = (1. - theta_norm ** (1. / vg_m)) ** vg_m - 2.
    endif
    
    
    

    res = part1 * part2 * part3 * (part4 + part5)
            
    
end function calc_diffusivity_vg

!
! Calculate hydraulic conductivity using van Genuchten parameterisation.
!
pure function calc_conductivity_vg(theta_norm, vg_l, vg_m, gamma_sat) result(res)
    implicit none
    real, intent(in) :: theta_norm, vg_l, vg_m, gamma_sat
    real :: res

	
    if((theta_norm ** (1. / vg_m))>=1.0) then  	!some values above 1 are appeared and it turns into 
    		  							      	!a complex number (due to rising a negative value to a fraction power) that couses segmentation fault:
    											!here I put 0.9 to avoid negative values !Arseni
    
    	res = gamma_sat * theta_norm**vg_l * ( 1.- (1.-0.9)**vg_m )**2.
    else
    
    	res = gamma_sat * theta_norm**vg_l * ( 1.- (1.-theta_norm**(1./vg_m))**vg_m )**2.
    	
    endif

end function calc_conductivity_vg

!
! Check NetCDF calls
!
subroutine check(status)
    integer, intent (in) :: status
    if(status /= nf90_noerr) then
        print *,'NetCDF error: ', trim(nf90_strerror(status))
        stop
    end if
end subroutine check

end module modlsm