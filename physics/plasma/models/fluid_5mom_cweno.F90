
! This file is part of the muphyII multiphysics plasma simulation project.
! 
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at https://mozilla.org/MPL/2.0/.

! Centrally weighted essentially non-oscillating (CWENO) method

! fluid(:,:,:,1) = n
! fluid(:,:,:,2:4) = n u
! fluid(:,:,:,5) = eps = tr(Eps_ij) = 3 n k T / m + n u \dot u (scalar energy density)
! Eps_ij = \int v_i v_j f dv (energy density tensor)
! P_ij = m \int (v_i - u_i) (v_j - u_j) f dv (pressure tensor)
! p = tr(p_ij)/N (scalar pressure; N dimensionality)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CWENO part of Euler step t^n -> t^n+1, to be used in Runge-Kutta                     !
! CWENO follows Shu-Osher, with Kurganov-Levwy reconstruction                          !
!                                                                                      !
! General remarks:                                                                     !
! - Single-fluid calculations (e.g. j = j_s needs to be summed over to obtain fluid-j) !
! - E and B are presumed to hold the values at the right time t^{n+a} for this step    !
! - fluid holds values at t_n, fluid_interim holds the result from the previous step   !
! - j will be filled with face-centered values obtained from fluxes                    !
! - alpha and beta are RK3 parameters, see Shu-Osher for details                       !
!                                                                                      !
! Arguments:                                                                           !
!     dimX          : Dimensions of domain [nx-1,ny-1,nz-1] excluding ghostcells       !
!     BD            : Number of ghostcells in each direction                           !
!     fluid         : aray of [n,nu,nuu] at time t^n, where nu = {nu_x,nu_y,nu_z}      !
!     fluid_interim : result from last RK3 step, same form as 'fluid'                  !
!     j             : Current density, to be calculated at face centers from flux and  !
!                     used in Maxwell's equations, j = [j_x,j_y,j_z]                   !
!     E             : Electric field, E = [E_x,E_y,E_z]                                !
!     B             : Magnetic field, B = [B_x,B_y,B_z]                                !
!     alpha0, alpha1, beta                                                             !
!                   : Runge-Kutta prefactors, see 'Result'                             !
!     dt            : Time step width, dt = t^{n+1} - t^n                              !
!                                                                                      !
! Note that (n, nu, nuu) are the fluid fields (density, velocity, energy) as           !
! non-normalized integrals, using fluid density n, fluid velocity u                    !
!                                                                                      !
! Result:                                                                              !
!     fluid_interim = alpha0 * fluid + alpha1 * fluid_interim + beta * dt * RHS(fluid) !
!                                                                                      !
! Note that two ghostcells are needed in the flux calculation and reconstruction.      !
!                                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine five_moment_single_step_rk3(fluid, fluid_interim, E, B, q_m, alpha0, alpha1, beta, &
                   cweno_limiter, cweno_epsilon, dimX, BD, dx, dt, dimensionality_x, dimensionality_v)
    use definitions
    implicit none

    integer,  intent(in) :: dimX(1:3), BD(3), dimensionality_x, dimensionality_v
    real(kind=PRC),  intent(in) :: alpha0, alpha1, cweno_limiter, cweno_epsilon, beta, dx(1:3), dt, q_m

    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:5), intent(in) :: fluid
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:5), intent(inout) :: fluid_interim
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:3), intent(in) :: E, B

    real(kind=PRC), allocatable, dimension(:,:,:,:,:,:) :: num_flux
    real(kind=PRC), allocatable, dimension(:,:,:,:,:,:) :: reco
    real(kind=PRC), allocatable, dimension(:,:,:,:,:) :: flux
    real(kind=PRC), allocatable, dimension(:,:,:,:) :: total_flux, source, speed

    integer :: dir, fix, x, y, z

    allocate(num_flux(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3, 1:2))
    allocate(reco(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3, 1:2))
    allocate(flux(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3))
    allocate(total_flux(0:dimX(1), 0:dimX(2), 0:dimX(3), 1:5))
    allocate(source(0:dimX(1), 0:dimX(2), 0:dimX(3), 1:5))
    allocate(speed(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:3))

    !$acc enter data create(num_flux, reco, flux, total_flux, source, speed)

    call five_moment_cweno_reconstruct(dimX, dimensionality_x, BD, fluid_interim, reco, 5, cweno_epsilon)

    call five_moment_calc_flux(dimX, dimensionality_x, dimensionality_v, reco, num_flux, speed, cweno_limiter)

    !$acc data present(num_flux, reco, flux, speed, E, B, fluid_interim, total_flux, source)

    do dir = 1,dimensionality_x
        !$acc parallel
#ifdef ACC_2D_OPT
        !$acc loop gang collapse(3)
        do fix = 1,5 ! field index
            do z = 0,dimX(3)+1
                do y = 0,dimX(2)+1
                    !$acc loop vector
                    do x = 0,dimX(1)+1
#else
        !$acc loop gang collapse(2)
        do fix = 1,5 ! field index
            do z = 0,dimX(3)+1
                !$acc loop vector collapse(2)
                do y = 0,dimX(2)+1
                    do x = 0,dimX(1)+1
#endif
                        flux(x,y,z,fix,dir) = .5_PRC*( num_flux(x,y,z,fix,dir,2) + num_flux(x,y,z,fix,dir,1) &
                                + speed(x,y,z,dir) * ( reco(x,y,z,fix,dir,1) - reco(x,y,z,fix,dir,2) ) ) ! tht, 5.30 ! 1: _minus, 2: _plus
                    end do
                end do
            end do
        end do
        !$acc end parallel
    end do

    !$acc parallel
#ifdef ACC_2D_OPT
    !$acc loop gang collapse(3)
    do fix = 1,5
        do z = 0,dimX(3)
            do y = 0,dimX(2)
                !$acc loop vector
                do x = 0,dimX(1)
#else
    !$acc loop gang collapse(2)
    do fix = 1,5
        do z = 0,dimX(3)
            !$acc loop vector collapse(2)
            do y = 0,dimX(2)
                do x = 0,dimX(1)
#endif
                    total_flux(x,y,z,fix) = (flux(x+1,y,z,fix,1) - flux(x,y,z,fix,1)) / dx(1)
                end do
            end do
        end do
    end do
    !$acc end parallel

    if (dimensionality_x > 1) then
        !$acc parallel
#ifdef ACC_2D_OPT
        !$acc loop gang collapse(3)
        do fix = 1,5
            do z = 0,dimX(3)
                do y = 0,dimX(2)
                    !$acc loop vector
                    do x = 0,dimX(1)
#else
        !$acc loop gang collapse(2)
        do fix = 1,5
            do z = 0,dimX(3)
                !$acc loop vector collapse(2)
                do y = 0,dimX(2)
                    do x = 0,dimX(1)
#endif
                        total_flux(x,y,z,fix) = total_flux(x,y,z,fix) + (flux(x,y+1,z,fix,2) - flux(x,y,z,fix,2)) / dx(2)
                    end do
                end do
            end do
        end do
        !$acc end parallel
    end if

    if (dimensionality_x > 2) then
        !$acc parallel
#ifdef ACC_2D_OPT
        !$acc loop gang collapse(3)
        do fix = 1,5
            do z = 0,dimX(3)
                do y = 0,dimX(2)
                    !$acc loop vector
                    do x = 0,dimX(1)
#else
        !$acc loop gang collapse(2)
        do fix = 1,5
            do z = 0,dimX(3)
                !$acc loop vector collapse(2)
                do y = 0,dimX(2)
                    do x = 0,dimX(1)
#endif
                        total_flux(x,y,z,fix) = total_flux(x,y,z,fix) + (flux(x,y,z+1,fix,3) - flux(x,y,z,fix,3)) / dx(3)
                    end do
                end do
            end do
        end do
        !$acc end parallel
    end if

    !$acc parallel
#ifdef ACC_2D_OPT
    !$acc loop gang collapse(2)
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            !$acc loop vector
            do x = 0,dimX(1)
#else
    !$acc loop gang
    do z = 0,dimX(3)
        !$acc loop vector collapse(2)
        do y = 0,dimX(2)
            do x = 0,dimX(1)
#endif
                source(x,y,z,1)  = 0._PRC
                ! q_m * ( Ex * n + nuy * Bz - nuz * By )
                source(x,y,z,2)  = q_m * (E(x,y,z,1)*fluid_interim(x,y,z,1) + fluid_interim(x,y,z,3)*B(x,y,z,3) - fluid_interim(x,y,z,4)*B(x,y,z,2))
                ! q_m * ( Ey * n + nuz * Bx - nux * Bz )
                source(x,y,z,3)  = q_m * (E(x,y,z,2)*fluid_interim(x,y,z,1) + fluid_interim(x,y,z,4)*B(x,y,z,1) - fluid_interim(x,y,z,2)*B(x,y,z,3))
                ! q_m * ( Ez * n + nux * By - nuy * Bx )
                source(x,y,z,4)  = q_m * (E(x,y,z,3)*fluid_interim(x,y,z,1) + fluid_interim(x,y,z,2)*B(x,y,z,2) - fluid_interim(x,y,z,3)*B(x,y,z,1))
                ! 2._PRC * q_m * ( Ex * nux + Ey * nuy + Ez * nuz )
                source(x,y,z,5) = 2._PRC*q_m * (E(x,y,z,1)*fluid_interim(x,y,z,2) + E(x,y,z,2)*fluid_interim(x,y,z,3) + E(x,y,z,3)*fluid_interim(x,y,z,4))
            end do
        end do
    end do
    !$acc end parallel

    ! runge-kutta-3 step (Shu-Osher)
    !$acc parallel
#ifdef ACC_2D_OPT
    !$acc loop gang collapse(3)
    do fix = 1,5
        do z = 0,dimX(3)
            do y = 0,dimX(2)
                !$acc loop vector
                do x = 0,dimX(1)
#else
    !$acc loop gang collapse(2)
    do fix = 1,5
        do z = 0,dimX(3)
            !$acc loop vector collapse(2)
            do y = 0,dimX(2)
                do x = 0,dimX(1)
#endif
                    fluid_interim(x,y,z,fix) = &
                        alpha0 * fluid(x,y,z,fix) + &
                        alpha1 * fluid_interim(x,y,z,fix) + &
                        beta * dt * (source(x,y,z,fix) - total_flux(x,y,z,fix))
                end do
            end do
        end do
    end do
    !$acc end parallel

    !$acc end data

    !$acc exit data delete(num_flux, reco, flux, total_flux, source, speed)

    deallocate(num_flux)
    deallocate(reco)
    deallocate(flux)
    deallocate(total_flux)
    deallocate(source)
    deallocate(speed)

end subroutine five_moment_single_step_rk3


subroutine five_moment_calc_flux(dimX, dimensionality_x, dimensionality_v, reco, flux, speed, cweno_limiter)
    use definitions
    implicit none
    integer,  intent(in) :: dimX(3), dimensionality_x, dimensionality_v
    real(kind=PRC), intent(in) :: cweno_limiter
    real(kind=PRC), dimension(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3, 1:2), intent(in)    :: reco
    real(kind=PRC), dimension(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3, 1:2), intent(inout) :: flux
    real(kind=PRC), dimension(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1,      1:3),      intent(inout) :: speed

    integer :: dir, side, x, y, z
    real(kind=PRC) :: r_nunu, r_n_inv, r_nux, r_nuy, r_nuz, r_nudir, r_nuu

    !$acc data present(reco, flux, speed)

    !$acc parallel
    !$acc loop seq
    do side = 1,2
        !$acc loop seq
        do dir = 1,dimensionality_x
#ifdef ACC_2D_OPT
            !$acc loop gang collapse(2) independent
            do z = 0,dimX(3)+1
                do y = 0,dimX(2)+1
                    !$acc loop vector independent
                    do x = 0,dimX(1)+1
#else
            !$acc loop gang independent
            do z = 0,dimX(3)+1
                !$acc loop vector collapse(2) independent
                do y = 0,dimX(2)+1
                    do x = 0,dimX(1)+1
#endif
                        ! what follows represents get_flux() in muphy1
                        r_nux = reco(x,y,z,2,dir,side) ! just for cache
                        r_nuy = reco(x,y,z,3,dir,side)
                        r_nuz = reco(x,y,z,4,dir,side)
                        r_nuu = reco(x,y,z,5,dir,side)
                        if (dir==1) then
                           r_nudir = r_nux
                        else if (dir==2) then
                           r_nudir = r_nuy
                        else
                           r_nudir = r_nuz
                        end if
                        r_nunu = r_nux**2 + r_nuy**2 + r_nuz**2
                        r_n_inv = 1._PRC/max(reco(x,y,z, 1, dir, side), cweno_limiter)

                        flux(x,y,z,1,dir,side) = r_nudir  ! flux(n) = reco(u_dir) , etc.

                        flux(x,y,z,2,dir,side) = r_nux * r_nudir * r_n_inv
                        flux(x,y,z,3,dir,side) = r_nuy * r_nudir * r_n_inv
                        flux(x,y,z,4,dir,side) = r_nuz * r_nudir * r_n_inv

                        flux(x,y,z,1+dir,dir,side) = flux(x,y,z,1+dir,dir,side) + &
                             (r_nuu - r_nunu * r_n_inv)/dimensionality_v

                        flux(x,y,z,5,dir,side) = r_nudir * r_n_inv * &
                             ((2._PRC+dimensionality_v)*r_nuu - 2._PRC*r_nunu*r_n_inv)/dimensionality_v

                        speed(x,y,z,dir) = abs(r_nudir) * r_n_inv + &
                             sqrt(max(r_nuu - r_nunu*r_n_inv, 0._PRC) * r_n_inv)
                    end do
                end do
            end do
        end do
    end do
    !$acc end parallel
    !$acc end data

end subroutine five_moment_calc_flux


subroutine five_moment_calc_un_from_flux(un, fluid, cweno_limiter, cweno_epsilon, dimX, BD, dimensionality_v)
    use definitions
    implicit none

    integer,  intent(in) :: dimX(1:3), BD(3), dimensionality_v
    real(kind=PRC), intent(in) :: cweno_limiter, cweno_epsilon
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:5), intent(in) :: fluid
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:3), intent(inout) :: un

    real(kind=PRC), allocatable, dimension(:,:,:,:,:,:) :: num_flux
    real(kind=PRC), allocatable, dimension(:,:,:,:,:,:) :: reco
    real(kind=PRC), allocatable, dimension(:,:,:,:) :: speed
    integer :: dir, x, y, z

    allocate(num_flux(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3, 1:2))
    allocate(reco(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3, 1:2))
    allocate(speed(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:3))

    !$acc enter data create(num_flux, reco, speed)

    !$acc data present(num_flux, reco, speed, un, fluid)

    call five_moment_cweno_reconstruct(dimX, 3, BD, fluid, reco, 5, cweno_epsilon)

    ! TODO there is some overhead due to reusing the function
    call five_moment_calc_flux(dimX, 3, dimensionality_v, reco, num_flux, speed, cweno_limiter)

    do dir = 1,3
        !$acc parallel
#ifdef ACC_2D_OPT
        !$acc loop gang collapse(2)
        do z = 0,dimX(3)+1
            do y = 0,dimX(2)+1
                !$acc loop vector
                do x = 0,dimX(1)+1
#else
        !$acc loop gang
        do z = 0,dimX(3)+1
            !$acc loop vector collapse(2)
            do y = 0,dimX(2)+1
                do x = 0,dimX(1)+1
#endif
                    un(x,y,z,dir) = & ! tht, 5.30 ! 1: _minus, 2: _plus
                            .5_PRC*( num_flux(x,y,z,1,dir,2) + num_flux(x,y,z,1,dir,1) &
                            + speed(x,y,z,dir) * ( reco(x,y,z,1,dir,1) - reco(x,y,z,1,dir,2) ) )
                end do
            end do
        end do
        !$acc end parallel
    end do

    !$acc end data

    !$acc exit data delete(num_flux, reco, speed)

    deallocate(num_flux)
    deallocate(reco)
    deallocate(speed)

end subroutine five_moment_calc_un_from_flux


subroutine five_moment_cweno_reconstruct(dimX, dimensionality_x, BD, fields, reco, num_fields, cweno_epsilon)
    use definitions
    implicit none

    integer,  intent(in) :: dimX(1:3), dimensionality_x, BD(3), num_fields
    real(kind=PRC), intent(in) :: cweno_epsilon
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),num_fields), &
            intent(in) :: fields
    real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, num_fields, 1:3, 1:2), &
            intent(inout) :: reco

    real(kind=PRC), allocatable, dimension(:,:,:,:) :: a, b
    real(kind=PRC) :: data_down, data_up, w_c, w_l, w_r, w, c1
    integer :: dir, fix, x, y, z
    integer, dimension(1:3,1:3) :: shift

    allocate(a(-1:dimX(1)+1, -1:dimX(2)+1, -1:dimX(3)+1, num_fields))
    allocate(b(-1:dimX(1)+1, -1:dimX(2)+1, -1:dimX(3)+1, num_fields))
    !$acc enter data create(a,b)

    shift = reshape((/1,0,0, 0,1,0, 0,0,1/), (/3,3/))

    !$acc data present(a, b, fields, reco) copyin(shift)

    do dir = 1,dimensionality_x
        !$acc parallel
#ifdef ACC_2D_OPT
        !$acc loop gang collapse(3)
        do fix = 1,num_fields
            do z = -1,dimX(3)+1
                do y = -1,dimX(2)+1
                    !$acc loop vector
                    do x = -1,dimX(1)+1
#else
        !$acc loop gang collapse(2)
        do fix = 1,num_fields
            do z = -1,dimX(3)+1
                !$acc loop vector collapse(2)
                do y = -1,dimX(2)+1
                    do x = -1,dimX(1)+1
#endif
                        data_down = fields(x-shift(dir,1),y-shift(dir,2),z-shift(dir,3),fix)
                        data_up   = fields(x+shift(dir,1),y+shift(dir,2),z+shift(dir,3),fix)
                        c1 = (data_down - 2._PRC*fields(x,y,z,fix) + data_up)
                        w_l = .25_PRC / (cweno_epsilon + (data_down - fields(x,y,z,fix))**2)**2
                        w_r = .25_PRC / (cweno_epsilon + (fields(x,y,z,fix) - data_up)**2)**2
                        w_c = .5_PRC  / (cweno_epsilon + (data_up - data_down)**2 * .25_PRC + c1**2 * 13._PRC/3._PRC)**2
                        w = w_l + w_r + w_c
                        w_l = w_l / w
                        w_r = w_r / w
                        w_c = w_c / w
                        c1 = w_c * c1
                        a(x,y,z,fix) = fields(x,y,z,fix) + c1/6._PRC
                        b(x,y,z,fix) = .5_PRC * (((.5_PRC*w_c + w_r)*data_up + (-w_r + w_l)*fields(x,y,z,fix) - (.5_PRC*w_c + w_l)*data_down))
                    end do
                end do
            end do
        end do
        !$acc end parallel

        !$acc parallel
#ifdef ACC_2D_OPT
        !$acc loop gang collapse(3)
        do fix = 1,num_fields
            do z = 0,dimX(3)+1
                do y = 0,dimX(2)+1
                    !$acc loop vector
                    do x = 0,dimX(1)+1
#else
        !$acc loop gang collapse(2)
        do fix = 1,num_fields
            do z = 0,dimX(3)+1
                !$acc loop vector collapse(2)
                do y = 0,dimX(2)+1
                    do x = 0,dimX(1)+1
#endif
                        reco(x,y,z,fix,dir,1) = a(x-shift(dir,1),y-shift(dir,2),z-shift(dir,3),fix) + &
                                                b(x-shift(dir,1),y-shift(dir,2),z-shift(dir,3),fix)
                        reco(x,y,z,fix,dir,2) = a(x,y,z,fix) - b(x,y,z,fix)
                    end do
                end do
            end do
        end do
        !$acc end parallel
    end do

    !$acc end data

    !$acc exit data delete(a,b)
    deallocate(a)
    deallocate(b)

end subroutine five_moment_cweno_reconstruct


function five_moment_max_dt(fluid, cfl, cweno_limiter, dimX, BD, dx, dimensionality_v) result(dt)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), dimensionality_v
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),5) :: fluid
    real(kind=PRC), intent(in) :: dx(3), cfl, cweno_limiter

    real(kind=PRC), allocatable, dimension(:,:,:,:) :: v
    real(kind=PRC) :: n_inv, nunu, v_max
    integer :: x, y, z, dir
    real(kind=PRC) :: dt

    allocate(v(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3))
    !$acc enter data create(v)

    !$acc data present(fluid, v)
    !$acc parallel
#ifdef ACC_2D_OPT
    !$acc loop gang collapse(2)
    do z = -BD(3),dimX(3)+BD(3)
        do y = -BD(2),dimX(2)+BD(2)
            !$acc loop vector
            do x = -BD(1),dimX(1)+BD(1)
#else
    !$acc loop gang
    do z = -BD(3),dimX(3)+BD(3)
        !$acc loop vector collapse(2)
        do y = -BD(2),dimX(2)+BD(2)
            do x = -BD(1),dimX(1)+BD(1)
#endif
                n_inv = 1._PRC/max(fluid(x,y,z,1), cweno_limiter)
                nunu = fluid(x,y,z,2)**2 + fluid(x,y,z,3)**2 + fluid(x,y,z,4)**2

                do dir = 1,dimensionality_v
                    v(x,y,z,dir) = abs(fluid(x,y,z,dir+1)) * n_inv &
                         + sqrt((2._PRC+dimensionality_v)/dimensionality_v * max(fluid(x,y,z,5) - nunu * n_inv, 0._PRC) * n_inv)
                end do
            end do
        end do
    end do
    !$acc end parallel

    !$acc kernels
    v_max = maxval(v)
    !$acc end kernels

    !$acc end data

    dt = minval( dx / v_max ) * cfl

    !$acc exit data delete(v)
    deallocate(v)

end function five_moment_max_dt

