
! This file is part of the muphyII multiphysics plasma simulation project.
! 
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at https://mozilla.org/MPL/2.0/.

! TODO MHD to OpenAcc, unite with other fluid files

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

subroutine five_moment_mhd_single_step_rk3(fluid, fluid_interim, E, B, J, &
        alpha0, alpha1, beta, &
        cweno_limiter, cweno_epsilon, dimX, BD, dx, dt, dimensionality_x)
    use definitions
    implicit none
    interface
        subroutine five_moment_mhd_calc_single_step(fluid_interim, num_flux, reco, &
                flux, total_flux, speed, &
                cweno_limiter, cweno_epsilon, dimX, BD, dx, dimensionality_x)
            use definitions
            implicit none
            integer,  intent(in) :: dimX(1:3), BD(3), dimensionality_x
            real(kind=PRC),  intent(in) :: cweno_limiter, cweno_epsilon, dx(1:3)
            real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:5), target :: fluid_interim
            real(kind=PRC), intent(out), dimension(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3, 1:2)         :: num_flux
            real(kind=PRC), intent(out), dimension(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3, 1:2), target :: reco
            real(kind=PRC), intent(out), dimension(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3)              :: flux
            real(kind=PRC), intent(out), dimension(0:dimX(1), 0:dimX(2), 0:dimX(3), 1:5)                         :: total_flux
            real(kind=PRC), intent(out), dimension(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:3)                   :: speed
        end subroutine five_moment_mhd_calc_single_step
        subroutine five_moment_mhd_calc_source(source, E, B, J, rho, dimX, BD)
            use definitions
            implicit none
            integer :: dimX(3), BD(3)
            real(kind=PRC), dimension(0:dimX(1),0:dimX(2),0:dimX(3),1:5), intent(out) :: source
            real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:3), intent(in), target :: E, B, J
            real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)),     intent(in) :: rho
        end subroutine five_moment_mhd_calc_source
    end interface

    integer,  intent(in) :: dimX(1:3), BD(3), dimensionality_x
    real(kind=PRC),  intent(in) :: alpha0, alpha1, cweno_limiter, cweno_epsilon, beta, dx(1:3), dt

    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:5), intent(in) :: fluid
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:5), intent(inout), target :: fluid_interim
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:3), intent(in), target :: E, B, J

    real(kind=PRC), allocatable, dimension(:,:,:,:,:,:) :: num_flux
    real(kind=PRC), allocatable, dimension(:,:,:,:,:,:), target :: reco
    real(kind=PRC), allocatable, dimension(:,:,:,:,:) :: flux
    real(kind=PRC), allocatable, dimension(:,:,:,:) :: total_flux, source, speed

    allocate(num_flux(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3, 1:2))
    allocate(reco(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3, 1:2))
    allocate(flux(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3))
    allocate(total_flux(0:dimX(1), 0:dimX(2), 0:dimX(3), 1:5))
    allocate(speed(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:3))
    allocate(source(0:dimX(1), 0:dimX(2), 0:dimX(3), 1:5))

    call five_moment_mhd_calc_single_step(fluid_interim, num_flux, reco, &
            flux, total_flux, speed, &
            cweno_limiter, cweno_epsilon, dimX, BD, dx, dimensionality_x)
    call five_moment_mhd_calc_source(source, E, B, J, fluid_interim(:,:,:,1), dimX, BD)
    call five_moment_mhd_apply_single_step(fluid, fluid_interim, source, total_flux, &
            alpha0, alpha1, beta, dt, dimX, BD)

    deallocate(num_flux)
    deallocate(reco)
    deallocate(flux)
    deallocate(total_flux)
    deallocate(speed)
    deallocate(source)
end subroutine five_moment_mhd_single_step_rk3


subroutine five_moment_mhd_calc_single_step(fluid_interim, num_flux, reco, &
        flux, total_flux, speed, &
        cweno_limiter, cweno_epsilon, dimX, BD, dx, dimensionality_x)
    use definitions
    implicit none
    interface
        subroutine five_moment_mhd_calc_flux(dimX, dimensionality_x, reco, flux, speed, cweno_limiter)
            use definitions
            implicit none
            integer,  intent(in) :: dimX(3), dimensionality_x
            real(kind=PRC), intent(in) :: cweno_limiter
            real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3, 1:2), intent(in), target :: reco
            real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3, 1:2), intent(inout)        :: flux
            real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1,      1:3),      intent(inout)        :: speed
        end subroutine five_moment_mhd_calc_flux

        subroutine five_moment_mhd_cweno_reconstruct(dimX, dimensionality_x, BD, fields, reco, num_fields, cweno_epsilon)
            use definitions
            implicit none
            integer,  intent(in) :: dimX(1:3), dimensionality_x, BD(3), num_fields
            real(kind=PRC), intent(in) :: cweno_epsilon
            real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),num_fields),&
                    intent(in),  target :: fields
            real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, num_fields, 1:3, 1:2),&
                    intent(inout) :: reco
        end subroutine five_moment_mhd_cweno_reconstruct
    end interface

    integer,  intent(in) :: dimX(1:3), BD(3), dimensionality_x
    real(kind=PRC),  intent(in) :: cweno_limiter, cweno_epsilon, dx(1:3)
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:5), target :: fluid_interim
    real(kind=PRC), intent(out), dimension(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3, 1:2)         :: num_flux
    real(kind=PRC), intent(out), dimension(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3, 1:2), target :: reco
    real(kind=PRC), intent(out), dimension(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3)              :: flux
    real(kind=PRC), intent(out), dimension(0:dimX(1), 0:dimX(2), 0:dimX(3), 1:5)                         :: total_flux
    real(kind=PRC), intent(out), dimension(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:3)                   :: speed

    integer :: dir, fix

    call five_moment_mhd_cweno_reconstruct(dimX, dimensionality_x, BD, fluid_interim, reco, 5, cweno_epsilon)

    call five_moment_mhd_calc_flux(dimX, dimensionality_x, reco, num_flux, speed, cweno_limiter)
    do dir = 1,dimensionality_x
        do fix = 1,5
            flux(:,:,:,fix,dir) = .5_PRC*( num_flux(:,:,:,fix,dir,2) + num_flux(:,:,:,fix,dir,1) &
                    + speed(:,:,:,dir) * ( reco(:,:,:,fix,dir,1) - reco(:,:,:,fix,dir,2) ) ) ! tht, 5.30 ! 1: _minus, 2: _plus
        end do
    end do

    ! what follows represents apply_flux() in muphy1
    ! tht, 5.29
    total_flux = (flux(1:dimX(1)+1,0:dimX(2),0:dimX(3),:,1) - flux(0:dimX(1),0:dimX(2),0:dimX(3),:,1)) / dx(1)
    if (dimensionality_x > 1) then
        total_flux = total_flux + (flux(0:dimX(1),1:dimX(2)+1,0:dimX(3),:,2) - flux(0:dimX(1),0:dimX(2),0:dimX(3),:,2)) / dx(2)
    end if
    if (dimensionality_x > 2) then
        total_flux = total_flux + (flux(0:dimX(1),0:dimX(2),1:dimX(3)+1,:,3) - flux(0:dimX(1),0:dimX(2),0:dimX(3),:,3)) / dx(3)
    end if
end subroutine five_moment_mhd_calc_single_step


subroutine five_moment_mhd_calc_source(source, E, B, J, rho, dimX, BD)
    use definitions
    implicit none

    integer :: dimX(3), BD(3)
    real(kind=PRC), dimension(0:dimX(1),0:dimX(2),0:dimX(3),1:5), intent(out) :: source
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:3), intent(in), target :: E, B, J
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)),     intent(in) :: rho

    real(kind=PRC), dimension(:,:,:), pointer :: Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz

    ! what follows represents get_background() in muphy1
    Ex  => E(0:dimX(1), 0:dimX(2), 0:dimX(3), 1)
    Ey  => E(0:dimX(1), 0:dimX(2), 0:dimX(3), 2)
    Ez  => E(0:dimX(1), 0:dimX(2), 0:dimX(3), 3)
    Bx  => B(0:dimX(1), 0:dimX(2), 0:dimX(3), 1)
    By  => B(0:dimX(1), 0:dimX(2), 0:dimX(3), 2)
    Bz  => B(0:dimX(1), 0:dimX(2), 0:dimX(3), 3)
    Jx  => J(0:dimX(1), 0:dimX(2), 0:dimX(3), 1)
    Jy  => J(0:dimX(1), 0:dimX(2), 0:dimX(3), 2)
    Jz  => J(0:dimX(1), 0:dimX(2), 0:dimX(3), 3)

    source(:,:,:, 1) = 0._PRC
    source(:,:,:, 2) = ( rho * Ex + Jy * Bz - Jz * By )
    source(:,:,:, 3) = ( rho * Ey + Jz * Bx - Jx * Bz )
    source(:,:,:, 4) = ( rho * Ez + Jx * By - Jy * Bx )
    source(:,:,:, 5) = 2._PRC * ( Jx * Ex + Jy * Ey + Jz * Ez )

end subroutine five_moment_mhd_calc_source


subroutine five_moment_mhd_apply_single_step(fluid, fluid_interim, source, total_flux, &
        alpha0, alpha1, beta, dt, dimX, BD)
    use definitions
    implicit none

    real(kind=PRC),  intent(in) :: alpha0, alpha1, beta, dt
    integer :: dimX(3), BD(3)

    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:5), intent(in) :: fluid
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:5), intent(inout) :: fluid_interim
    real(kind=PRC), dimension(0:dimX(1),0:dimX(2),0:dimX(3),1:5), intent(in) :: source, total_flux

    ! apply flux
    fluid_interim(                 0:dimX(1), 0:dimX(2), 0:dimX(3), :) = &
            alpha0 * fluid(        0:dimX(1), 0:dimX(2), 0:dimX(3), :) + &
            alpha1 * fluid_interim(0:dimX(1), 0:dimX(2), 0:dimX(3), :) + &
            beta * dt * (source - total_flux)
end subroutine five_moment_mhd_apply_single_step


! subroutine five_moment_mhd_calc_source_2fluid(source, fluid_interim, E, B, q_m)
!     use definitions
!     implicit none
!
!     real(kind=PRC) :: q_m
!      real(kind=PRC), dimension(0:dimX(1),0:dimX(2),0:dimX(3),1:5), intent(out) :: source
!      real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:5), intent(in), target :: fluid_interim     real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:3), intent(in), target :: E, B
!
!      real(kind=PRC), dimension(:,:,:), pointer :: n, nux, nuy, nuz, Ex, Ey, Ez, Bx, By, Bz
!
!      ! what follows represents get_background() in muphy1
!
!      ! use pointers for better readability
!      n   => fluid_interim(0:dimX(1), 0:dimX(2), 0:dimX(3), 1)
!      nux => fluid_interim(0:dimX(1), 0:dimX(2), 0:dimX(3), 2)
!      nuy => fluid_interim(0:dimX(1), 0:dimX(2), 0:dimX(3), 3)
!      nuz => fluid_interim(0:dimX(1), 0:dimX(2), 0:dimX(3), 4)
!
!      Ex  => E(0:dimX(1), 0:dimX(2), 0:dimX(3), 1)
!      Ey  => E(0:dimX(1), 0:dimX(2), 0:dimX(3), 2)
!      Ez  => E(0:dimX(1), 0:dimX(2), 0:dimX(3), 3)
!      Bx  => B(0:dimX(1), 0:dimX(2), 0:dimX(3), 1)
!      By  => B(0:dimX(1), 0:dimX(2), 0:dimX(3), 2)
!      Bz  => B(0:dimX(1), 0:dimX(2), 0:dimX(3), 3)
!
!      source(:,:,:, 1) = 0._PRC
!      source(:,:,:, 2) = q_m * ( Ex * n + nuy * Bz - nuz * By )
!      source(:,:,:, 3) = q_m * ( Ey * n + nuz * Bx - nux * Bz )
!      source(:,:,:, 4) = q_m * ( Ez * n + nux * By - nuy * Bx )
!      source(:,:,:, 5) = 2._PRC * q_m * ( Ex * nux + Ey * nuy + Ez * nuz )
!  end subroutine five_moment_mhd_calc_source_2fluid




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flux calculation from reconstruction                                              !
!                                                                                   !
! Arguments as in parent function                                                   !
!                                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine five_moment_mhd_calc_flux(dimX, dimensionality_x, reco, flux, speed, cweno_limiter)
    use definitions
    implicit none
    ! reconstruction has form
    ! [DATA, field, direction], that is,
    ! [DATA, [n, nux, nuy, nuz, nuu], [X,Y,Z]]
    integer,  intent(in) :: dimX(3), dimensionality_x
    real(kind=PRC), intent(in) :: cweno_limiter
    real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3, 1:2), intent(in), target :: reco
    real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3, 1:2), intent(inout)        :: flux
    real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1,      1:3),      intent(inout)        :: speed

    integer :: dir, side
    real(kind=PRC), allocatable, dimension(:,:,:) :: r_nunu, r_n_inv
    real(kind=PRC), dimension(:,:,:), pointer :: r_n, r_nux, r_nuy, r_nuz, r_nuu, r_nudir

    allocate(r_nunu(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1))
    allocate(r_n_inv(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1))

    speed = 0._PRC

    do side = 1,2
        do dir = 1,dimensionality_x
            ! use pointers for better readability
            r_n     => reco(:, :, :, 1,     dir, side)
            r_nux   => reco(:, :, :, 2,     dir, side)
            r_nuy   => reco(:, :, :, 3,     dir, side)
            r_nuz   => reco(:, :, :, 4,     dir, side)
            r_nudir => reco(:, :, :, 1+dir, dir, side)   ! field(1+dir) = nu(dir)
            r_nuu   => reco(:, :, :, 5,     dir, side)

            ! what follows represents get_flux() in muphy1
            r_nunu = r_nux**2 + r_nuy**2 + r_nuz**2
            r_n_inv = 1._PRC/max(r_n, cweno_limiter)

            flux(:,:,:,1,dir,side) = r_nudir  ! flux(n) = reco(u_dir) , etc.

            flux(:,:,:,2,dir,side) = r_nux * r_nudir * r_n_inv
            flux(:,:,:,3,dir,side) = r_nuy * r_nudir * r_n_inv
            flux(:,:,:,4,dir,side) = r_nuz * r_nudir * r_n_inv
            flux(:,:,:,1+dir,dir,side) = flux(:,:,:,1+dir,dir,side) + (r_nuu - r_nunu * r_n_inv)/3._PRC

            flux(:,:,:,5,dir,side) = r_nudir * r_n_inv/3._PRC * (5._PRC*r_nuu - 2._PRC*r_nunu*r_n_inv)

            speed(:,:,:,dir) = max( speed(:,:,:,dir), &
                    abs(r_nudir) * r_n_inv + sqrt(5._PRC/3._PRC * max(r_nuu - r_nunu*r_n_inv, 0._PRC) * r_n_inv) )
        end do
    end do

    deallocate(r_nunu)
    deallocate(r_n_inv)

end subroutine five_moment_mhd_calc_flux


subroutine five_moment_mhd_calc_un_from_flux(un, fluid, cweno_limiter, cweno_epsilon, dimX, BD)
    use definitions
    implicit none
    interface
        subroutine five_moment_mhd_calc_flux(dimX, dimensionality_x, reco, flux, speed, cweno_limiter)
            use definitions
            implicit none
            integer,  intent(in) :: dimX(3), dimensionality_x
            real(kind=PRC), intent(in) :: cweno_limiter
            real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3, 1:2), intent(in), target :: reco
            real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3, 1:2), intent(inout)        :: flux
            real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1,      1:3),      intent(inout)        :: speed
        end subroutine five_moment_mhd_calc_flux

        subroutine five_moment_mhd_cweno_reconstruct(dimX, dimensionality_x, BD, fields, reco, num_fields, cweno_epsilon)
            use definitions
            implicit none
            integer,  intent(in) :: dimX(1:3), dimensionality_x, BD(3), num_fields
            real(kind=PRC), intent(in) :: cweno_epsilon
            real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),num_fields),&
                    intent(in),  target :: fields
            real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, num_fields, 1:3, 1:2),&
                    intent(inout) :: reco
        end subroutine five_moment_mhd_cweno_reconstruct
    end interface

    integer,  intent(in) :: dimX(1:3), BD(3)
    real(kind=PRC), intent(in) :: cweno_limiter, cweno_epsilon
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:5), intent(in) :: fluid
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:3), intent(inout) :: un

    real(kind=PRC), allocatable, dimension(:,:,:,:,:,:) :: num_flux
    real(kind=PRC), allocatable, dimension(:,:,:,:,:,:) :: reco
    real(kind=PRC), allocatable, dimension(:,:,:,:) :: speed
    integer :: dir

    allocate(num_flux(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3, 1:2))
    allocate(reco(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:5, 1:3, 1:2))
    allocate(speed(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:3))

    call five_moment_mhd_cweno_reconstruct(dimX, 3, BD, fluid(:,:,:,1), reco, 5, cweno_epsilon)

    ! TODO there is some overhead due to reusing the function
    call five_moment_mhd_calc_flux(dimX, 3, reco, num_flux, speed, cweno_limiter)
    do dir = 1,3
        un(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1,dir) = & ! tht, 5.30 ! 1: _minus, 2: _plus
                .5_PRC*( num_flux(:,:,:,1,dir,2) + num_flux(:,:,:,1,dir,1) &
                + speed(:,:,:,dir) * ( reco(:,:,:,1,dir,1) - reco(:,:,:,1,dir,2) ) )
    end do

    deallocate(num_flux)
    deallocate(reco)
    deallocate(speed)

end subroutine five_moment_mhd_calc_un_from_flux


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate reconstruction                                                          !
!                                                                                   !
! Arguments as in parent function                                                   !
!                                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine five_moment_mhd_cweno_reconstruct(dimX, dimensionality_x, BD, fields, reco, num_fields, cweno_epsilon)
    use definitions
    implicit none

    integer,  intent(in) :: dimX(1:3), dimensionality_x, BD(3), num_fields
    real(kind=PRC), intent(in) :: cweno_epsilon
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),num_fields),&
            intent(in),  target :: fields
    real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, num_fields, 1:3, 1:2),&
            intent(inout) :: reco

    real(kind=PRC), dimension(:,:,:,:), pointer :: data, data_down, data_up
    real(kind=PRC), allocatable, dimension(:,:,:,:) :: w_c, w_l, w_r, w, c1, a, b

    integer :: dir
    integer, dimension(1:3,1:3) :: shift

    allocate(w_c(-1:dimX(1)+1, -1:dimX(2)+1, -1:dimX(3)+1, num_fields))
    allocate(w_l(-1:dimX(1)+1, -1:dimX(2)+1, -1:dimX(3)+1, num_fields))
    allocate(w_r(-1:dimX(1)+1, -1:dimX(2)+1, -1:dimX(3)+1, num_fields))
    allocate(w  (-1:dimX(1)+1, -1:dimX(2)+1, -1:dimX(3)+1, num_fields))
    allocate(c1 (-1:dimX(1)+1, -1:dimX(2)+1, -1:dimX(3)+1, num_fields))
    allocate(a  (-1:dimX(1)+1, -1:dimX(2)+1, -1:dimX(3)+1, num_fields))
    allocate(b  (-1:dimX(1)+1, -1:dimX(2)+1, -1:dimX(3)+1, num_fields))

    shift = reshape((/1,0,0, 0,1,0, 0,0,1/), (/3,3/))

    data => fields(-1:dimX(1)+1, -1:dimX(2)+1, -1:dimX(3)+1, :)
    do dir = 1,dimensionality_x
        data_down => fields(-1-shift(dir,1):dimX(1)+1-shift(dir,1), -1-shift(dir,2):dimX(2)+1-shift(dir,2), &
                -1-shift(dir,3):dimX(3)+1-shift(dir,3), :)
        data_up   => fields(-1+shift(dir,1):dimX(1)+1+shift(dir,1), -1+shift(dir,2):dimX(2)+1+shift(dir,2), &
                -1+shift(dir,3):dimX(3)+1+shift(dir,3), :)
        c1 = ( data_down - 2._PRC* data + data_up )
        call five_moment_mhd_compute_weight_alpha(w_l, .25_PRC, (data_down - data)**2, dimX, num_fields, cweno_epsilon)
        call five_moment_mhd_compute_weight_alpha(w_r, .25_PRC, (data - data_up)**2, dimX, num_fields, cweno_epsilon)
        call five_moment_mhd_compute_weight_alpha(w_c, .5_PRC, (data_up - data_down)**2 * .25_PRC + (c1)**2 *13._PRC/3._PRC, &
                dimX, num_fields, cweno_epsilon)
        w = w_l + w_r + w_c
        w_l = w_l / w
        w_r = w_r / w
        w_c = w_c / w
        c1 = w_c * c1
        a = data + c1/6._PRC
        b = ( ( (.5_PRC*w_c + w_r) * data_up + (-w_r + w_l) * data - (.5_PRC*w_c + w_l) * data_down ) ) / 2._PRC

        reco(:,:,:, :, dir, 1) = &
                a(-shift(dir,1):dimX(1)+1-shift(dir,1), -shift(dir,2):dimX(2)+1-shift(dir,2), &
                -  shift(dir,3):dimX(3)+1-shift(dir,3), :) + &
                b(-shift(dir,1):dimX(1)+1-shift(dir,1), -shift(dir,2):dimX(2)+1-shift(dir,2), &
                -  shift(dir,3):dimX(3)+1-shift(dir,3), :)
        reco(:,:,:, :, dir, 2) = &
                a(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, :) - b(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, :)
        ! reco_minus[pos_up] = a + b
        ! reco_plus[pos]     = a - b
        ! 1: _minus, 2: _plus
    end do

    deallocate(w_c)
    deallocate(w_l)
    deallocate(w_r)
    deallocate(w)
    deallocate(c1)
    deallocate(a)
    deallocate(b)

end subroutine five_moment_mhd_cweno_reconstruct


subroutine five_moment_mhd_compute_weight_alpha(alpha, C, I, dimX, num_fields, cweno_epsilon)
    use definitions
    implicit none

    integer, intent(in) :: dimX(1:3), num_fields
    real(kind=PRC), intent(in) :: C, cweno_epsilon
    real(kind=PRC), intent(in) :: I(-1:dimX(1)+1, -1:dimX(2)+1, -1:dimX(3)+1, num_fields)
    real(kind=PRC), intent(out) :: alpha(-1:dimX(1)+1, -1:dimX(2)+1, -1:dimX(3)+1, num_fields)

    alpha = C / (cweno_epsilon + I)**2

end subroutine five_moment_mhd_compute_weight_alpha


function five_moment_mhd_max_dt(fluid, cfl, cweno_limiter, dimX, BD, dx) result(dt)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:5) :: fluid
    real(kind=PRC), intent(in) :: dx(3), cfl, cweno_limiter
    real(kind=PRC) :: dt

    real(kind=PRC), allocatable, dimension(:,:,:) :: n_inv, nunu
    real(kind=PRC) :: v_max(3)
    integer :: dir

    allocate(n_inv(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)))
    allocate(nunu(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)))

    n_inv = 1._PRC/max(fluid(:,:,:,1), cweno_limiter)
    nunu = fluid(:,:,:,2)**2 + fluid(:,:,:,3)**2 + fluid(:,:,:,4)**2
    do dir = 1,3
        v_max(dir) = maxval( abs(fluid(:,:,:,dir+1)) * n_inv &
                + sqrt(5._PRC/3._PRC * max(fluid(:,:,:,5) - nunu * n_inv, 0._PRC) * n_inv ) )
    end do

    ! stability
    dt = minval( dx / v_max ) * cfl

    ! resolve gyrofrequency (three steps per gyration)
    ! TODO check
    ! dt = min(dt, 0.3/maxval(sqrt(B(:,:,:,1)**2+B(:,:,:,2)**2+B(:,:,:,3)**2))*m)

    deallocate(n_inv)
    deallocate(nunu)

end function five_moment_mhd_max_dt

