
! This file is part of the muphyII multiphysics plasma simulation project.
! 
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at https://mozilla.org/MPL/2.0/.


! Centrally weighted essentially non-oscillating (CWENO) method

! fluid(:,:,:,1) = n
! fluid(:,:,:,2:4) = n u
! fluid(:,:,:,5:10) = Eps_ij = \int v_i v_j f dv (energy density tensor)
! eps = tr(Eps_ij) = 3 n k T / m + n u \dot u (scalar energy density)
! P_ij = m \int (v_i - u_i) (v_j - u_j) f dv (pressure tensor)
!      = m Eps_ij - m n u_i u_j
! p = tr(p_ij)/N (scalar pressure; N dimensionality)

subroutine ten_moment_calc_single_step_rk3(fluid_interim, total_flux, source, E, B, q_m, &
                  cweno_limiter, cweno_epsilon, dimensionality_x, dimX, BD, dx)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), dimensionality_x
    real(kind=PRC), intent(in) :: cweno_limiter, cweno_epsilon, q_m, dx(3)

    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:10), intent(inout) :: fluid_interim
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:3 ), intent(in) :: E, B
    real(kind=PRC), dimension(0:dimX(1),0:dimX(2),0:dimX(3),1:10), intent(inout) :: total_flux, source

    real(kind=PRC), allocatable, dimension(:,:,:,:,:,:) :: num_flux, reco
    real(kind=PRC), allocatable, dimension(:,:,:,:,:) :: flux
    real(kind=PRC), allocatable, dimension(:,:,:,:) :: speed

    allocate(num_flux(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:10, 1:3, 1:2))
    allocate(reco(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:10, 1:3, 1:2))
    allocate(flux(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:10, 1:3))
    allocate(speed(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:3))

    !$acc enter data create(num_flux, reco, flux, speed)

    call ten_moment_cweno_reconstruct(dimX, dimensionality_x, BD, fluid_interim, reco, 10, cweno_epsilon)

    call ten_moment_calc_flux(dimX, dimensionality_x, reco, num_flux, speed, cweno_limiter)

    call ten_moment_apply_flux(dimX, dx, dimensionality_x, reco, flux, num_flux, speed, total_flux)

    call ten_moment_calc_source(dimX, BD, fluid_interim, E, B, source, q_m)

    !$acc exit data delete(num_flux, reco, flux, speed)

    deallocate(num_flux)
    deallocate(reco)
    deallocate(flux)
    deallocate(speed)

end subroutine ten_moment_calc_single_step_rk3


subroutine ten_moment_deltaf_calc_single_step_rk3(fluid_interim, total_flux, source, Q_raw, E, B, q_m, &
                  cweno_limiter, cweno_epsilon, dimensionality_x, dimX, BD, dx)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), dimensionality_x
    real(kind=PRC), intent(in) :: cweno_limiter, cweno_epsilon, q_m, dx(3)

    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:10), intent(inout) :: fluid_interim
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:10), intent(in) :: Q_raw
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:3 ), intent(in) :: E, B
    real(kind=PRC), dimension(0:dimX(1),0:dimX(2),0:dimX(3),1:10), intent(inout) :: total_flux, source

    real(kind=PRC), allocatable, dimension(:,:,:,:,:,:) :: num_flux, reco, reco_Q
    real(kind=PRC), allocatable, dimension(:,:,:,:,:) :: flux
    real(kind=PRC), allocatable, dimension(:,:,:,:) :: speed

    allocate(num_flux(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:10, 1:3, 1:2))
    allocate(reco(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:10, 1:3, 1:2))
    allocate(reco_Q(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:10, 1:3, 1:2))
    allocate(flux(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:10, 1:3))
    allocate(speed(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:3))

    !$acc enter data create(num_flux, reco, reco_Q, flux, speed)

    call ten_moment_cweno_reconstruct(dimX, dimensionality_x, BD, fluid_interim, reco, 10, cweno_epsilon)

    call ten_moment_cweno_reconstruct(dimX, dimensionality_x, BD, Q_raw, reco_Q, 10, cweno_epsilon)

    call ten_moment_deltaf_calc_flux(dimX, dimensionality_x, reco, reco_Q, num_flux, speed, cweno_limiter)

    call ten_moment_apply_flux(dimX, dx, dimensionality_x, reco, flux, num_flux, speed, total_flux)

    call ten_moment_calc_source(dimX, BD, fluid_interim, E, B, source, q_m)

    !$acc exit data delete(num_flux, reco, reco_Q, flux, speed)

    deallocate(num_flux)
    deallocate(reco)
    deallocate(reco_Q)
    deallocate(flux)
    deallocate(speed)

end subroutine ten_moment_deltaf_calc_single_step_rk3


subroutine ten_moment_apply_single_step_rk3(fluid, fluid_interim, total_flux, source, alpha0, alpha1, beta, &
                  cweno_limiter, dimensionality_v, dimX, BD, dt)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), dimensionality_v
    real(kind=PRC) :: alpha0, alpha1, beta, cweno_limiter, dt
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:10), &
            intent(in) :: fluid
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:10), &
            intent(inout) :: fluid_interim
    real(kind=PRC), dimension(0:dimX(1),0:dimX(2),0:dimX(3),1:10), intent(in) :: source, total_flux

    integer :: fix, x, y, z

    ! ensure positiveness of temperature
    ! TODO check implications
    ! TODO treatment of non-diagonal elements (positive-semi-definiteness of the tensor) needed?

    !$acc data present(fluid, fluid_interim, total_flux, source)

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
                fluid_interim(x,y,z,1) = max(fluid_interim(x,y,z,1), cweno_limiter)
                fluid_interim(x,y,z,5) = max(fluid_interim(x,y,z,5), fluid_interim(x,y,z,2)**2/fluid_interim(x,y,z,1) + cweno_limiter)
                if (dimensionality_v > 1) then
                    fluid_interim(x,y,z,8) = max(fluid_interim(x,y,z,8), fluid_interim(x,y,z,3)**2/fluid_interim(x,y,z,1) + cweno_limiter)
                end if
                if (dimensionality_v > 2) then
                    fluid_interim(x,y,z,10) = max(fluid_interim(x,y,z,10), fluid_interim(x,y,z,4)**2/fluid_interim(x,y,z,1) + cweno_limiter)
                end if
            end do
        end do
    end do
    !$acc end parallel

    ! runge-kutta-3 step (Shu-Osher)
    !$acc parallel
#ifdef ACC_2D_OPT
    !$acc loop gang collapse(3)
    do fix = 1,10
        do z = 0,dimX(3)
            do y = 0,dimX(2)
                !$acc loop vector
                do x = 0,dimX(1)
#else
    !$acc loop gang collapse(2)
    do fix = 1,10
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

end subroutine ten_moment_apply_single_step_rk3


subroutine ten_moment_calc_flux(dimX, dimensionality_x, reco, flux, speed, cweno_limiter)
    use definitions
    implicit none

    integer,  intent(in) :: dimX(3), dimensionality_x
    real(kind=PRC), intent(in) :: cweno_limiter
    real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:10, 1:3, 1:2), intent(in)    :: reco
    real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:10, 1:3, 1:2), intent(inout) :: flux
    real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1,       1:3),      intent(inout) :: speed

    integer :: s, d, x, y, z
    real(kind=PRC) :: r_n_inv

    !$acc data present(reco, flux, speed)

    do s = 1,2 ! side
        do d = 1,dimensionality_x ! direction
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
                        ! what follows represents get_flux() in muphy1
                        r_n_inv = 1._PRC/max(reco(x,y,z,1,d,s), cweno_limiter)

                        if (d == 1) then
                            flux(x,y,z,1,d,s)  = reco(x,y,z,2,d,s) ! r_nux
                            flux(x,y,z,2,d,s)  = reco(x,y,z,5,d,s) ! r_Epsxx
                            flux(x,y,z,3,d,s)  = reco(x,y,z,6,d,s) ! r_Epsxy
                            flux(x,y,z,4,d,s)  = reco(x,y,z,7,d,s) ! r_Epsxz
                            ! 3._PRC*r_nux*r_Epsxx*r_n_inv - 2._PRC*r_nux*r_nux*r_nux*r_n_inv*r_n_inv
                            flux(x,y,z,5,d,s)  = 3._PRC*reco(x,y,z,2,d,s)*reco(x,y,z,5,d,s)*r_n_inv &
                                                 - 2._PRC*reco(x,y,z,2,d,s)*reco(x,y,z,2,d,s)*reco(x,y,z,2,d,s)*r_n_inv*r_n_inv ! xx
                            ! (r_nuy*r_Epsxx + 2._PRC*r_nux*r_Epsxy)*r_n_inv - 2._PRC*r_nux*r_nux*r_nuy*r_n_inv*r_n_inv
                            flux(x,y,z,6,d,s)  = (reco(x,y,z,3,d,s)*reco(x,y,z,5,d,s) + 2._PRC*reco(x,y,z,2,d,s)*reco(x,y,z,6,d,s))*r_n_inv &
                                                 - 2._PRC*reco(x,y,z,2,d,s)*reco(x,y,z,2,d,s)*reco(x,y,z,3,d,s)*r_n_inv*r_n_inv ! xy
                            ! (r_nuz*r_Epsxx + 2._PRC*r_nux*r_Epsxz)*r_n_inv - 2._PRC*r_nux*r_nux*r_nuz*r_n_inv*r_n_inv
                            flux(x,y,z,7,d,s)  = (reco(x,y,z,4,d,s)*reco(x,y,z,5,d,s) + 2._PRC*reco(x,y,z,2,d,s)*reco(x,y,z,7,d,s))*r_n_inv &
                                                 - 2._PRC*reco(x,y,z,2,d,s)*reco(x,y,z,2,d,s)*reco(x,y,z,4,d,s)*r_n_inv*r_n_inv ! xz
                            ! (r_nux*r_Epsyy + 2._PRC*r_nuy*r_Epsxy)*r_n_inv - 2._PRC*r_nux*r_nuy*r_nuy*r_n_inv*r_n_inv
                            flux(x,y,z,8,d,s)  = (reco(x,y,z,2,d,s)*reco(x,y,z,8,d,s) + 2._PRC*reco(x,y,z,3,d,s)*reco(x,y,z,6,d,s))*r_n_inv &
                                                 - 2._PRC*reco(x,y,z,2,d,s)*reco(x,y,z,3,d,s)*reco(x,y,z,3,d,s)*r_n_inv*r_n_inv ! yy
                            ! (r_nux*r_Epsyz + r_nuy*r_Epsxz + r_nuz*r_Epsxy)*r_n_inv - 2._PRC*r_nux*r_nuy*r_nuz*r_n_inv*r_n_inv
                            flux(x,y,z,9,d,s)  = (reco(x,y,z,2,d,s)*reco(x,y,z,9,d,s) + reco(x,y,z,3,d,s)*reco(x,y,z,7,d,s) + reco(x,y,z,4,d,s)*reco(x,y,z,6,d,s))*r_n_inv &
                                                 - 2._PRC*reco(x,y,z,2,d,s)*reco(x,y,z,3,d,s)*reco(x,y,z,4,d,s)*r_n_inv*r_n_inv ! yz
                            ! (r_nux*r_Epszz + 2._PRC*r_nuz*r_Epsxz)*r_n_inv - 2._PRC*r_nux*r_nuz*r_nuz*r_n_inv*r_n_inv
                            flux(x,y,z,10,d,s) = (reco(x,y,z,2,d,s)*reco(x,y,z,10,d,s) + 2._PRC*reco(x,y,z,4,d,s)*reco(x,y,z,7,d,s))*r_n_inv &
                                                 - 2._PRC*reco(x,y,z,2,d,s)*reco(x,y,z,4,d,s)*reco(x,y,z,4,d,s)*r_n_inv*r_n_inv ! zz
                            ! abs(r_nux) * r_n_inv + sqrt(3._PRC * max(r_Epsxx - r_nux*r_nux*r_n_inv, 0._PRC) * r_n_inv)
                            speed(x,y,z,d) = abs(reco(x,y,z,2,d,s)) * r_n_inv + &
                                             sqrt(3._PRC * max(reco(x,y,z,5,d,s) - reco(x,y,z,2,d,s)*reco(x,y,z,2,d,s)*r_n_inv, 0._PRC) * r_n_inv)
                        else if (d == 2) then
                            flux(x,y,z,1,d,s)  = reco(x,y,z,3,d,s) ! r_nuy
                            flux(x,y,z,2,d,s)  = reco(x,y,z,6,d,s) ! r_Epsxy
                            flux(x,y,z,3,d,s)  = reco(x,y,z,8,d,s) ! r_Epsyy
                            flux(x,y,z,4,d,s)  = reco(x,y,z,9,d,s) ! r_Epsyz
                            ! (r_nuy*r_Epsxx + 2._PRC*r_nux*r_Epsxy)*r_n_inv - 2._PRC*r_nuy*r_nux*r_nux*r_n_inv*r_n_inv
                            flux(x,y,z,5,d,s)  = (reco(x,y,z,3,d,s)*reco(x,y,z,5,d,s) + 2._PRC*reco(x,y,z,2,d,s)*reco(x,y,z,6,d,s))*r_n_inv &
                                                 - 2._PRC*reco(x,y,z,3,d,s)*reco(x,y,z,2,d,s)*reco(x,y,z,2,d,s)*r_n_inv*r_n_inv ! xx

                            ! (r_nux*r_Epsyy + 2._PRC*r_nuy*r_Epsxy)*r_n_inv - 2._PRC*r_nuy*r_nux*r_nuy*r_n_inv*r_n_inv
                            flux(x,y,z,6,d,s)  = (reco(x,y,z,2,d,s)*reco(x,y,z,8,d,s) + 2._PRC*reco(x,y,z,3,d,s)*reco(x,y,z,6,d,s))*r_n_inv &
                                                 - 2._PRC*reco(x,y,z,3,d,s)*reco(x,y,z,2,d,s)*reco(x,y,z,3,d,s)*r_n_inv*r_n_inv ! xy
                            ! (r_nux*r_Epsyz + r_nuy*r_Epsxz + r_nuz*r_Epsxy)*r_n_inv - 2._PRC*r_nuy*r_nux*r_nuz*r_n_inv*r_n_inv
                            flux(x,y,z,7,d,s)  = (reco(x,y,z,2,d,s)*reco(x,y,z,9,d,s) + reco(x,y,z,3,d,s)*reco(x,y,z,7,d,s) + reco(x,y,z,4,d,s)*reco(x,y,z,6,d,s))*r_n_inv &
                                                 - 2._PRC*reco(x,y,z,3,d,s)*reco(x,y,z,2,d,s)*reco(x,y,z,4,d,s)*r_n_inv*r_n_inv ! xz
                            ! 3._PRC*r_nuy*r_Epsyy*r_n_inv - 2._PRC*r_nuy*r_nuy*r_nuy*r_n_inv*r_n_inv
                            flux(x,y,z,8,d,s)  = 3._PRC*reco(x,y,z,3,d,s)*reco(x,y,z,8,d,s)*r_n_inv &
                                                 - 2._PRC*reco(x,y,z,3,d,s)*reco(x,y,z,3,d,s)*reco(x,y,z,3,d,s)*r_n_inv*r_n_inv ! yy
                            ! (r_nuz*r_Epsyy + 2._PRC*r_nuy*r_Epsyz)*r_n_inv - 2._PRC*r_nuy*r_nuy*r_nuz*r_n_inv*r_n_inv
                            flux(x,y,z,9,d,s)  = (reco(x,y,z,4,d,s)*reco(x,y,z,8,d,s) + 2._PRC*reco(x,y,z,3,d,s)*reco(x,y,z,9,d,s))*r_n_inv &
                                                 - 2._PRC*reco(x,y,z,3,d,s)*reco(x,y,z,3,d,s)*reco(x,y,z,4,d,s)*r_n_inv*r_n_inv ! yz
                            ! (r_nuy*r_Epszz + 2._PRC*r_nuz*r_Epsyz)*r_n_inv - 2._PRC*r_nuy*r_nuz*r_nuz*r_n_inv*r_n_inv
                            flux(x,y,z,10,d,s) = (reco(x,y,z,3,d,s)*reco(x,y,z,10,d,s) + 2._PRC*reco(x,y,z,4,d,s)*reco(x,y,z,9,d,s))*r_n_inv &
                                                 - 2._PRC*reco(x,y,z,3,d,s)*reco(x,y,z,4,d,s)*reco(x,y,z,4,d,s)*r_n_inv*r_n_inv ! zz
                            ! abs(r_nuy) * r_n_inv + sqrt(3._PRC * max(r_Epsyy - r_nuy*r_nuy*r_n_inv, 0._PRC) * r_n_inv)
                            speed(x,y,z,d) = abs(reco(x,y,z,3,d,s)) * r_n_inv + &
                                             sqrt(3._PRC * max(reco(x,y,z,8,d,s) - reco(x,y,z,3,d,s)*reco(x,y,z,3,d,s)*r_n_inv, 0._PRC) * r_n_inv)
                        else if (d == 3) then
                            flux(x,y,z,1,d,s)  = reco(x,y,z,4,d,s)  ! r_nuz
                            flux(x,y,z,2,d,s)  = reco(x,y,z,7,d,s)  ! r_Epsxz
                            flux(x,y,z,3,d,s)  = reco(x,y,z,9,d,s)  ! r_Epsyz
                            flux(x,y,z,4,d,s)  = reco(x,y,z,10,d,s) ! r_Epszz
                            ! (r_nuz*r_Epsxx + 2._PRC*r_nux*r_Epsxz)*r_n_inv - 2._PRC*r_nuz*r_nux*r_nux*r_n_inv*r_n_inv
                            flux(x,y,z,5,d,s)  = (reco(x,y,z,4,d,s)*reco(x,y,z,5,d,s) + 2._PRC*reco(x,y,z,2,d,s)*reco(x,y,z,7,d,s))*r_n_inv &
                                                 - 2._PRC*reco(x,y,z,4,d,s)*reco(x,y,z,2,d,s)*reco(x,y,z,2,d,s)*r_n_inv*r_n_inv ! xx
                            ! (r_nux*r_Epsyz + r_nuy*r_Epsxz + r_nuz*r_Epsxy)*r_n_inv - 2._PRC*r_nuz*r_nux*r_nuy*r_n_inv*r_n_inv
                            flux(x,y,z,6,d,s)  = (reco(x,y,z,2,d,s)*reco(x,y,z,9,d,s) + reco(x,y,z,3,d,s)*reco(x,y,z,7,d,s) + reco(x,y,z,4,d,s)*reco(x,y,z,6,d,s))*r_n_inv &
                                                 - 2._PRC*reco(x,y,z,4,d,s)*reco(x,y,z,2,d,s)*reco(x,y,z,3,d,s)*r_n_inv*r_n_inv ! xy
                            ! (r_nux*r_Epszz + 2._PRC*r_nuz*r_Epsxz)*r_n_inv - 2._PRC*r_nuz*r_nux*r_nuz*r_n_inv*r_n_inv
                            flux(x,y,z,7,d,s)  = (reco(x,y,z,2,d,s)*reco(x,y,z,10,d,s) + 2._PRC*reco(x,y,z,4,d,s)*reco(x,y,z,7,d,s))*r_n_inv &
                                                 - 2._PRC*reco(x,y,z,4,d,s)*reco(x,y,z,2,d,s)*reco(x,y,z,4,d,s)*r_n_inv*r_n_inv ! xz
                            ! (r_nuz*r_Epsyy + 2._PRC*r_nuy*r_Epsyz)*r_n_inv - 2._PRC*r_nuz*r_nuy*r_nuy*r_n_inv*r_n_inv
                            flux(x,y,z,8,d,s)  = (reco(x,y,z,4,d,s)*reco(x,y,z,8,d,s) + 2._PRC*reco(x,y,z,3,d,s)*reco(x,y,z,9,d,s))*r_n_inv &
                                                 - 2._PRC*reco(x,y,z,4,d,s)*reco(x,y,z,3,d,s)*reco(x,y,z,3,d,s)*r_n_inv*r_n_inv ! yy
                            ! (r_nuy*r_Epszz + 2._PRC*r_nuz*r_Epsyz)*r_n_inv - 2._PRC*r_nuz*r_nuy*r_nuz*r_n_inv*r_n_inv
                            flux(x,y,z,9,d,s)  = (reco(x,y,z,3,d,s)*reco(x,y,z,10,d,s) + 2._PRC*reco(x,y,z,4,d,s)*reco(x,y,z,9,d,s))*r_n_inv &
                                                 - 2._PRC*reco(x,y,z,4,d,s)*reco(x,y,z,3,d,s)*reco(x,y,z,4,d,s)*r_n_inv*r_n_inv ! yz
                            ! 3._PRC*r_nuz*r_Epszz*r_n_inv - 2._PRC*r_nuz*r_nuz*r_nuz*r_n_inv*r_n_inv
                            flux(x,y,z,10,d,s) = 3._PRC*reco(x,y,z,4,d,s)*reco(x,y,z,10,d,s)*r_n_inv &
                                                 - 2._PRC*reco(x,y,z,4,d,s)*reco(x,y,z,4,d,s)*reco(x,y,z,4,d,s)*r_n_inv*r_n_inv ! zz
                            ! max( speed(:,:,:,d), abs(r_nuz) * r_n_inv + sqrt(3._PRC * max(r_Epszz - r_nuz*r_nuz*r_n_inv, 0._PRC) * r_n_inv) )
                            speed(x,y,z,d) = abs(reco(x,y,z,4,d,s)) * r_n_inv + &
                                             sqrt(3._PRC * max(reco(x,y,z,10,d,s) - reco(x,y,z,4,d,s)*reco(x,y,z,4,d,s)*r_n_inv, 0._PRC) * r_n_inv)
                        end if
                    end do
                end do
            end do
            !$acc end parallel
        end do
    end do

    !$acc end data

end subroutine ten_moment_calc_flux


subroutine ten_moment_deltaf_calc_flux(dimX, dimensionality_x, reco, reco_Q, flux, speed, cweno_limiter)
    use definitions
    implicit none

    integer,  intent(in) :: dimX(3), dimensionality_x
    real(kind=PRC), intent(in) :: cweno_limiter
    real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:10, 1:3, 1:2), intent(in)    :: reco, reco_Q
    real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:10, 1:3, 1:2), intent(inout) :: flux
    real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1,       1:3),      intent(inout) :: speed

    integer :: s, d, x, y, z
    real(kind=PRC) :: r_n_inv

    !$acc data present(reco, reco_Q, flux, speed)

    do s = 1,2 ! side
        do d = 1,dimensionality_x ! direction
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
                        r_n_inv = 1._PRC/max(reco(x,y,z,1,d,s), cweno_limiter)

                        if (d == 1) then
                            flux(x,y,z,1, d,s) = reco(x,y,z,2,d,s)    ! r_nux
                            flux(x,y,z,2, d,s) = reco(x,y,z,5,d,s)    ! r_Epsxx
                            flux(x,y,z,3, d,s) = reco(x,y,z,6,d,s)    ! r_Epsxy
                            flux(x,y,z,4, d,s) = reco(x,y,z,7,d,s)    ! r_Epsxz
                            flux(x,y,z,5, d,s) = reco_Q(x,y,z,1,d,s)  ! r_Q_raw_xxx
                            flux(x,y,z,6, d,s) = reco_Q(x,y,z,2,d,s)  ! r_Q_raw_xxy
                            flux(x,y,z,7, d,s) = reco_Q(x,y,z,3,d,s)  ! r_Q_raw_xxz
                            flux(x,y,z,8, d,s) = reco_Q(x,y,z,4,d,s)  ! r_Q_raw_xyy
                            flux(x,y,z,9, d,s) = reco_Q(x,y,z,5,d,s)  ! r_Q_raw_xyz
                            flux(x,y,z,10,d,s) = reco_Q(x,y,z,6,d,s)  ! r_Q_raw_xzz
                            ! abs(r_nux) * r_n_inv + sqrt(3._PRC * max(r_Epsxx - r_nux*r_nux*r_n_inv, 0._PRC) * r_n_inv)
                            speed(x,y,z,d) = abs(reco(x,y,z,2,d,s)) * r_n_inv + &
                                             sqrt(3._PRC * max(reco(x,y,z,5,d,s) - reco(x,y,z,2,d,s)*reco(x,y,z,2,d,s)*r_n_inv, 0._PRC) * r_n_inv)
                        else if (d == 2) then
                            flux(x,y,z,1, d,s) = reco(x,y,z,3,d,s)    ! r_nuy
                            flux(x,y,z,2, d,s) = reco(x,y,z,6,d,s)    ! r_Epsxy
                            flux(x,y,z,3, d,s) = reco(x,y,z,8,d,s)    ! r_Epsyy
                            flux(x,y,z,4, d,s) = reco(x,y,z,9,d,s)    ! r_Epsyz
                            flux(x,y,z,5, d,s) = reco_Q(x,y,z,2,d,s)  ! r_Q_raw_xxy
                            flux(x,y,z,6, d,s) = reco_Q(x,y,z,4,d,s)  ! r_Q_raw_xyy
                            flux(x,y,z,7, d,s) = reco_Q(x,y,z,5,d,s)  ! r_Q_raw_xyz
                            flux(x,y,z,8, d,s) = reco_Q(x,y,z,7,d,s)  ! r_Q_raw_yyy
                            flux(x,y,z,9, d,s) = reco_Q(x,y,z,8,d,s)  ! r_Q_raw_yyz
                            flux(x,y,z,10,d,s) = reco_Q(x,y,z,9,d,s)  ! r_Q_raw_yzz
                            ! abs(r_nuy) * r_n_inv + sqrt(3._PRC * max(r_Epsyy - r_nuy*r_nuy*r_n_inv, 0._PRC) * r_n_inv)
                            speed(x,y,z,d) = abs(reco(x,y,z,3,d,s)) * r_n_inv + &
                                             sqrt(3._PRC * max(reco(x,y,z,8,d,s) - reco(x,y,z,3,d,s)*reco(x,y,z,3,d,s)*r_n_inv, 0._PRC) * r_n_inv)
                        else if (d == 3) then
                            flux(x,y,z,1, d,s) = reco(x,y,z,4,d,s)    ! r_nuz
                            flux(x,y,z,2, d,s) = reco(x,y,z,7,d,s)    ! r_Epsxz
                            flux(x,y,z,3, d,s) = reco(x,y,z,9,d,s)    ! r_Epsyz
                            flux(x,y,z,4, d,s) = reco(x,y,z,10,d,s)   ! r_Epszz
                            flux(x,y,z,5, d,s) = reco_Q(x,y,z,3,d,s)  ! r_Q_raw_xxz
                            flux(x,y,z,6, d,s) = reco_Q(x,y,z,5,d,s)  ! r_Q_raw_xyz
                            flux(x,y,z,7, d,s) = reco_Q(x,y,z,6,d,s)  ! r_Q_raw_xzz
                            flux(x,y,z,8, d,s) = reco_Q(x,y,z,8,d,s)  ! r_Q_raw_yyz
                            flux(x,y,z,9, d,s) = reco_Q(x,y,z,9,d,s)  ! r_Q_raw_yzz
                            flux(x,y,z,10,d,s) = reco_Q(x,y,z,10,d,s) ! r_Q_raw_zzz
                            ! max( speed(:,:,:,d), abs(r_nuz) * r_n_inv + sqrt(3._PRC * max(r_Epszz - r_nuz*r_nuz*r_n_inv, 0._PRC) * r_n_inv) )
                            speed(x,y,z,d) = abs(reco(x,y,z,4,d,s)) * r_n_inv + &
                                             sqrt(3._PRC * max(reco(x,y,z,10,d,s) - reco(x,y,z,4,d,s)*reco(x,y,z,4,d,s)*r_n_inv, 0._PRC) * r_n_inv)
                        end if
                    end do
                end do
            end do
            !$acc end parallel
        end do
    end do

    !$acc end data

end subroutine ten_moment_deltaf_calc_flux


subroutine ten_moment_apply_flux(dimX, dx, dimensionality_x, reco, flux, num_flux, speed, total_flux)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimensionality_x
    real(kind=PRC), intent(in) :: dx(3)
    real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:10, 1:3, 1:2), intent(in) :: reco
    real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:10, 1:3),      intent(inout) :: flux
    real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:10, 1:3, 1:2), intent(in) :: num_flux
    real(kind=PRC), dimension( 0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1,       1:3),      intent(in) :: speed
    real(kind=PRC), dimension( 0:dimX(1)  , 0:dimX(2)  , 0:dimX(3),   1:10),           intent(inout) :: total_flux

    integer :: dir, fix, x, y, z

    !$acc data present(reco, flux, num_flux, speed, total_flux)

    do dir = 1,dimensionality_x
        !$acc parallel
#ifdef ACC_2D_OPT
        !$acc loop gang collapse(3)
        do fix = 1,10 ! field index
            do z = 0,dimX(3)+1
                do y = 0,dimX(2)+1
                    !$acc loop vector
                    do x = 0,dimX(1)+1
#else
        !$acc loop gang collapse(2)
        do fix = 1,10 ! field index
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
    do fix = 1,10
        do z = 0,dimX(3)
            do y = 0,dimX(2)
                !$acc loop vector
                do x = 0,dimX(1)
#else
    !$acc loop gang collapse(2)
    do fix = 1,10
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
        do fix = 1,10
            do z = 0,dimX(3)
                do y = 0,dimX(2)
                    !$acc loop vector
                    do x = 0,dimX(1)
#else
        !$acc loop gang collapse(2)
        do fix = 1,10
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
        do fix = 1,10
            do z = 0,dimX(3)
                do y = 0,dimX(2)
                    !$acc loop vector
                    do x = 0,dimX(1)
#else
        !$acc loop gang collapse(2)
        do fix = 1,10
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

    !$acc end data

end subroutine ten_moment_apply_flux


subroutine ten_moment_calc_source(dimX, BD, fluid_interim, E, B, source, q_m)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:10), intent(in) :: fluid_interim
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:3 ), intent(in) :: E, B
    real(kind=PRC), dimension(0:dimX(1),0:dimX(2),0:dimX(3),1:10), intent(inout) :: source
    real(kind=PRC), intent(in) :: q_m

    integer :: x, y, z

    !$acc data present(fluid_interim, E, B, source)

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
                ! 2 * q_m * ( Ex * nux + Bz * Epsxy - By * Epsxz )
                source(x,y,z,5)  = 2._PRC*q_m * (E(x,y,z,1)*fluid_interim(x,y,z,2) + B(x,y,z,3)*fluid_interim(x,y,z,6) - B(x,y,z,2)*fluid_interim(x,y,z,7) ) ! xx
                ! 2 * q_m * ( Ey * nuy + Bx * Epsyz - Bz * Epsxy )
                source(x,y,z,8)  = 2._PRC*q_m * (E(x,y,z,2)*fluid_interim(x,y,z,3) + B(x,y,z,1)*fluid_interim(x,y,z,9) - B(x,y,z,3)*fluid_interim(x,y,z,6) ) ! yy
                ! 2 * q_m * ( Ez * nuz + By * Epsxz - Bx * Epsyz )
                source(x,y,z,10) = 2._PRC*q_m * (E(x,y,z,3)*fluid_interim(x,y,z,4) + B(x,y,z,2)*fluid_interim(x,y,z,7) - B(x,y,z,1)*fluid_interim(x,y,z,9) ) ! zz
                ! q_m * ( Ex * nuy + Ey * nux + Bz * (Epsyy - Epsxx) + Bx * Epsxz - By * Epsyz )
                source(x,y,z,6)  = q_m * (E(x,y,z,1)*fluid_interim(x,y,z,3) + E(x,y,z,2)*fluid_interim(x,y,z,2) + &
                                          B(x,y,z,3)*(fluid_interim(x,y,z,8) - fluid_interim(x,y,z,5)) + &
                                          B(x,y,z,1)*fluid_interim(x,y,z,7) - B(x,y,z,2)*fluid_interim(x,y,z,9)) ! xy
                ! q_m * ( Ex * nuz + Ez * nux + By * (Epsxx - Epszz) + Bz * Epsyz - Bx * Epsxy )
                source(x,y,z,7)  = q_m * (E(x,y,z,1)*fluid_interim(x,y,z,4) + E(x,y,z,3)*fluid_interim(x,y,z,2) + &
                                          B(x,y,z,2)*(fluid_interim(x,y,z,5) - fluid_interim(x,y,z,10)) + &
                                          B(x,y,z,3)*fluid_interim(x,y,z,9) - B(x,y,z,1)*fluid_interim(x,y,z,6)) ! xz
                ! q_m * ( Ey * nuz + Ez * nuy + Bx * (Epszz - Epsyy) + By * Epsxy - Bz * Epsxz )
                source(x,y,z,9)  = q_m * (E(x,y,z,2)*fluid_interim(x,y,z,4) + E(x,y,z,3)*fluid_interim(x,y,z,3) + &
                                          B(x,y,z,1)*(fluid_interim(x,y,z,10) - fluid_interim(x,y,z,8)) + &
                                          B(x,y,z,2)*fluid_interim(x,y,z,6) - B(x,y,z,3)*fluid_interim(x,y,z,7)) ! yz
            end do
        end do
    end do
    !$acc end parallel

    !$acc end data

end subroutine ten_moment_calc_source


subroutine ten_moment_calc_un_from_flux(un, fluid, cweno_limiter, cweno_epsilon, dimX, BD)
    use definitions
    implicit none

    integer,  intent(in) :: dimX(1:3), BD(3)
    real(kind=PRC), intent(in) :: cweno_limiter, cweno_epsilon
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:10), intent(in) :: fluid
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:3), intent(inout) :: un

    real(kind=PRC), allocatable, dimension(:,:,:,:,:,:) :: num_flux
    real(kind=PRC), allocatable, dimension(:,:,:,:,:,:) :: reco
    real(kind=PRC), allocatable, dimension(:,:,:,:) :: speed
    integer :: dir, x, y, z

    allocate(num_flux(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:10, 1:3, 1:2))
    allocate(reco(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:10, 1:3, 1:2))
    allocate(speed(0:dimX(1)+1, 0:dimX(2)+1, 0:dimX(3)+1, 1:3))

    !$acc enter data create(num_flux, reco, speed)

    !$acc data present(num_flux, reco, speed, un, fluid)

    call ten_moment_cweno_reconstruct(dimX, 3, BD, fluid, reco, 10, cweno_epsilon)

    ! TODO there is some overhead due to reusing the function
    call ten_moment_calc_flux(dimX, 3, reco, num_flux, speed, cweno_limiter)

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

end subroutine ten_moment_calc_un_from_flux


function ten_moment_max_dt(fluid, cfl, cweno_limiter, dimX, BD, dx) result(dt)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),10) :: fluid
    real(kind=PRC), intent(in) :: dx(3), cfl, cweno_limiter

    real(kind=PRC), allocatable, dimension(:,:,:,:) :: v
    real(kind=PRC) :: n_inv, v_max
    integer :: x, y, z
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

                v(x,y,z,1) = abs(fluid(x,y,z,2)) * n_inv &
                     + sqrt(3._PRC * max(fluid(x,y,z,5) - fluid(x,y,z,2)*fluid(x,y,z,2)*n_inv, 0._PRC) * n_inv )
                v(x,y,z,2) = abs(fluid(x,y,z,3)) * n_inv &
                     + sqrt(3._PRC * max(fluid(x,y,z,8) - fluid(x,y,z,3)*fluid(x,y,z,3)*n_inv, 0._PRC) * n_inv )
                v(x,y,z,3) = abs(fluid(x,y,z,4)) * n_inv &
                     + sqrt(3._PRC * max(fluid(x,y,z,10) - fluid(x,y,z,4)*fluid(x,y,z,4)*n_inv, 0._PRC) * n_inv )
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

end function ten_moment_max_dt


subroutine ten_moment_cweno_reconstruct(dimX, dimensionality_x, BD, fields, reco, num_fields, cweno_epsilon)
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

    !$acc kernels
    reco = 0._PRC
    !$acc end kernels

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

end subroutine ten_moment_cweno_reconstruct


subroutine ten_moment_add_divheatflux_to_source(divQ_h, source, m, dimX)
    use definitions
    implicit none

    integer :: dimX(3)
    real(kind=PRC) :: m
    real(kind=PRC), intent(in), dimension(0:dimX(1),0:dimX(2),0:dimX(3),6) :: divQ_h
    real(kind=PRC), intent(inout), dimension(0:dimX(1),0:dimX(2),0:dimX(3),10) :: source

    integer :: i, x, y, z

    !$acc data present(divQ_h, source)

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            !$acc loop vector collapse(2)
            do x = 0,dimX(1)
                do i = 1,6
                    source(x,y,z,i+4) = source(x,y,z,i+4) - divQ_h(x,y,z,i)/m
                end do
            end do
        end do
    end do
    !$acc end parallel

    !$acc end data

end subroutine ten_moment_add_divheatflux_to_source


subroutine ten_moment_isotropization_heat_flux_closure(fluid_interim, source, k0, cweno_limiter, dimensionality_v, dimX, BD)
    use definitions
    implicit none

    integer :: dimensionality_v, dimX(3), BD(3)
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:10), &
            intent(in) :: fluid_interim
    real(kind=PRC), dimension(0:dimX(1),0:dimX(2),0:dimX(3),1:10), intent(inout) :: source
    real(kind=PRC) :: k0, cweno_limiter

    real(kind=PRC) :: tau, n_inv, P_m_xx, P_m_xy, P_m_xz, P_m_yy, P_m_yz, P_m_zz, tr_P_m
    integer :: x, y, z

    ! isotropization heat flux closure
    ! Wang, Hakim, Bhattacharjee, Germaschewski 2015
    ! https://arxiv.org/abs/1409.0262
    ! https://doi.org/10.1063/1.4906063

    ! typical choice of parameter::ten_moment_k0 is
    ! ten_moment_k0 = 1. / sqrt(m_s)


    !$acc data present(fluid_interim, source)

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
                n_inv = 1._PRC/max(fluid_interim(x,y,z,1), cweno_limiter)

                ! Epsxx - nux * nux * n_inv
                P_m_xx = fluid_interim(x,y,z,5) - fluid_interim(x,y,z,2)*fluid_interim(x,y,z,2)*n_inv
                ! Epsxy - nux * nuy * n_inv
                P_m_xy = fluid_interim(x,y,z,6) - fluid_interim(x,y,z,2)*fluid_interim(x,y,z,3)*n_inv
                ! Epsxz - nux * nuz * n_inv
                P_m_xz = fluid_interim(x,y,z,7) - fluid_interim(x,y,z,2)*fluid_interim(x,y,z,4)*n_inv
                ! Epsyy - nuy * nuy * n_inv
                P_m_yy = fluid_interim(x,y,z,8) - fluid_interim(x,y,z,3)*fluid_interim(x,y,z,3)*n_inv
                ! Epsyz - nuy * nuz * n_inv
                P_m_yz = fluid_interim(x,y,z,9) - fluid_interim(x,y,z,3)*fluid_interim(x,y,z,4)*n_inv
                ! Epszz - nuz * nuz * n_inv
                P_m_zz = fluid_interim(x,y,z,10) - fluid_interim(x,y,z,4)*fluid_interim(x,y,z,4)*n_inv

                if (dimensionality_v == 3) then
                    tr_P_m = (P_m_xx + P_m_yy + P_m_zz) / 3._PRC
                else if (dimensionality_v == 2) then
                    tr_P_m = (P_m_xx + P_m_yy) / 2._PRC
                else
                    tr_P_m = P_m_xx
                end if

                tau = k0 * sqrt(tr_P_m * n_inv)

                source(x,y,z,5)  = source(x,y,z,5)  - tau * (P_m_xx - tr_P_m)
                source(x,y,z,6)  = source(x,y,z,6)  - tau *  P_m_xy
                source(x,y,z,7)  = source(x,y,z,7)  - tau *  P_m_xz
                if (dimensionality_v > 1) then
                    source(x,y,z,8)  = source(x,y,z,8)  - tau * (P_m_yy - tr_P_m)
                end if
                source(x,y,z,9)  = source(x,y,z,9)  - tau *  P_m_yz
                if (dimensionality_v > 2) then
                    source(x,y,z,10) = source(x,y,z,10) - tau * (P_m_zz - tr_P_m)
                end if
            end do
        end do
    end do
    !$acc end parallel

    !$acc end data

end subroutine ten_moment_isotropization_heat_flux_closure


subroutine ten_moment_gradient_p_heat_flux_closure_cycle(fluid_interim, source, P_m_interim, k0, current_cycle, ncycles, &
                                                cweno_limiter, dimX, BD, dx, dimensionality_x, dimensionality_v, dt)
    use definitions
    implicit none

    integer :: dimX(3), BD(3), dimensionality_x, dimensionality_v, current_cycle, ncycles
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:10), &
            intent(in) :: fluid_interim
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:6), &
            intent(inout) :: P_m_interim
    real(kind=PRC), dimension(0:dimX(1),0:dimX(2),0:dimX(3),1:10), intent(inout) :: source
    real(kind=PRC) :: k0, cweno_limiter, dx(3), dt

    real(kind=PRC), allocatable, dimension(:,:,:) :: tau, n_inv, tr_P_m
    real(kind=PRC), allocatable, dimension(:,:,:,:) :: change

    integer :: i, x, y, z

    ! gradient heat flux closure
    ! https://arxiv.org/abs/2008.06440
    ! also see:
    ! https://doi.org/10.1017/S002237781800048X
    ! https://arxiv.org/abs/1801.07628
    ! for the k_0 values in the latter paper read
    ! k_0s = ... d_0i, not k_0s = ... d_0s

    ! typical choice of parameter::ten_moment_k0 is
    ! ten_moment_k0 = 1./3. / sqrt(m_s)

    allocate(tau(0:dimX(1),0:dimX(2),0:dimX(3)))
    allocate(change(0:dimX(1),0:dimX(2),0:dimX(3),6))
    allocate(n_inv (-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)))
    allocate(tr_P_m(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)))

    !$acc enter data create(tau, change, n_inv, tr_P_m)

    !$acc data present(tau, change, n_inv, tr_P_m, P_m_interim, fluid_interim, source)

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
                n_inv(x,y,z) = 1._PRC/max(fluid_interim(x,y,z,1), cweno_limiter)

                if (current_cycle == 0) then
                    ! Epsxx - nux * nux * n_inv
                    P_m_interim(x,y,z,1) = fluid_interim(x,y,z,5) - fluid_interim(x,y,z,2)*fluid_interim(x,y,z,2)*n_inv(x,y,z)
                    ! Epsxy - nux * nuy * n_inv
                    P_m_interim(x,y,z,2) = fluid_interim(x,y,z,6) - fluid_interim(x,y,z,2)*fluid_interim(x,y,z,3)*n_inv(x,y,z)
                    ! Epsxz - nux * nuz * n_inv
                    P_m_interim(x,y,z,3) = fluid_interim(x,y,z,7) - fluid_interim(x,y,z,2)*fluid_interim(x,y,z,4)*n_inv(x,y,z)
                    ! Epsyy - nuy * nuy * n_inv
                    P_m_interim(x,y,z,4) = fluid_interim(x,y,z,8) - fluid_interim(x,y,z,3)*fluid_interim(x,y,z,3)*n_inv(x,y,z)
                    ! Epsyz - nuy * nuz * n_inv
                    P_m_interim(x,y,z,5) = fluid_interim(x,y,z,9) - fluid_interim(x,y,z,3)*fluid_interim(x,y,z,4)*n_inv(x,y,z)
                    ! Epszz - nuz * nuz * n_inv
                    P_m_interim(x,y,z,6) = fluid_interim(x,y,z,10) - fluid_interim(x,y,z,4)*fluid_interim(x,y,z,4)*n_inv(x,y,z)
                end if

                if (dimensionality_v == 3) then
                    tr_P_m(x,y,z) = (P_m_interim(x,y,z,1)+P_m_interim(x,y,z,4)+P_m_interim(x,y,z,6))/3._PRC
                else if (dimensionality_v == 2) then
                    tr_P_m(x,y,z) = (P_m_interim(x,y,z,1)+P_m_interim(x,y,z,4))/2._PRC
                else
                    tr_P_m(x,y,z) = P_m_interim(x,y,z,1)
                end if
            end do
        end do
    end do
    !$acc end parallel

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
                tau(x,y,z) = - sqrt(tr_P_m(x,y,z)*n_inv(x,y,z)) / k0 / ncycles
            end do
        end do
    end do
    !$acc end parallel

    ! calculate laplacian
    !$acc parallel
#ifdef ACC_2D_OPT
    !$acc loop gang collapse(3)
    do i = 1,6
        do z = 0,dimX(3)
            do y = 0,dimX(2)
                !$acc loop vector
                do x = 0,dimX(1)
#else
    !$acc loop gang collapse(2)
    do i = 1,6
        do z = 0,dimX(3)
            !$acc loop vector collapse(2)
            do y = 0,dimX(2)
                do x = 0,dimX(1)
#endif
                    ! diagonal entries
                    if (i == 1 .or. (i == 4 .and. dimensionality_v > 1) .or. (i == 6 .and. dimensionality_v > 2)) then
                        change(x,y,z,i) = (-2._PRC*(P_m_interim(x  ,y,z,i)-tr_P_m(x  ,y,z)) + &
                                                   (P_m_interim(x-1,y,z,i)-tr_P_m(x-1,y,z)) + &
                                                   (P_m_interim(x+1,y,z,i)-tr_P_m(x+1,y,z)))/(dx(1)*dx(1))
                    else 
                        change(x,y,z,i) = (-2._PRC*P_m_interim(x  ,y,z,i) + &
                                                   P_m_interim(x-1,y,z,i) + &
                                                   P_m_interim(x+1,y,z,i))/(dx(1)*dx(1))
                    end if
                end do
            end do
        end do
    end do
    !$acc end parallel
    if (dimensionality_x > 1) then
        !$acc parallel
#ifdef ACC_2D_OPT
        !$acc loop gang collapse(3)
        do i = 1,6
            do z = 0,dimX(3)
                do y = 0,dimX(2)
                    !$acc loop vector
                    do x = 0,dimX(1)
#else
        !$acc loop gang collapse(2)
        do i = 1,6
            do z = 0,dimX(3)
                !$acc loop vector collapse(2)
                do y = 0,dimX(2)
                    do x = 0,dimX(1)
#endif
                        ! diagonal entries
                        if (i == 1 .or. (i == 4 .and. dimensionality_v > 1) .or. (i == 6 .and. dimensionality_v > 2)) then
                            change(x,y,z,i) = change(x,y,z,i) + (-2._PRC*(P_m_interim(x  ,y,z,i)-tr_P_m(x  ,y,z)) + &
                                                                         (P_m_interim(x,y-1,z,i)-tr_P_m(x,y-1,z)) + &
                                                                         (P_m_interim(x,y+1,z,i)-tr_P_m(x,y+1,z)))/(dx(2)*dx(2))
                        else 
                            change(x,y,z,i) = change(x,y,z,i) + (-2._PRC*P_m_interim(x  ,y,z,i) + &
                                                                         P_m_interim(x,y-1,z,i) + &
                                                                         P_m_interim(x,y+1,z,i))/(dx(2)*dx(2))
                        end if
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
        do i = 1,6
            do z = 0,dimX(3)
                do y = 0,dimX(2)
                    !$acc loop vector
                    do x = 0,dimX(1)
#else
        !$acc loop gang collapse(2)
        do i = 1,6
            do z = 0,dimX(3)
                !$acc loop vector collapse(2)
                do y = 0,dimX(2)
                    do x = 0,dimX(1)
#endif
                        ! diagonal entries
                        if (i == 1 .or. (i == 4 .and. dimensionality_v > 1) .or. (i == 6 .and. dimensionality_v > 2)) then
                            change(x,y,z,i) = change(x,y,z,i) + (-2._PRC*(P_m_interim(x  ,y,z,i)-tr_P_m(x  ,y,z)) + &
                                                                         (P_m_interim(x,y,z-1,i)-tr_P_m(x,y,z-1)) + &
                                                                         (P_m_interim(x,y,z+1,i)-tr_P_m(x,y,z+1)))/(dx(3)*dx(3))
                        else 
                            change(x,y,z,i) = change(x,y,z,i) + (-2._PRC*P_m_interim(x  ,y,z,i) + &
                                                                         P_m_interim(x,y,z-1,i) + &
                                                                         P_m_interim(x,y,z+1,i))/(dx(3)*dx(3))
                        end if
                    end do
                end do
            end do
        end do
        !$acc end parallel
    end if
    ! end calculate laplacian

    !$acc parallel
#ifdef ACC_2D_OPT
    !$acc loop gang collapse(3)
    do i = 1,6
        do z = 0,dimX(3)
            do y = 0,dimX(2)
                !$acc loop vector
                do x = 0,dimX(1)
#else
    !$acc loop gang collapse(2)
    do i = 1,6
        do z = 0,dimX(3)
            !$acc loop vector collapse(2)
            do y = 0,dimX(2)
                do x = 0,dimX(1)
#endif
                    change(x,y,z,i) = -tau(x,y,z)*change(x,y,z,i)
                    source(x,y,z,i+4) = source(x,y,z,i+4) + change(x,y,z,i)
                end do
            end do
        end do
    end do
    !$acc end parallel

    if (current_cycle < ncycles-1) then
        !$acc parallel
#ifdef ACC_2D_OPT
        !$acc loop gang collapse(3)
        do i = 1,6
            do z = 0,dimX(3)
                do y = 0,dimX(2)
                    !$acc loop vector
                    do x = 0,dimX(1)
#else
        !$acc loop gang collapse(2)
        do i = 1,6
            do z = 0,dimX(3)
                !$acc loop vector collapse(2)
                do y = 0,dimX(2)
                    do x = 0,dimX(1)
#endif
                        P_m_interim(x,y,z,i) = P_m_interim(x,y,z,i) + change(x,y,z,i)*dt
                    end do
                end do
            end do
        end do
        !$acc end parallel
    end if

    !$acc end data

    !$acc exit data delete(tau, change, n_inv, tr_P_m)

    deallocate(n_inv)
    deallocate(tr_P_m)
    deallocate(tau)
    deallocate(change)

end subroutine ten_moment_gradient_p_heat_flux_closure_cycle


subroutine ten_moment_gradient_t_heat_flux_closure_cycle(fluid_interim, source, T_m_interim, k0, &
                        current_cycle, ncycles, cweno_limiter, dimX, BD, dx, dimensionality_x, dt)
    use definitions
    implicit none

    integer :: dimX(3), BD(3), dimensionality_x, current_cycle, ncycles
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:10), &
            intent(in) :: fluid_interim
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:6), &
            intent(inout) :: T_m_interim
    real(kind=PRC), dimension(0:dimX(1),0:dimX(2),0:dimX(3),1:10), intent(inout) :: source
    real(kind=PRC) :: k0, cweno_limiter, dx(3), dt

    real(kind=PRC), allocatable, dimension(:,:,:) :: tau, n_inv, tr_T_m
    real(kind=PRC), allocatable, dimension(:,:,:,:) :: change

    integer :: i, x, y, z

    ! gradient T heat flux closure

    ! typical choice of parameter::ten_moment_k0 is
    ! ten_moment_k0 = 1./(2.*sqrt(2./M_PI)) / sqrt(m_s)
    
    allocate(tau(0:dimX(1),0:dimX(2),0:dimX(3)))
    allocate(change(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),6))
    allocate(n_inv (-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)))
    allocate(tr_T_m(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)))

    !$acc enter data create(tau, change, n_inv, tr_T_m)

    !$acc data present(tau, change, n_inv, tr_T_m, T_m_interim, fluid_interim, source)

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
                n_inv(x,y,z) = 1._PRC/max(fluid_interim(x,y,z,1), cweno_limiter)

                if (current_cycle == 0) then
                    ! Epsxx - nux * nux * n_inv
                    T_m_interim(x,y,z,1) = (fluid_interim(x,y,z,5) - fluid_interim(x,y,z,2)*fluid_interim(x,y,z,2)*n_inv(x,y,z))*n_inv(x,y,z)
                    ! Epsxy - nux * nuy * n_inv
                    T_m_interim(x,y,z,2) = (fluid_interim(x,y,z,6) - fluid_interim(x,y,z,2)*fluid_interim(x,y,z,3)*n_inv(x,y,z))*n_inv(x,y,z)
                    ! Epsxz - nux * nuz * n_inv
                    T_m_interim(x,y,z,3) = (fluid_interim(x,y,z,7) - fluid_interim(x,y,z,2)*fluid_interim(x,y,z,4)*n_inv(x,y,z))*n_inv(x,y,z)
                    ! Epsyy - nuy * nuy * n_inv
                    T_m_interim(x,y,z,4) = (fluid_interim(x,y,z,8) - fluid_interim(x,y,z,3)*fluid_interim(x,y,z,3)*n_inv(x,y,z))*n_inv(x,y,z)
                    ! Epsyz - nuy * nuz * n_inv
                    T_m_interim(x,y,z,5) = (fluid_interim(x,y,z,9) - fluid_interim(x,y,z,3)*fluid_interim(x,y,z,4)*n_inv(x,y,z))*n_inv(x,y,z)
                    ! Epszz - nuz * nuz * n_inv
                    T_m_interim(x,y,z,6) = (fluid_interim(x,y,z,10) - fluid_interim(x,y,z,4)*fluid_interim(x,y,z,4)*n_inv(x,y,z))*n_inv(x,y,z)
                end if

                tr_T_m(x,y,z) = (T_m_interim(x,y,z,1)+T_m_interim(x,y,z,4)+T_m_interim(x,y,z,6))/3._PRC
            end do
        end do
    end do
    !$acc end parallel

    ! calculate laplacian
    !$acc parallel
#ifdef ACC_2D_OPT
    !$acc loop gang collapse(3)
    do i = 1,6
        do z = 0,dimX(3)
            do y = 0,dimX(2)
                !$acc loop vector
                do x = 0,dimX(1)
#else
    !$acc loop gang collapse(2)
    do i = 1,6
        do z = 0,dimX(3)
            !$acc loop vector collapse(2)
            do y = 0,dimX(2)
                do x = 0,dimX(1)
#endif
                    change(x,y,z,i) = (-2._PRC*T_m_interim(x  ,y,z,i) + &
                                               T_m_interim(x-1,y,z,i) + &
                                               T_m_interim(x+1,y,z,i))/(dx(1)*dx(1))
                end do
            end do
        end do
    end do
    !$acc end parallel
    if (dimensionality_x > 1) then
        !$acc parallel
#ifdef ACC_2D_OPT
        !$acc loop gang collapse(3)
        do i = 1,6
            do z = 0,dimX(3)
                do y = 0,dimX(2)
                    !$acc loop vector
                    do x = 0,dimX(1)
#else
        !$acc loop gang collapse(2)
        do i = 1,6
            do z = 0,dimX(3)
                !$acc loop vector collapse(2)
                do y = 0,dimX(2)
                    do x = 0,dimX(1)
#endif
                        change(x,y,z,i) = change(x,y,z,i) + (-2._PRC*T_m_interim(x,y  ,z,i) + &
                                                                     T_m_interim(x,y-1,z,i) + &
                                                                     T_m_interim(x,y+1,z,i))/(dx(2)*dx(2))
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
        do i = 1,6
            do z = 0,dimX(3)
                do y = 0,dimX(2)
                    !$acc loop vector
                    do x = 0,dimX(1)
#else
        !$acc loop gang collapse(2)
        do i = 1,6
            do z = 0,dimX(3)
                !$acc loop vector collapse(2)
                do y = 0,dimX(2)
                    do x = 0,dimX(1)
#endif
                        change(x,y,z,i) = change(x,y,z,i) + (-2._PRC*T_m_interim(x,y,z  ,i) + &
                                                                     T_m_interim(x,y,z-1,i) + &
                                                                     T_m_interim(x,y,z+1,i))/(dx(3)*dx(3))
                    end do
                end do
            end do
        end do
        !$acc end parallel
    end if
    ! end calculate laplacian

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
                !tau(x,y,z) = - fluid_interim(x,y,z,1) / k0 / ncycles
                tau(x,y,z) = - sqrt(tr_T_m(x,y,z)*fluid_interim(x,y,z,1)) / k0 / ncycles
            end do
        end do
    end do
    !$acc end parallel

    !$acc parallel
#ifdef ACC_2D_OPT
    !$acc loop gang collapse(3)
    do i = 1,6
        do z = 0,dimX(3)
            do y = 0,dimX(2)
                !$acc loop vector
                do x = 0,dimX(1)
#else
    !$acc loop gang collapse(2)
    do i = 1,6
        do z = 0,dimX(3)
            !$acc loop vector collapse(2)
            do y = 0,dimX(2)
                do x = 0,dimX(1)
#endif
                    change(x,y,z,i) = -tau(x,y,z)*change(x,y,z,i)
                    source(x,y,z,i+4) = source(x,y,z,i+4) + change(x,y,z,i)
                end do
            end do
        end do
    end do
    !$acc end parallel

    if (current_cycle < ncycles-1) then
        !$acc parallel
#ifdef ACC_2D_OPT
        !$acc loop gang collapse(3)
        do i = 1,6
            do z = 0,dimX(3)
                do y = 0,dimX(2)
                    !$acc loop vector
                    do x = 0,dimX(1)
#else
        !$acc loop gang collapse(2)
        do i = 1,6
            do z = 0,dimX(3)
                !$acc loop vector collapse(2)
                do y = 0,dimX(2)
                    do x = 0,dimX(1)
#endif
                        T_m_interim(x,y,z,i) = T_m_interim(x,y,z,i) + change(x,y,z,i)*n_inv(x,y,z)*dt
                    end do
                end do
            end do
        end do
        !$acc end parallel
    end if

    !$acc end data

    !$acc exit data delete(tau, change, n_inv, tr_T_m)

    deallocate(tau)
    deallocate(change)
    deallocate(n_inv)
    deallocate(tr_T_m)

end subroutine ten_moment_gradient_T_heat_flux_closure_cycle


subroutine ten_moment_gradient_sym_heat_flux_closure_cycle(fluid_interim, source, T_m_interim, k0, &
                        current_cycle, ncycles, cweno_limiter, dimX, BD, dx, dt)
    use definitions
    implicit none

    integer :: dimX(3), BD(3), current_cycle, ncycles
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:10), &
            intent(in) :: fluid_interim
    real(kind=PRC), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),1:6), &
            intent(inout) :: T_m_interim
    real(kind=PRC), dimension(0:dimX(1),0:dimX(2),0:dimX(3),1:10), intent(inout) :: source
    real(kind=PRC) :: k0, cweno_limiter, dx(3), dt

    real(kind=PRC), allocatable, dimension(:,:,:) :: tau, n_inv, tr_T_m
    real(kind=PRC), allocatable, dimension(:,:,:,:) :: change, Q_h

    integer :: i, x, y, z

    ! gradient heat flux closure that symmetrizes grad(T)
    ! Ng, Hakim, Wang, Bhattacharjee 2020
    ! https://doi.org/10.1063/5.0012067

    ! typical choice of parameter::ten_moment_k0 is
    ! ten_moment_k0 = 1. / sqrt(4./(9.*M_PI)) / sqrt(m_s)
    
    allocate(tau(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)))
    allocate(change(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),6))
    allocate(Q_h(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),10))
    allocate(n_inv (-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)))
    allocate(tr_T_m(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)))

    !$acc enter data create(tau, change, Q_h, n_inv, tr_T_m)

    !$acc data present(tau, change, Q_h, n_inv, tr_T_m, T_m_interim, fluid_interim, source)

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
                n_inv(x,y,z) = 1._PRC/max(fluid_interim(x,y,z,1), cweno_limiter)

                if (current_cycle == 0) then
                    ! Epsxx - nux * nux * n_inv
                    T_m_interim(x,y,z,1) = (fluid_interim(x,y,z,5) - fluid_interim(x,y,z,2)*fluid_interim(x,y,z,2)*n_inv(x,y,z))*n_inv(x,y,z)
                    ! Epsxy - nux * nuy * n_inv
                    T_m_interim(x,y,z,2) = (fluid_interim(x,y,z,6) - fluid_interim(x,y,z,2)*fluid_interim(x,y,z,3)*n_inv(x,y,z))*n_inv(x,y,z)
                    ! Epsxz - nux * nuz * n_inv
                    T_m_interim(x,y,z,3) = (fluid_interim(x,y,z,7) - fluid_interim(x,y,z,2)*fluid_interim(x,y,z,4)*n_inv(x,y,z))*n_inv(x,y,z)
                    ! Epsyy - nuy * nuy * n_inv
                    T_m_interim(x,y,z,4) = (fluid_interim(x,y,z,8) - fluid_interim(x,y,z,3)*fluid_interim(x,y,z,3)*n_inv(x,y,z))*n_inv(x,y,z)
                    ! Epsyz - nuy * nuz * n_inv
                    T_m_interim(x,y,z,5) = (fluid_interim(x,y,z,9) - fluid_interim(x,y,z,3)*fluid_interim(x,y,z,4)*n_inv(x,y,z))*n_inv(x,y,z)
                    ! Epszz - nuz * nuz * n_inv
                    T_m_interim(x,y,z,6) = (fluid_interim(x,y,z,10) - fluid_interim(x,y,z,4)*fluid_interim(x,y,z,4)*n_inv(x,y,z))*n_inv(x,y,z)
                end if

                tr_T_m(x,y,z) = (T_m_interim(x,y,z,1)+T_m_interim(x,y,z,4)+T_m_interim(x,y,z,6))/3._PRC
                tau(x,y,z) = -2._PRC*sqrt(2._PRC*tr_T_m(x,y,z)*fluid_interim(x,y,z,1)) / k0
            end do
        end do
    end do
    !$acc end parallel


    !$acc parallel
#ifdef ACC_2D_OPT
    !$acc loop gang collapse(2)
    do z = -BD(3)+1,dimX(3)+BD(3)-1
        do y = -BD(2)+1,dimX(2)+BD(2)-1
            !$acc loop vector
            do x = -BD(1)+1,dimX(1)+BD(1)-1
#else
    !$acc loop gang
    do z = -BD(3)+1,dimX(3)+BD(3)-1
        !$acc loop vector collapse(2)
        do y = -BD(2)+1,dimX(2)+BD(2)-1
            do x = -BD(1)+1,dimX(1)+BD(1)-1
#endif
                ! Q_xxx = 3 d_x T_xx
                Q_h(x,y,z,1 ) = tau(x,y,z)*( &
                                3._PRC*(T_m_interim(x+1,y,z,1) - T_m_interim(x-1,y,z,1))/(2._PRC*dx(1)))
                ! Q_xxy = 2 d_x T_xy + d_y Txx
                Q_h(x,y,z,2 ) = tau(x,y,z)*( &
                                2._PRC*(T_m_interim(x+1,y,z,2) - T_m_interim(x-1,y,z,2))/(2._PRC*dx(1)) + &
                                       (T_m_interim(x,y+1,z,1) - T_m_interim(x,y-1,z,1))/(2._PRC*dx(2)))
                ! Q_xxz = 2 d_x T_xz + d_z Txx
                Q_h(x,y,z,3 ) = tau(x,y,z)*( &
                                2._PRC*(T_m_interim(x+1,y,z,3) - T_m_interim(x-1,y,z,3))/(2._PRC*dx(1)) + &
                                       (T_m_interim(x,y,z+1,1) - T_m_interim(x,y,z-1,1))/(2._PRC*dx(3)))
                ! Q_xyy = 2 d_y T_xy + d_x Tyy
                Q_h(x,y,z,4 ) = tau(x,y,z)*( &
                                2._PRC*(T_m_interim(x,y+1,z,2) - T_m_interim(x,y-1,z,2))/(2._PRC*dx(2)) + &
                                       (T_m_interim(x+1,y,z,4) - T_m_interim(x-1,y,z,4))/(2._PRC*dx(1)))
                ! Q_xyz = d_x T_yz + d_y T_xz + d_z T_xy
                Q_h(x,y,z,5 ) = tau(x,y,z)*( &
                                       (T_m_interim(x+1,y,z,5) - T_m_interim(x-1,y,z,5))/(2._PRC*dx(1)) + &
                                       (T_m_interim(x,y+1,z,3) - T_m_interim(x,y-1,z,3))/(2._PRC*dx(2)) + &
                                       (T_m_interim(x,y,z+1,2) - T_m_interim(x,y,z-1,2))/(2._PRC*dx(3)))
                ! Q_xzz = 2 d_z T_xz + d_x Tzz
                Q_h(x,y,z,6 ) = tau(x,y,z)*( &
                                2._PRC*(T_m_interim(x,y,z+1,3) - T_m_interim(x,y,z-1,3))/(2._PRC*dx(3)) + &
                                       (T_m_interim(x+1,y,z,6) - T_m_interim(x-1,y,z,6))/(2._PRC*dx(1)))
                ! Q_yyy = 3 d_y T_yy
                Q_h(x,y,z,7 ) = tau(x,y,z)*( &
                                3._PRC*(T_m_interim(x,y+1,z,4) - T_m_interim(x,y-1,z,4))/(2._PRC*dx(2)))
                ! Q_yyz = 2 d_y T_yz + d_z Tyy
                Q_h(x,y,z,8 ) = tau(x,y,z)*( &
                                2._PRC*(T_m_interim(x,y+1,z,5) - T_m_interim(x,y-1,z,5))/(2._PRC*dx(2)) + &
                                       (T_m_interim(x,y,z+1,4) - T_m_interim(x,y,z-1,4))/(2._PRC*dx(3)))
                ! Q_yzz = 2 d_z T_yz + d_y Tzz
                Q_h(x,y,z,9 ) = tau(x,y,z)*( &
                                2._PRC*(T_m_interim(x,y,z+1,5) - T_m_interim(x,y,z-1,5))/(2._PRC*dx(3)) + &
                                       (T_m_interim(x,y+1,z,6) - T_m_interim(x,y-1,z,6))/(2._PRC*dx(2)))
                ! Q_zzz = 3 d_z T_zz
                Q_h(x,y,z,10) = tau(x,y,z)*( &
                                3._PRC*(T_m_interim(x,y,z+1,6) - T_m_interim(x,y,z-1,6))/(2._PRC*dx(3)))
            end do
        end do
    end do
    !$acc end parallel


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
                ! (divQ)xx = d/dx Qxxx + d/dy Qxxy + d/dz Qxxz
                change(x,y,z,1) = -((Q_h(x+1,y,z,1) - Q_h(x-1,y,z,1))/(2._PRC*dx(1)) + &
                                    (Q_h(x,y+1,z,2) - Q_h(x,y-1,z,2))/(2._PRC*dx(2)) + &
                                    (Q_h(x,y,z+1,3) - Q_h(x,y,z-1,3))/(2._PRC*dx(3))) / ncycles
                ! (divQ)xy = d/dx Qxyx + d/dy Qxyy + d/dz Qxyz
                change(x,y,z,2) = -((Q_h(x+1,y,z,2) - Q_h(x-1,y,z,2))/(2._PRC*dx(1)) + &
                                    (Q_h(x,y+1,z,4) - Q_h(x,y-1,z,4))/(2._PRC*dx(2)) + &
                                    (Q_h(x,y,z+1,5) - Q_h(x,y,z-1,5))/(2._PRC*dx(3))) / ncycles
                ! (divQ)xz = d/dx Qxzx + d/dy Qxzy + d/dz Qxzz
                change(x,y,z,3) = -((Q_h(x+1,y,z,3) - Q_h(x-1,y,z,3))/(2._PRC*dx(1)) + &
                                    (Q_h(x,y+1,z,5) - Q_h(x,y-1,z,5))/(2._PRC*dx(2)) + &
                                    (Q_h(x,y,z+1,6) - Q_h(x,y,z-1,6))/(2._PRC*dx(3))) / ncycles
                ! (divQ)yy = d/dx Qyyx + d/dy Qyyy + d/dz Qyyz
                change(x,y,z,4) = -((Q_h(x+1,y,z,4) - Q_h(x-1,y,z,4))/(2._PRC*dx(1)) + &
                                    (Q_h(x,y+1,z,7) - Q_h(x,y-1,z,7))/(2._PRC*dx(2)) + &
                                    (Q_h(x,y,z+1,8) - Q_h(x,y,z-1,8))/(2._PRC*dx(3))) / ncycles
                ! (divQ)yz = d/dx Qyzx + d/dy Qyzy + d/dz Qyzz
                change(x,y,z,5) = -((Q_h(x+1,y,z,5) - Q_h(x-1,y,z,5))/(2._PRC*dx(1)) + &
                                    (Q_h(x,y+1,z,8) - Q_h(x,y-1,z,8))/(2._PRC*dx(2)) + &
                                    (Q_h(x,y,z+1,9) - Q_h(x,y,z-1,9))/(2._PRC*dx(3))) / ncycles
                ! (divQ)zz = d/dx Qzzx + d/dy Qzzy + d/dz Qzzz
                change(x,y,z,6) = -((Q_h(x+1,y,z,6) - Q_h(x-1,y,z,6))/(2._PRC*dx(1)) + &
                                    (Q_h(x,y+1,z,9) - Q_h(x,y-1,z,9))/(2._PRC*dx(2)) + &
                                    (Q_h(x,y,z+1,10) - Q_h(x,y,z-1,10))/(2._PRC*dx(3))) / ncycles

            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc parallel
#ifdef ACC_2D_OPT
    !$acc loop gang collapse(3)
    do i = 1,6
        do z = 0,dimX(3)
            do y = 0,dimX(2)
                !$acc loop vector
                do x = 0,dimX(1)
#else
    !$acc loop gang collapse(2)
    do i = 1,6
        do z = 0,dimX(3)
            !$acc loop vector collapse(2)
            do y = 0,dimX(2)
                do x = 0,dimX(1)
#endif
                    source(x,y,z,i+4) = source(x,y,z,i+4) + change(x,y,z,i)
                end do
            end do
        end do
    end do
    !$acc end parallel

    if (current_cycle < ncycles-1) then
        !$acc parallel
#ifdef ACC_2D_OPT
        !$acc loop gang collapse(3)
        do i = 1,6
            do z = 0,dimX(3)
                do y = 0,dimX(2)
                    !$acc loop vector
                    do x = 0,dimX(1)
#else
        !$acc loop gang collapse(2)
        do i = 1,6
            do z = 0,dimX(3)
                !$acc loop vector collapse(2)
                do y = 0,dimX(2)
                    do x = 0,dimX(1)
#endif
                        T_m_interim(x,y,z,i) = T_m_interim(x,y,z,i) + change(x,y,z,i)*n_inv(x,y,z)*dt
                    end do
                end do
            end do
        end do
        !$acc end parallel
    end if

    !$acc end data

    !$acc exit data delete(tau, change, Q_h, n_inv, tr_T_m)

    deallocate(tau)
    deallocate(change)
    deallocate(Q_h)
    deallocate(n_inv)
    deallocate(tr_T_m)

end subroutine ten_moment_gradient_sym_heat_flux_closure_cycle

