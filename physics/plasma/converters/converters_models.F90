
! This file is part of the muphyII multiphysics plasma simulation project.
! 
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at https://mozilla.org/MPL/2.0/.


!_______________ model converters _______________

subroutine convert_five_moments_to_ten_moments(ten_moments, five_moments, array_size, dimensionality_v, BD, skip_boundary)
    use definitions
    implicit none

    integer, intent(in) :: array_size(3), BD(3), dimensionality_v
    real(kind=PRC), intent(out), dimension(array_size(1),array_size(2),array_size(3),10) :: ten_moments
    real(kind=PRC), intent(in), dimension(array_size(1),array_size(2),array_size(3),5) :: five_moments
    logical :: skip_boundary

    integer :: i, x, y, z, BDX, BDY, BDZ

    if (skip_boundary) then
        BDX = BD(1)
        BDY = BD(2)
        BDZ = BD(3)
    else
        BDX = 0
        BDY = 0
        BDZ = 0
    end if

    !$acc data present_or_copyout(ten_moments) present_or_copyin(five_moments)

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 1+BDZ,array_size(3)-BDZ
        do y = 1+BDY,array_size(2)-BDY
            !$acc loop vector
            do x = 1+BDX,array_size(1)-BDX
                do i = 1,4
                    ten_moments(x,y,z,i) = five_moments(x,y,z,i)                ! n, u n
                end do
                ten_moments(x,y,z,5) = five_moments(x,y,z,5) / dimensionality_v ! Eps_xx
                ten_moments(x,y,z,6) = 0._PRC                                   ! Eps_xy
                ten_moments(x,y,z,7) = 0._PRC                                   ! Eps_xz
                if (dimensionality_v > 1) then
                    ten_moments(x,y,z,8) = ten_moments(x,y,z,5)                 ! Eps_yy
                else
                    ten_moments(x,y,z,8) = 0._PRC
                end if
                ten_moments(x,y,z,9) = 0._PRC                                   ! Eps_yz
                if (dimensionality_v > 2) then
                    ten_moments(x,y,z,10) = ten_moments(x,y,z,5)                ! Eps_zz
                else
                    ten_moments(x,y,z,10) = 0._PRC
                end if
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc end data

end subroutine convert_five_moments_to_ten_moments


subroutine convert_ten_moments_to_five_moments(five_moments, ten_moments, array_size, BD, skip_boundary)
    use definitions
    implicit none

    integer, intent(in) :: array_size(3), BD(3)
    real(kind=PRC), intent(out), dimension(array_size(1),array_size(2),array_size(3),5) :: five_moments
    real(kind=PRC), intent(in), dimension(array_size(1),array_size(2),array_size(3),10) :: ten_moments
    logical :: skip_boundary

    integer :: i, x, y, z, BDX, BDY, BDZ

    if (skip_boundary) then
        BDX = BD(1)
        BDY = BD(2)
        BDZ = BD(3)
    else
        BDX = 0
        BDY = 0
        BDZ = 0
    end if

    !$acc data present_or_copyout(five_moments) present_or_copyin(ten_moments)

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 1+BDZ,array_size(3)-BDZ
        do y = 1+BDY,array_size(2)-BDY
            !$acc loop vector
            do x = 1+BDX,array_size(1)-BDX
                do i = 1,4
                    five_moments(x,y,z,i) = ten_moments(x,y,z,i) ! n, u n
                end do
                five_moments(x,y,z,5) = ten_moments(x,y,z,5) + & ! eps = tr(Eps)
                                        ten_moments(x,y,z,8) + &
                                        ten_moments(x,y,z,10)
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc end data

end subroutine convert_ten_moments_to_five_moments


subroutine convert_five_moments_to_f(f, five_moments, m, array_size, dimV, dv, vb, dimensionality_v, BD, skip_boundary)
    use definitions
    implicit none

    integer, intent(in) :: array_size(3), dimV(3), dimensionality_v, BD(3)
    real(kind=PRC), intent(in) :: m, dv(3), vb(3)
    real(kind=PRC), intent(out), dimension(array_size(1),array_size(2),array_size(3),dimV(1),dimV(2),dimV(3)) :: f
    real(kind=PRC), intent(in), dimension(array_size(1),array_size(2),array_size(3),5) :: five_moments
    logical :: skip_boundary

    integer :: x, y, z, vx, vy, vz, BDX, BDY, BDZ
    real(kind=PRC) :: v1, v2, v3, ux, uy, uz, T

    if (skip_boundary) then
        BDX = BD(1)
        BDY = BD(2)
        BDZ = BD(3)
    else
        BDX = 0
        BDY = 0
        BDZ = 0
    end if

    !$acc data present(f) present_or_copyin(five_moments)
    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 1+BDZ,array_size(3)-BDZ
        do y = 1+BDY,array_size(2)-BDY
            !$acc loop vector
            do x = 1+BDX,array_size(1)-BDX
                ux = five_moments(x,y,z,2)/five_moments(x,y,z,1)
                uy = five_moments(x,y,z,3)/five_moments(x,y,z,1)
                uz = five_moments(x,y,z,4)/five_moments(x,y,z,1)
                T = m * (five_moments(x,y,z,5)/five_moments(x,y,z,1) - &
                    (ux*ux + uy*uy + uz*uz)) / 3._PRC
                !$acc loop seq
                do vz = 1,dimV(3)
                    v3 = vb(3) + (vz-0.5_PRC)*dv(3) - uz
                    !$acc loop seq
                    do vy = 1,dimV(2)
                        v2 = vb(2) + (vy-0.5_PRC)*dv(2) - uy
                        !$acc loop seq
                        do vx = 1,dimV(1)
                            v1 = vb(1) + (vx-0.5_PRC)*dv(1) - ux

                            f(x,y,z,vx,vy,vz) = sqrt(m/(T*2._PRC*M_PI))**dimensionality_v * &
                                five_moments(x,y,z,1) * exp(-m*(v1*v1 + v2*v2 + v3*v3)/(2._PRC*T))
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel
    !$acc end data

end subroutine convert_five_moments_to_f


subroutine convert_ten_moments_to_f(f, ten_moments, m, array_size, dimV, dv, vb, dimensionality_v, BD, skip_boundary)
    use definitions
    implicit none

    integer, intent(in) :: array_size(3), dimV(3), dimensionality_v, BD(3)
    real(kind=PRC), intent(in) :: m, dv(3), vb(3)
    real(kind=PRC), intent(out), dimension(array_size(1),array_size(2),array_size(3),dimV(1),dimV(2),dimV(3)) :: f
    real(kind=PRC), intent(in), dimension(array_size(1),array_size(2),array_size(3),10) :: ten_moments
    logical :: skip_boundary

    integer :: x, y, z, vx, vy, vz, BDX, BDY, BDZ
    real(kind=PRC) :: v1, v2, v3, factor
    real(kind=PRC), allocatable, dimension(:,:,:,:) :: u
    real(kind=PRC), allocatable, dimension(:,:,:,:) :: T, inv_T
    real(kind=PRC), allocatable, dimension(:,:,:) :: det_T

    ! f_Maxwell(x,v) = n(x) * (m/(2*pi*k_B))^(3/2) / sqrt(det(T(x))) *
    !                  exp( -m/(2*k_B) * transpose(v-u(x)) inv(T) (v-u(x)) )

    ! TODO: limiter for 1/n?

    if (skip_boundary) then
        BDX = BD(1)
        BDY = BD(2)
        BDZ = BD(3)
    else
        BDX = 0
        BDY = 0
        BDZ = 0
    end if

    allocate(u(array_size(1),array_size(2),array_size(3),3))
    allocate(T(array_size(1),array_size(2),array_size(3),6))
    allocate(inv_T(array_size(1),array_size(2),array_size(3),6))
    allocate(det_T(array_size(1),array_size(2),array_size(3)))

    !$acc enter data create(u, T, inv_T, det_T)

    !$acc data present(f, u, T, inv_T, det_T) present_or_copyin(ten_moments)

    factor = sqrt(m/(2._PRC*M_PI))**dimensionality_v

    if (dimensionality_v == 3) then
        !$acc parallel
        !$acc loop gang collapse(2)
        do z = 1+BDZ,array_size(3)-BDZ
            do y = 1+BDY,array_size(2)-BDY
                !$acc loop vector
                do x = 1+BDX,array_size(1)-BDX
                    u(x,y,z,1) = ten_moments(x,y,z,2)/ten_moments(x,y,z,1) ! u_x
                    u(x,y,z,2) = ten_moments(x,y,z,3)/ten_moments(x,y,z,1) ! u_y
                    u(x,y,z,3) = ten_moments(x,y,z,4)/ten_moments(x,y,z,1) ! u_z

                    T(x,y,z,1) = m*(ten_moments(x,y,z,5 )/ten_moments(x,y,z,1) - u(x,y,z,1)*u(x,y,z,1)) ! T_xx = T_11 = T(1)
                    T(x,y,z,2) = m*(ten_moments(x,y,z,6 )/ten_moments(x,y,z,1) - u(x,y,z,1)*u(x,y,z,2)) ! T_xy = T_12 = T(2)
                    T(x,y,z,3) = m*(ten_moments(x,y,z,7 )/ten_moments(x,y,z,1) - u(x,y,z,1)*u(x,y,z,3)) ! T_xz = T_13 = T(3)
                    T(x,y,z,4) = m*(ten_moments(x,y,z,8 )/ten_moments(x,y,z,1) - u(x,y,z,2)*u(x,y,z,2)) ! T_yy = T_22 = T(4)
                    T(x,y,z,5) = m*(ten_moments(x,y,z,9 )/ten_moments(x,y,z,1) - u(x,y,z,2)*u(x,y,z,3)) ! T_yz = T_23 = T(5)
                    T(x,y,z,6) = m*(ten_moments(x,y,z,10)/ten_moments(x,y,z,1) - u(x,y,z,3)*u(x,y,z,3)) ! T_zz = T_33 = T(6)

                    ! calculate inverse
                    ! inv_T_11 = T_33 * T_22 - T_23^2  
                    ! inv_T_12 = T_13 * T_23 - T_33 * T_12  
                    ! inv_T_13 = T_12 * T_23 - T_13 * T_22  
                    ! inv_T_22 = T_33 * T_11 - T_13^2  
                    ! inv_T_23 = T_12 * T_13 - T_11 * T_23
                    ! inv_T_33 = T_11 * T_22 - T_12^2
                    ! det_T = (T_11 * inv_T_11) + (T_12 * inv_T_12) + (T_13 * inv_T_13)
                    ! inv_T = inv_T / det_T

                    inv_T(x,y,z,1) = T(x,y,z,6)*T(x,y,z,4) - T(x,y,z,5)**2
                    inv_T(x,y,z,2) = T(x,y,z,3)*T(x,y,z,5) - T(x,y,z,6)*T(x,y,z,2)
                    inv_T(x,y,z,3) = T(x,y,z,2)*T(x,y,z,5) - T(x,y,z,3)*T(x,y,z,4)
                    inv_T(x,y,z,4) = T(x,y,z,6)*T(x,y,z,1) - T(x,y,z,3)**2
                    inv_T(x,y,z,5) = T(x,y,z,2)*T(x,y,z,3) - T(x,y,z,1)*T(x,y,z,5)
                    inv_T(x,y,z,6) = T(x,y,z,1)*T(x,y,z,4) - T(x,y,z,2)**2

                    det_T(x,y,z) = T(x,y,z,1)*inv_T(x,y,z,1) + T(x,y,z,2)*inv_T(x,y,z,2) + T(x,y,z,3)*inv_T(x,y,z,3)

                    inv_T(x,y,z,1) = inv_T(x,y,z,1)/det_T(x,y,z)
                    inv_T(x,y,z,2) = inv_T(x,y,z,2)/det_T(x,y,z)
                    inv_T(x,y,z,3) = inv_T(x,y,z,3)/det_T(x,y,z)
                    inv_T(x,y,z,4) = inv_T(x,y,z,4)/det_T(x,y,z)
                    inv_T(x,y,z,5) = inv_T(x,y,z,5)/det_T(x,y,z)
                    inv_T(x,y,z,6) = inv_T(x,y,z,6)/det_T(x,y,z)
                enddo
            enddo
        enddo
        !$acc end parallel

        !$acc parallel
        !$acc loop gang collapse(3)
        do z = 1+BDZ,array_size(3)-BDZ
            do y = 1+BDY,array_size(2)-BDY
                do x = 1+BDX,array_size(1)-BDX
                    !$acc loop vector collapse(3)
                    do vz = 1,dimV(3)
                        do vy = 1,dimV(2)
                            do vx = 1,dimV(1)
                                v3 = vb(3) + (vz-0.5_PRC)*dv(3) - u(x,y,z,3)
                                v2 = vb(2) + (vy-0.5_PRC)*dv(2) - u(x,y,z,2)
                                v1 = vb(1) + (vx-0.5_PRC)*dv(1) - u(x,y,z,1)

                                f(x,y,z,vx,vy,vz) = factor * ten_moments(x,y,z,1) / sqrt(det_T(x,y,z)) * &
                                    exp(-.5_PRC*m*( &
                                    v1*(v1*inv_T(x,y,z,1) + v2*inv_T(x,y,z,2) + v3*inv_T(x,y,z,3)) + &
                                    v2*(v1*inv_T(x,y,z,2) + v2*inv_T(x,y,z,4) + v3*inv_T(x,y,z,5)) + &
                                    v3*(v1*inv_T(x,y,z,3) + v2*inv_T(x,y,z,5) + v3*inv_T(x,y,z,6))))
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        !$acc end parallel

    else if (dimensionality_v == 2) then
        !$acc parallel
        !$acc loop gang collapse(2)
        do z = 1+BDZ,array_size(3)-BDZ
            do y = 1+BDY,array_size(2)-BDY
                !$acc loop vector
                do x = 1+BDX,array_size(1)-BDX
                    u(x,y,z,1) = ten_moments(x,y,z,2)/ten_moments(x,y,z,1) ! u_x
                    u(x,y,z,2) = ten_moments(x,y,z,3)/ten_moments(x,y,z,1) ! u_y

                    T(x,y,z,1) = m*(ten_moments(x,y,z,5 )/ten_moments(x,y,z,1) - u(x,y,z,1)*u(x,y,z,1)) ! T_xx = T_11 = T(1)
                    T(x,y,z,2) = m*(ten_moments(x,y,z,6 )/ten_moments(x,y,z,1) - u(x,y,z,1)*u(x,y,z,2)) ! T_xy = T_12 = T(2)
                    T(x,y,z,3) = m*(ten_moments(x,y,z,8 )/ten_moments(x,y,z,1) - u(x,y,z,2)*u(x,y,z,2)) ! T_yy = T_22 = T(3)

                    ! calculate inverse
                    ! det_T = T_11*T_22 - T_21*T_12
                    det_T(x,y,z) = T(x,y,z,1)*T(x,y,z,3) - T(x,y,z,2)**2
                    ! inv_T_11 =  T_22 / det_T
                    inv_T(x,y,z,1) = T(x,y,z,3) / det_T(x,y,z)
                    ! inv_T_12 = -T_12 / det_T
                    inv_T(x,y,z,2) = -T(x,y,z,2) / det_T(x,y,z)
                    ! inv_T_22 =  T_11 / det_T
                    inv_T(x,y,z,3) = T(x,y,z,1) / det_T(x,y,z)
                enddo
            enddo
        enddo
        !$acc end parallel

        !$acc parallel
        !$acc loop gang collapse(3)
        do z = 1+BDZ,array_size(3)-BDZ
            do y = 1+BDY,array_size(2)-BDY
                do x = 1+BDX,array_size(1)-BDX
                    !$acc loop vector collapse(3)
                    do vz = 1,dimV(3)
                        do vy = 1,dimV(2)
                            do vx = 1,dimV(1)
                                v2 = vb(2) + (vy-0.5_PRC)*dv(2) - u(x,y,z,2)
                                v1 = vb(1) + (vx-0.5_PRC)*dv(1) - u(x,y,z,1)

                                f(x,y,z,vx,vy,vz) = factor * ten_moments(x,y,z,1) / sqrt(det_T(x,y,z)) * &
                                    exp(-.5_PRC*m*( &
                                    v1*(v1*inv_T(x,y,z,1) + v2*inv_T(x,y,z,2)) + &
                                    v2*(v1*inv_T(x,y,z,2) + v2*inv_T(x,y,z,3))))
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        !$acc end parallel

    else if (dimensionality_v == 1) then
        !$acc parallel
        !$acc loop gang collapse(2)
        do z = 1+BDZ,array_size(3)-BDZ
            do y = 1+BDY,array_size(2)-BDY
                !$acc loop vector
                do x = 1+BDX,array_size(1)-BDX
                    u(x,y,z,1) = ten_moments(x,y,z,2)/ten_moments(x,y,z,1) ! u_x
                    T(x,y,z,1) = m*(ten_moments(x,y,z,5)/ten_moments(x,y,z,1) - u(x,y,z,1)*u(x,y,z,1))
                enddo
            enddo
        enddo
        !$acc end parallel

        !$acc parallel
        !$acc loop gang collapse(3)
        do z = 1+BDZ,array_size(3)-BDZ
            do y = 1+BDY,array_size(2)-BDY
                do x = 1+BDX,array_size(1)-BDX
                    !$acc loop vector collapse(3)
                    do vz = 1,dimV(3)
                        do vy = 1,dimV(2)
                            do vx = 1,dimV(1)
                                v1 = vb(1) + (vx-0.5_PRC)*dv(1) - u(x,y,z,1)

                                f(x,y,z,vx,vy,vz) = factor * ten_moments(x,y,z,1) / sqrt(T(x,y,z,1)) * &
                                    exp(-.5_PRC*m*v1*v1/T(x,y,z,1))
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        !$acc end parallel
    end if

    !$acc end data

    !$acc exit data delete(u, T, inv_T, det_T)
    deallocate(u)
    deallocate(T)
    deallocate(inv_T)
    deallocate(det_T)

end subroutine convert_ten_moments_to_f


subroutine convert_f_to_five_moments(five_moments, f, array_size, dimV, dv, vb, dimensionality_v, BD, skip_boundary)
    use definitions
    implicit none

    integer, intent(in) :: array_size(3), dimV(3), dimensionality_v, BD(3)
    real(kind=PRC), intent(in) :: dv(3), vb(3)
    real(kind=PRC), intent(out), dimension(array_size(1),array_size(2),array_size(3),5) :: five_moments
    real(kind=PRC), intent(in), dimension(array_size(1),array_size(2),array_size(3),dimV(1),dimV(2),dimV(3)) :: f
    logical :: skip_boundary

    integer :: x, y, z, vx, vy, vz, BDX, BDY, BDZ
    real(kind=PRC) :: vxVal, vyVal, vzVal, v2, d3v

    d3v = dv(1)
    if (dimensionality_v > 1) d3v = d3v*dv(2)
    if (dimensionality_v > 2) d3v = d3v*dv(3)

    if (skip_boundary) then
        BDX = BD(1)
        BDY = BD(2)
        BDZ = BD(3)
    else
        BDX = 0
        BDY = 0
        BDZ = 0
    end if

    !$acc data present(f) present_or_copyout(five_moments)

    !$acc kernels
    five_moments = 0._PRC
    !$acc end kernels

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 1+BDZ,array_size(3)-BDZ
        do y = 1+BDY,array_size(2)-BDY
            !$acc loop vector
            do x = 1+BDX,array_size(1)-BDX
                !$acc loop seq
                do vz = 1,dimV(3)
                    vzVal = (vb(3)+(vz-.5_PRC)*dv(3))
                    !$acc loop seq
                    do vy = 1,dimV(2)
                        vyVal = (vb(2)+(vy-.5_PRC)*dv(2))
                        !$acc loop seq
                        do vx = 1,dimV(1)
                            vxVal = (vb(1)+(vx-.5_PRC)*dv(1))
                            v2 = vxVal*vxVal + vyVal*vyVal + vzVal*vzVal

                            five_moments(x,y,z,1) = five_moments(x,y,z,1) + f(x,y,z,vx,vy,vz)       ! n
                            five_moments(x,y,z,2) = five_moments(x,y,z,2) + vxVal*f(x,y,z,vx,vy,vz) ! ux
                            five_moments(x,y,z,3) = five_moments(x,y,z,3) + vyVal*f(x,y,z,vx,vy,vz) ! uy
                            five_moments(x,y,z,4) = five_moments(x,y,z,4) + vzVal*f(x,y,z,vx,vy,vz) ! uz
                            five_moments(x,y,z,5) = five_moments(x,y,z,5) + v2*f(x,y,z,vx,vy,vz)    ! eps
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc kernels
    five_moments = five_moments*d3v
    !$acc end kernels

    !$acc end data

end subroutine convert_f_to_five_moments


subroutine convert_f_to_ten_moments(ten_moments, f, array_size, dimV, dv, vb, dimensionality_v, BD, skip_boundary)
    use definitions
    implicit none

    integer, intent(in) :: array_size(3), dimV(3), dimensionality_v, BD(3)
    real(kind=PRC), intent(in) :: dv(3), vb(3)
    real(kind=PRC), intent(out), dimension(array_size(1),array_size(2),array_size(3),10) :: ten_moments
    real(kind=PRC), intent(in), dimension(array_size(1),array_size(2),array_size(3),dimV(1),dimV(2),dimV(3)) :: f
    logical :: skip_boundary

    integer :: x, y, z, vx, vy, vz, BDX, BDY, BDZ
    real(kind=PRC) :: vxVal, vyVal, vzVal, d3v

    d3v = dv(1)
    if (dimensionality_v > 1) d3v = d3v*dv(2)
    if (dimensionality_v > 2) d3v = d3v*dv(3)

    if (skip_boundary) then
        BDX = BD(1)
        BDY = BD(2)
        BDZ = BD(3)
    else
        BDX = 0
        BDY = 0
        BDZ = 0
    end if

    !$acc data present(f) present_or_copyout(ten_moments)

    !$acc kernels
    ten_moments = 0._PRC
    !$acc end kernels

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 1+BDZ,array_size(3)-BDZ
        do y = 1+BDY,array_size(2)-BDY
            !$acc loop vector
            do x = 1+BDX,array_size(1)-BDX
                !$acc loop seq
                do vz = 1,dimV(3)
                    vzVal = (vb(3)+(vz-.5_PRC)*dv(3))
                    !$acc loop seq
                    do vy = 1,dimV(2)
                        vyVal = (vb(2)+(vy-.5_PRC)*dv(2))
                        !$acc loop seq
                        do vx = 1,dimV(1)
                            vxVal = (vb(1)+(vx-.5_PRC)*dv(1))

                            ten_moments(x,y,z,1 ) = ten_moments(x,y,z,1 ) + f(x,y,z,vx,vy,vz)             ! n
                            ten_moments(x,y,z,2 ) = ten_moments(x,y,z,2 ) + vxVal*f(x,y,z,vx,vy,vz)       ! ux
                            ten_moments(x,y,z,3 ) = ten_moments(x,y,z,3 ) + vyVal*f(x,y,z,vx,vy,vz)       ! uy
                            ten_moments(x,y,z,4 ) = ten_moments(x,y,z,4 ) + vzVal*f(x,y,z,vx,vy,vz)       ! uz
                            ten_moments(x,y,z,5 ) = ten_moments(x,y,z,5 ) + vxVal*vxVal*f(x,y,z,vx,vy,vz) ! Eps_xx
                            ten_moments(x,y,z,6 ) = ten_moments(x,y,z,6 ) + vxVal*vyVal*f(x,y,z,vx,vy,vz) ! Eps_xy
                            ten_moments(x,y,z,7 ) = ten_moments(x,y,z,7 ) + vxVal*vzVal*f(x,y,z,vx,vy,vz) ! Eps_xz
                            ten_moments(x,y,z,8 ) = ten_moments(x,y,z,8 ) + vyVal*vyVal*f(x,y,z,vx,vy,vz) ! Eps_yy
                            ten_moments(x,y,z,9 ) = ten_moments(x,y,z,9 ) + vyVal*vzVal*f(x,y,z,vx,vy,vz) ! Eps_yz
                            ten_moments(x,y,z,10) = ten_moments(x,y,z,10) + vzVal*vzVal*f(x,y,z,vx,vy,vz) ! Eps_zz
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc kernels
    ten_moments = ten_moments*d3v
    !$acc end kernels

    !$acc end data

end subroutine convert_f_to_ten_moments


subroutine convert_delta_f_fit(f, n, array_size, m, dimX, dimV, BD, dv, vb, dimensionality_v, direction)
    use definitions
    implicit none

    integer, intent(in) :: array_size(3), dimX(3), dimV(3), BD(3), dimensionality_v, direction
    real(kind=PRC), intent(in) :: m, dv(3), vb(3)
    real(kind=PRC), intent(in), dimension(array_size(1),array_size(2),array_size(3)) :: n
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3), &
                                              dimV(1),dimV(2),dimV(3)) :: f

    real(kind=PRC), allocatable, dimension(:,:,:,:,:,:) :: f_bd, f_maxwell_bd
    real(kind=PRC), allocatable, dimension(:,:,:,:) :: moments_bd
    real(kind=PRC), allocatable, dimension(:,:,:) :: n_new
    integer :: i, x, y, z, vx, vy, vz

    ! estimate value of distribution function in 10mom cells at vlasov-10mom coupling boundary

    ! input f has Vlasov distribution function in the non-boundary cells and
    ! 10mom-Maxwellian distribution function obtained from neighbour 10mom model in the
    ! boundary cells

    ! do constant extrapolation of Vlasov f into boundary cells, then remove
    ! the 10mom-Maxwellian part of the extrapolation to obtain a delta f,
    ! then add this delta f to the 10mom-Maxwellian of the fluid neighbour which
    ! is already present in the boundary cells

    allocate(f_bd(array_size(1),array_size(2),array_size(3),dimV(1),dimV(2),dimV(3)))
    allocate(f_maxwell_bd(array_size(1),array_size(2),array_size(3),dimV(1),dimV(2),dimV(3)))
    allocate(moments_bd(array_size(1),array_size(2),array_size(3),10))
    allocate(n_new(array_size(1),array_size(2),array_size(3)))
    !$acc enter data create(f_bd,f_maxwell_bd,moments_bd,n_new)
    
    !$acc data present(f,f_bd,f_maxwell_bd,moments_bd,n)

    select case (direction)
        case (0) ! x lower
            do i=1,BD(1)
                !$acc kernels
                f_bd(i,:,:,:,:,:) = f(0,:,:,:,:,:)
                !$acc end kernels
            end do
        case (1) ! x upper
            do i=1,BD(1)
                !$acc kernels
                f_bd(i,:,:,:,:,:) = f(dimX(1),:,:,:,:,:)
                !$acc end kernels
            end do
        case (2) ! y lower
            do i=1,BD(2)
                !$acc kernels
                f_bd(:,i,:,:,:,:) = f(:,0,:,:,:,:)
                !$acc end kernels
            end do
        case (3) ! y upper
            do i=1,BD(2)
                !$acc kernels
                f_bd(:,i,:,:,:,:) = f(:,dimX(2),:,:,:,:)
                !$acc end kernels
            end do
        case (4) ! z lower
            do i=1,BD(3)
                !$acc kernels
                f_bd(:,:,i,:,:,:) = f(:,:,0,:,:,:)
                !$acc end kernels
            end do
        case (5) ! z upper
            do i=1,BD(3)
                !$acc kernels
                f_bd(:,:,i,:,:,:) = f(:,:,dimX(3),:,:,:)
                !$acc end kernels
            end do
    end select

    ! TODO not necessary for both values, instead copy
    call convert_f_to_ten_moments(moments_bd, f_bd, array_size, dimV, dv, vb, dimensionality_v, BD, .false.)

    !$acc parallel
    !$acc loop gang collapse(3)
    do z = 1,array_size(3)
        do y = 1,array_size(2)
            do x = 1,array_size(1)
                !$acc loop vector collapse(3)
                do vz = 1,dimV(3)
                    do vy = 1,dimV(2)
                        do vx = 1,dimV(1)
                            f_bd(x,y,z,vx,vy,vz) = f_bd(x,y,z,vx,vy,vz)*n(x,y,z)/moments_bd(x,y,z,1)
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$acc end parallel
    
    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 1,array_size(3)
        do y = 1,array_size(2)
            !$acc loop vector
            do x = 1,array_size(1)
                do i=2,10
                    moments_bd(x,y,z,i) = moments_bd(x,y,z,i)*n(x,y,z)/moments_bd(x,y,z,1)
                end do
                moments_bd(x,y,z,1) = n(x,y,z)
            end do
        end do
    end do
    !$acc end parallel

    call convert_ten_moments_to_f(f_maxwell_bd, moments_bd, m, array_size, dimV, dv, vb, dimensionality_v, BD, .false.)

    !$acc kernels
    f_bd = f_bd - f_maxwell_bd
    !$acc end kernels

    select case (direction)
        case (0) ! x lower
            do i=1,BD(1)
                !$acc kernels
                f_bd(i,:,:,:,:,:) = max(f(-i,:,:,:,:,:) + f_bd(i,:,:,:,:,:)*(BD(1)+1._PRC-i)/(BD(1)+1._PRC), 0._PRC)
                !$acc end kernels
            end do
        case (1) ! x upper
            do i=1,BD(1)
                !$acc kernels
                f_bd(i,:,:,:,:,:) = max(f(dimX(1)+i,:,:,:,:,:) + f_bd(i,:,:,:,:,:)*(BD(1)+1._PRC-i)/(BD(1)+1._PRC), 0._PRC)
                !$acc end kernels
            end do
        case (2) ! y lower
            do i=1,BD(2)
                !$acc kernels
                f_bd(:,i,:,:,:,:) = max(f(:,-i,:,:,:,:) + f_bd(:,i,:,:,:,:)*(BD(2)+1._PRC-i)/(BD(2)+1._PRC), 0._PRC)
                !$acc end kernels
            end do
        case (3) ! y upper
            do i=1,BD(2)
                !$acc kernels
                f_bd(:,i,:,:,:,:) = max(f(:,dimX(2)+i,:,:,:,:) + f_bd(:,i,:,:,:,:)*(BD(2)+1._PRC-i)/(BD(2)+1._PRC), 0._PRC)
                !$acc end kernels
            end do
        case (4) ! z lower
            do i=1,BD(3)
                !$acc kernels
                f_bd(:,:,i,:,:,:) = max(f(:,:,-i,:,:,:) + f_bd(:,:,i,:,:,:)*(BD(3)+1._PRC-i)/(BD(3)+1._PRC), 0._PRC)
                !$acc end kernels
            end do
        case (5) ! z upper
            do i=1,BD(3)
                !$acc kernels
                f_bd(:,:,i,:,:,:) = max(f(:,:,dimX(3)+i,:,:,:) + f_bd(:,:,i,:,:,:)*(BD(3)+1._PRC-i)/(BD(3)+1._PRC), 0._PRC)
                !$acc end kernels
            end do
    end select

    ! rescale distribution function to conserve mass
    call convert_f_to_n(n_new, f_bd, array_size, dimV, dv, dimensionality_v, BD, .false.)

    !$acc parallel
    !$acc loop gang collapse(3)
    do z = 1,array_size(3)
        do y = 1,array_size(2)
            do x = 1,array_size(1)
                !$acc loop vector collapse(3)
                do vz = 1,dimV(3)
                    do vy = 1,dimV(2)
                        do vx = 1,dimV(1)
                            f_bd(x,y,z,vx,vy,vz) = f_bd(x,y,z,vx,vy,vz)*n(x,y,z)/n_new(x,y,z)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel

    select case (direction)
        case (0) ! x lower
            do i=1,BD(1)
                !$acc kernels
                f(-i,:,:,:,:,:) = f_bd(i,:,:,:,:,:)
                !$acc end kernels
            end do
        case (1) ! x upper
            do i=1,BD(1)
                !$acc kernels
                f(dimX(1)+i,:,:,:,:,:) = f_bd(i,:,:,:,:,:)
                !$acc end kernels
            end do
        case (2) ! y lower
            do i=1,BD(2)
                !$acc kernels
                f(:,-i,:,:,:,:) = f_bd(:,i,:,:,:,:)
                !$acc end kernels
            end do
        case (3) ! y upper
            do i=1,BD(2)
                !$acc kernels
                f(:,dimX(2)+i,:,:,:,:) = f_bd(:,i,:,:,:,:)
                !$acc end kernels
            end do
        case (4) ! z lower
            do i=1,BD(3)
                !$acc kernels
                f(:,:,-i,:,:,:) = f_bd(:,:,i,:,:,:)
                !$acc end kernels
            end do
        case (5) ! z upper
            do i=1,BD(3)
                !$acc kernels
                f(:,:,dimX(3)+i,:,:,:) = f_bd(:,:,i,:,:,:)
                !$acc end kernels
            end do
    end select

    !$acc end data

    !$acc exit data delete(f_bd,f_maxwell_bd,moments_bd,n_new)
    deallocate(n_new)
    deallocate(moments_bd)
    deallocate(f_maxwell_bd)
    deallocate(f_bd)

end subroutine convert_delta_f_fit


subroutine convert_delta_f_correction(f, ten_moments, m, dimX, dimV, BD, dv, vb, dimensionality_v)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), dimensionality_v
    real(kind=PRC), intent(in) :: m, dv(3), vb(3)
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3), &
                                              dimV(1),dimV(2),dimV(3)) :: f
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),10) :: ten_moments

    real(kind=PRC), allocatable, dimension(:,:,:,:,:,:) :: f_maxwell
    real(kind=PRC), allocatable, dimension(:,:,:,:) :: ten_moments_vlasov
    real(kind=PRC), allocatable, dimension(:,:,:) :: n_new
    integer :: array_size(3), x, y, z, vx, vy, vz

    ! correct the distribution function to match fluid moments by swapping out the
    ! respective Maxwellians (the maximum entropy part) according to some variation of
    ! f_corrected = f - f_Maxwell_Vlasov + f_Maxwell_Fluid
    ! used to conserve mass, momentum, temperature independent of the Vlasov scheme
    ! improves Vlasov-10mom coupling

    array_size = [dimX(1)+1+2*BD(1), dimX(2)+1+2*BD(2), dimX(3)+1+2*BD(3)]

    allocate(f_maxwell(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3)))
    allocate(ten_moments_vlasov(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),10))
    allocate(n_new(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)))
    
    !$acc data present(f,ten_moments) create(f_maxwell,ten_moments_vlasov,n_new)

    call convert_f_to_ten_moments(ten_moments_vlasov, f, array_size, dimV, dv, vb, dimensionality_v, BD, .true.)
    call convert_ten_moments_to_f(f_maxwell, ten_moments_vlasov, m, array_size, dimV, dv, vb, dimensionality_v, BD, .true.)

    !$acc parallel
    !$acc loop gang collapse(3)
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            do x = 0,dimX(1)
                !$acc loop vector collapse(3)
                do vz = 1,dimV(3)
                    do vy = 1,dimV(2)
                        do vx = 1,dimV(1)
                            ! f = f - f_Maxwell_Vlasov
                            ! force identical densities by rescaling
                            f(x,y,z,vx,vy,vz) = (f(x,y,z,vx,vy,vz) - f_maxwell(x,y,z,vx,vy,vz))*ten_moments(x,y,z,1)/ten_moments_vlasov(x,y,z,1)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel

    call convert_ten_moments_to_f(f_maxwell, ten_moments, m, array_size, dimV, dv, vb, dimensionality_v, BD, .true.)

    !$acc parallel
    !$acc loop gang collapse(3)
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            do x = 0,dimX(1)
                !$acc loop vector collapse(3)
                do vz = 1,dimV(3)
                    do vy = 1,dimV(2)
                        do vx = 1,dimV(1)
                            ! f = f + f_Maxwell_Fluid
                            f(x,y,z,vx,vy,vz) = max(f(x,y,z,vx,vy,vz) + f_maxwell(x,y,z,vx,vy,vz), 0._PRC)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel

    ! rescale distribution function to conserve mass
    call convert_f_to_n(n_new, f, array_size, dimV, dv, dimensionality_v, BD, .true.)

    !$acc parallel
    !$acc loop gang collapse(3)
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            do x = 0,dimX(1)
                !$acc loop vector collapse(3)
                do vz = 1,dimV(3)
                    do vy = 1,dimV(2)
                        do vx = 1,dimV(1)
                            f(x,y,z,vx,vy,vz) = f(x,y,z,vx,vy,vz)*ten_moments(x,y,z,1)/n_new(x,y,z)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc end data

    deallocate(n_new)
    deallocate(ten_moments_vlasov)
    deallocate(f_maxwell)

end subroutine convert_delta_f_correction


subroutine convert_delta_f_correction_keep_vlasov_n(f, ten_moments, m, dimX, dimV, BD, dv, vb, dimensionality_v)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), dimensionality_v
    real(kind=PRC), intent(in) :: m, dv(3), vb(3)
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3), &
                                              dimV(1),dimV(2),dimV(3)) :: f
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),10) :: ten_moments

    real(kind=PRC), allocatable, dimension(:,:,:,:,:,:) :: f_maxwell
    real(kind=PRC), allocatable, dimension(:,:,:,:) :: ten_moments_vlasov
    real(kind=PRC), allocatable, dimension(:,:,:) :: n_new
    integer :: array_size(3), x, y, z, vx, vy, vz

    ! This moment fitting version keeps the conserved Vlasov density (instead of the fluid density)
    ! as used in https://arxiv.org/abs/2109.06743. Thus, it cannot be used when j is calculated
    ! from the CWENO fluxes. Default fitting when j is calculated from interpolation and kept primarily
    ! for reproduction of the paper results (no actual advantages over using the fluid density).

    ! correct the distribution function to match fluid moments by swapping out the
    ! respective Maxwellians (the maximum entropy part) according to some variation of
    ! f_corrected = f - f_Maxwell_Vlasov + f_Maxwell_Fluid
    ! used to conserve mass, momentum, temperature independent of the Vlasov scheme
    ! improves Vlasov-10mom coupling

    array_size = [dimX(1)+1+2*BD(1), dimX(2)+1+2*BD(2), dimX(3)+1+2*BD(3)]

    allocate(f_maxwell(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3)))
    allocate(ten_moments_vlasov(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),10))
    allocate(n_new(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)))
    
    !$acc data present(f,ten_moments) create(f_maxwell,ten_moments_vlasov,n_new)

    call convert_f_to_ten_moments(ten_moments_vlasov, f, array_size, dimV, dv, vb, dimensionality_v, BD, .true.)
    call convert_ten_moments_to_f(f_maxwell, ten_moments_vlasov, m, array_size, dimV, dv, vb, dimensionality_v, BD, .true.)

    !$acc parallel
    !$acc loop gang collapse(3)
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            do x = 0,dimX(1)
                !$acc loop vector collapse(3)
                do vz = 1,dimV(3)
                    do vy = 1,dimV(2)
                        do vx = 1,dimV(1)
                            ! f = f - f_Maxwell_Vlasov
                            ! force identical densities by rescaling
                            f(x,y,z,vx,vy,vz) = (f(x,y,z,vx,vy,vz) - f_maxwell(x,y,z,vx,vy,vz))*ten_moments(x,y,z,1)/ten_moments_vlasov(x,y,z,1)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel

    call convert_ten_moments_to_f(f_maxwell, ten_moments, m, array_size, dimV, dv, vb, dimensionality_v, BD, .true.)

    !$acc parallel
    !$acc loop gang collapse(3)
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            do x = 0,dimX(1)
                !$acc loop vector collapse(3)
                do vz = 1,dimV(3)
                    do vy = 1,dimV(2)
                        do vx = 1,dimV(1)
                            ! f = f + f_Maxwell_Fluid
                            f(x,y,z,vx,vy,vz) = max(f(x,y,z,vx,vy,vz) + f_maxwell(x,y,z,vx,vy,vz), 0._PRC)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel

    ! rescale distribution function to conserve mass
    call convert_f_to_n(n_new, f, array_size, dimV, dv, dimensionality_v, BD, .true.)

    !$acc parallel
    !$acc loop gang collapse(3)
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            do x = 0,dimX(1)
                !$acc loop vector collapse(3)
                do vz = 1,dimV(3)
                    do vy = 1,dimV(2)
                        do vx = 1,dimV(1)
                            f(x,y,z,vx,vy,vz) = f(x,y,z,vx,vy,vz)*ten_moments(x,y,z,1)/n_new(x,y,z)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc end data

    deallocate(n_new)
    deallocate(ten_moments_vlasov)
    deallocate(f_maxwell)

end subroutine convert_delta_f_correction_keep_vlasov_n


subroutine convert_five_ten_moments_fit(ten_moments, dimX, BD, direction)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), direction
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),10) :: ten_moments

    integer :: i, x, y, z, fix
    real(kind=PRC) :: eps_extrapol, eps_bd

    ! estimate value of ten moments in boundary cells at five-ten moments
    ! coupling interface

    ! when entering the function, the ten_moments variable is supposed to
    ! have the neighbour's five moment data in the boundary cells
    ! do constant extrapolation of pressure tensor into boundary cells
    ! and rescale it to match the five moment scalar pressure

    
    !$acc data present(ten_moments)

    select case (direction)
        case (0) ! x lower
            !$acc kernels
            do z = 0,dimX(3)
                do y = 0,dimX(2)
                    do i = 1,BD(1)
                        eps_extrapol = ten_moments(0,y,z,5) + ten_moments(0,y,z,8) + ten_moments(0,y,z,10)
                        eps_bd = ten_moments(-i,y,z,5) + ten_moments(-i,y,z,8) + ten_moments(-i,y,z,10)
                        do fix = 5,10 ! field index
                            ten_moments(-i,y,z,fix) = ten_moments(0,y,z,fix)*eps_bd/eps_extrapol
                        end do
                    end do
                end do
            end do
            !$acc end kernels
        case (1) ! x upper
            !$acc kernels
            do z = 0,dimX(3)
                do y = 0,dimX(2)
                    do i = 1,BD(1)
                        eps_extrapol = ten_moments(dimX(1),y,z,5) + ten_moments(dimX(1),y,z,8) + ten_moments(dimX(1),y,z,10)
                        eps_bd = ten_moments(dimX(1)+i,y,z,5) + ten_moments(dimX(1)+i,y,z,8) + ten_moments(dimX(1)+i,y,z,10)
                        do fix = 5,10
                            ten_moments(dimX(1)+i,y,z,fix) = ten_moments(dimX(1),y,z,fix)*eps_bd/eps_extrapol
                        end do
                    end do
                end do
            end do
            !$acc end kernels
        case (2) ! y lower
            !$acc kernels
            do z = 0,dimX(3)
                do i = 1,BD(2)
                    do x = 0,dimX(1)
                        eps_extrapol = ten_moments(x,0,z,5) + ten_moments(x,0,z,8) + ten_moments(x,0,z,10)
                        eps_bd = ten_moments(x,-i,z,5) + ten_moments(x,-i,z,8) + ten_moments(x,-i,z,10)
                        do fix = 5,10
                            ten_moments(x,-i,z,fix) = ten_moments(x,0,z,fix)*eps_bd/eps_extrapol
                        end do
                    end do
                end do
            end do
            !$acc end kernels
        case (3) ! y upper
            !$acc kernels
            do z = 0,dimX(3)
                do i = 1,BD(2)
                    do x = 0,dimX(1)
                        eps_extrapol = ten_moments(x,dimX(2),z,5) + ten_moments(x,dimX(2),z,8) + ten_moments(x,dimX(2),z,10)
                        eps_bd = ten_moments(x,dimX(2)+i,z,5) + ten_moments(x,dimX(2)+i,z,8) + ten_moments(x,dimX(2)+i,z,10)
                        do fix = 5,10
                            ten_moments(x,dimX(2)+i,z,fix) = ten_moments(x,dimX(2),z,fix)*eps_bd/eps_extrapol
                        end do
                    end do
                end do
            end do
            !$acc end kernels
        case (4) ! z lower
            !$acc kernels
            do i = 1,BD(3)
                do y = 0,dimX(2)
                    do x = 0,dimX(1)
                        eps_extrapol = ten_moments(x,y,0,5) + ten_moments(x,y,0,8) + ten_moments(x,y,0,10)
                        eps_bd = ten_moments(x,y,-i,5) + ten_moments(x,y,-i,8) + ten_moments(x,y,-i,10)
                        do fix = 5,10
                            ten_moments(x,y,-i,fix) = ten_moments(x,y,0,fix)*eps_bd/eps_extrapol
                        end do
                    end do
                end do
            end do
            !$acc end kernels
        case (5) ! z upper
            !$acc kernels
            do i = 1,BD(3)
                do y = 0,dimX(2)
                    do x = 0,dimX(1)
                        eps_extrapol = ten_moments(x,y,dimX(3),5) + ten_moments(x,y,dimX(3),8) + ten_moments(x,y,dimX(3),10)
                        eps_bd = ten_moments(x,y,dimX(3)+i,5) + ten_moments(x,y,dimX(3)+i,8) + ten_moments(x,y,dimX(3)+i,10)
                        do fix = 5,10
                            ten_moments(x,y,dimX(3)+i,fix) = ten_moments(x,y,dimX(3),fix)*eps_bd/eps_extrapol
                        end do
                    end do
                end do
            end do
            !$acc end kernels
    end select

    !$acc end data

end subroutine convert_five_ten_moments_fit

