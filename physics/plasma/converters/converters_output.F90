
! This file is part of the muphyII multiphysics plasma simulation project.
! 
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at https://mozilla.org/MPL/2.0/.


!_______________ output converters _______________

subroutine convert_f_to_output(f,n,u,T,P,m,dimX,dimV,BD,dv,vb,dimensionality_v)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), dimensionality_v
    real(kind=PRC), intent(in) :: dv(3), vb(3), m
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3)) :: f
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n, T
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),6) :: P

    integer :: x,y,z,vx,vy,vz
    real(kind=PRC) :: vxVal, vyVal, vzVal, d3v

    d3v = dv(1)
    if (dimensionality_v > 1) d3v = d3v*dv(2)
    if (dimensionality_v > 2) d3v = d3v*dv(3)

    n = 0._PRC
    u = 0._PRC
    P = 0._PRC

    !$acc data present(f) copy(n,u,P,T)
    !$acc parallel
    !$acc loop gang collapse(2)
    do z = -BD(3),dimX(3)+BD(3)
        do y = -BD(2),dimX(2)+BD(2)
            !$acc loop vector
            do x = -BD(1),dimX(1)+BD(1)
                !$acc loop seq
                do vz = 1,dimV(3)
                    vzVal = vb(3) + (vz-0.5_PRC)*dv(3)
                    !$acc loop seq
                    do vy = 1,dimV(2)
                        vyVal = vb(2) + (vy-0.5_PRC)*dv(2)
                        !$acc loop seq
                        do vx = 1,dimV(1)
                            vxVal = vb(1) + (vx-0.5_PRC)*dv(1)

                            n(x,y,z) = n(x,y,z) + f(x,y,z,vx,vy,vz)
                            u(x,y,z,1) = u(x,y,z,1) + vxVal * f(x,y,z,vx,vy,vz)
                            u(x,y,z,2) = u(x,y,z,2) + vyVal * f(x,y,z,vx,vy,vz)
                            u(x,y,z,3) = u(x,y,z,3) + vzVal * f(x,y,z,vx,vy,vz)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc kernels
    do z = -BD(3),dimX(3)+BD(3)
        do y = -BD(2),dimX(2)+BD(2)
            do x = -BD(1),dimX(1)+BD(1)
                n(x,y,z) = n(x,y,z)*d3v
                u(x,y,z,1) = u(x,y,z,1)*d3v / n(x,y,z)
                u(x,y,z,2) = u(x,y,z,2)*d3v / n(x,y,z)
                u(x,y,z,3) = u(x,y,z,3)*d3v / n(x,y,z)
            enddo
        enddo
    enddo
    !$acc end kernels

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = -BD(3),dimX(3)+BD(3)
        do y = -BD(2),dimX(2)+BD(2)
            !$acc loop vector
            do x = -BD(1),dimX(1)+BD(1)
                !$acc loop seq
                do vz = 1,dimV(3)
                    vzVal = vb(3) + (vz-0.5_PRC)*dv(3)
                    !$acc loop seq
                    do vy = 1,dimV(2)
                        vyVal = vb(2) + (vy-0.5_PRC)*dv(2)
                        !$acc loop seq
                        do vx = 1,dimV(1)
                            vxVal = vb(1) + (vx-0.5_PRC)*dv(1)

                            P(x,y,z,1) = P(x,y,z,1) + (vxVal-u(x,y,z,1))*(vxVal-u(x,y,z,1)) * f(x,y,z,vx,vy,vz) ! xx
                            P(x,y,z,2) = P(x,y,z,2) + (vxVal-u(x,y,z,1))*(vyVal-u(x,y,z,2)) * f(x,y,z,vx,vy,vz) ! xy
                            P(x,y,z,3) = P(x,y,z,3) + (vxVal-u(x,y,z,1))*(vzVal-u(x,y,z,3)) * f(x,y,z,vx,vy,vz) ! xz
                            P(x,y,z,4) = P(x,y,z,4) + (vyVal-u(x,y,z,2))*(vyVal-u(x,y,z,2)) * f(x,y,z,vx,vy,vz) ! yy
                            P(x,y,z,5) = P(x,y,z,5) + (vyVal-u(x,y,z,2))*(vzVal-u(x,y,z,3)) * f(x,y,z,vx,vy,vz) ! yz
                            P(x,y,z,6) = P(x,y,z,6) + (vzVal-u(x,y,z,3))*(vzVal-u(x,y,z,3)) * f(x,y,z,vx,vy,vz) ! zz
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc kernels
    do z = -BD(3),dimX(3)+BD(3)
        do y = -BD(2),dimX(2)+BD(2)
            do x = -BD(1),dimX(1)+BD(1)
                P(x,y,z,1) = P(x,y,z,1)*m*d3v
                P(x,y,z,2) = P(x,y,z,2)*m*d3v
                P(x,y,z,3) = P(x,y,z,3)*m*d3v
                P(x,y,z,4) = P(x,y,z,4)*m*d3v
                P(x,y,z,5) = P(x,y,z,5)*m*d3v
                P(x,y,z,6) = P(x,y,z,6)*m*d3v

                T(x,y,z) = (P(x,y,z,1)+P(x,y,z,4)+P(x,y,z,6)) / (dimensionality_v*n(x,y,z))
            enddo
        enddo
    enddo
    !$acc end kernels
    !$acc end data

end subroutine convert_f_to_output


subroutine convert_five_moments_to_output(five_moments,n,u,T,P,m,dimX,BD,dimensionality_v)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), dimensionality_v
    real(kind=PRC), intent(in) :: m
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),5) :: five_moments
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n, T
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),6) :: P

    integer :: x,y,z

    !$acc data present(five_moments) copyout(n, u, P, T)

    !$acc parallel
    !$acc loop gang
    do z = -BD(3),dimX(3)+BD(3)
        !$acc loop vector collapse(2)
        do y = -BD(2),dimX(2)+BD(2)
            do x = -BD(1),dimX(1)+BD(1)
                n(x,y,z) = five_moments(x,y,z,1)
                u(x,y,z,1) = five_moments(x,y,z,2) / n(x,y,z)
                u(x,y,z,2) = five_moments(x,y,z,3) / n(x,y,z)
                u(x,y,z,3) = five_moments(x,y,z,4) / n(x,y,z)
                T(x,y,z) = m * (five_moments(x,y,z,5)/n(x,y,z) - &
                    (u(x,y,z,1)*u(x,y,z,1) + u(x,y,z,2)*u(x,y,z,2) + u(x,y,z,3)*u(x,y,z,3))) / dimensionality_v
                P(x,y,z,1) = m * (five_moments(x,y,z,5) - & ! P_xx
                    (u(x,y,z,1)*u(x,y,z,1) + u(x,y,z,2)*u(x,y,z,2) + u(x,y,z,3)*u(x,y,z,3))*n(x,y,z)) / dimensionality_v
                P(x,y,z,2) = 0._PRC         ! P_xy
                P(x,y,z,3) = 0._PRC         ! P_xz
                if (dimensionality_v > 1) then
                    P(x,y,z,4) = P(x,y,z,1) ! P_yy
                else
                    P(x,y,z,4) = 0._PRC
                end if
                P(x,y,z,5) = 0._PRC         ! P_yz
                if (dimensionality_v > 2) then
                    P(x,y,z,6) = P(x,y,z,1) ! P_zz
                else
                    P(x,y,z,6) = 0._PRC
                end if
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc end data

end subroutine convert_five_moments_to_output


subroutine convert_ten_moments_to_output(ten_moments,n,u,T,P,m,dimX,BD,dimensionality_v)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), dimensionality_v
    real(kind=PRC), intent(in) :: m
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),10) :: ten_moments
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n, T
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),6) :: P

    integer :: x,y,z

    !$acc data present(ten_moments) copyout(n, u, P, T)

    !$acc parallel
    !$acc loop gang
    do z = -BD(3),dimX(3)+BD(3)
        !$acc loop vector collapse(2)
        do y = -BD(2),dimX(2)+BD(2)
            do x = -BD(1),dimX(1)+BD(1)
                n(x,y,z) = ten_moments(x,y,z,1)
                u(x,y,z,1) = ten_moments(x,y,z,2) / n(x,y,z)
                u(x,y,z,2) = ten_moments(x,y,z,3) / n(x,y,z)
                u(x,y,z,3) = ten_moments(x,y,z,4) / n(x,y,z)
                P(x,y,z,1) = m * (ten_moments(x,y,z,5)  - n(x,y,z) * u(x,y,z,1)*u(x,y,z,1)) ! xx
                P(x,y,z,2) = m * (ten_moments(x,y,z,6)  - n(x,y,z) * u(x,y,z,1)*u(x,y,z,2)) ! xy
                P(x,y,z,3) = m * (ten_moments(x,y,z,7)  - n(x,y,z) * u(x,y,z,1)*u(x,y,z,3)) ! xz
                P(x,y,z,4) = m * (ten_moments(x,y,z,8)  - n(x,y,z) * u(x,y,z,2)*u(x,y,z,2)) ! yy
                P(x,y,z,5) = m * (ten_moments(x,y,z,9)  - n(x,y,z) * u(x,y,z,2)*u(x,y,z,3)) ! yz
                P(x,y,z,6) = m * (ten_moments(x,y,z,10) - n(x,y,z) * u(x,y,z,3)*u(x,y,z,3)) ! zz
                if (dimensionality_v == 1) then
                    T(x,y,z) = P(x,y,z,1) / n(x,y,z)
                else if (dimensionality_v == 2) then
                    T(x,y,z) = (P(x,y,z,1)+P(x,y,z,4))/2._PRC / n(x,y,z)
                else
                    T(x,y,z) = (P(x,y,z,1)+P(x,y,z,4)+P(x,y,z,6))/3._PRC / n(x,y,z)
                end if
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc end data

end subroutine convert_ten_moments_to_output


