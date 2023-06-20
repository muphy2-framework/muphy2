
! This file is part of the muphyII multiphysics plasma simulation project.
! 
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at https://mozilla.org/MPL/2.0/.


!_______________ converters for physical quantities _______________

subroutine convert_f_to_n(n, f, array_size, dimV, dv, dimensionality_v, BD, skip_boundary)
    use definitions
    implicit none

    integer, intent(in) :: array_size(3), dimV(3), dimensionality_v, BD(3)
    real(kind=PRC), intent(in) :: dv(3)
    real(kind=PRC), intent(out), dimension(array_size(1),array_size(2),array_size(3)) :: n
    real(kind=PRC), intent(in), dimension(array_size(1),array_size(2),array_size(3),dimV(1),dimV(2),dimV(3)) :: f
    logical :: skip_boundary

    integer :: x, y, z, vx, vy, vz, BDX, BDY, BDZ
    real(kind=PRC) :: d3v

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

    !$acc data present(f,n)

    !$acc kernels
    n = 0._PRC
    !$acc end kernels

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 1+BDZ,array_size(3)-BDZ
        do y = 1+BDY,array_size(2)-BDY
            !$acc loop vector
            do x = 1+BDX,array_size(1)-BDX
                !$acc loop seq
                do vz = 1,dimV(3)
                    !$acc loop seq
                    do vy = 1,dimV(2)
                        !$acc loop seq
                        do vx = 1,dimV(1)
                            n(x,y,z) = n(x,y,z) + f(x,y,z,vx,vy,vz)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc kernels
    n = n * d3v
    !$acc end kernels

    !$acc end data

end subroutine convert_f_to_n


subroutine convert_f_to_rho(rho,f_e,f_i,q_e,q_i,dimX,dimV_e,dimV_i,BD,dv_e,dv_i,dimensionality_v)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV_e(3), dimV_i(3), BD(3), dimensionality_v
    real(kind=PRC), intent(in) :: dv_e(3), dv_i(3), q_e, q_i
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),&
                                        -BD(3):dimX(3)+BD(3),dimV_e(1),dimV_e(2),dimV_e(3)) :: f_e
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),&
                                        -BD(3):dimX(3)+BD(3),dimV_i(1),dimV_i(2),dimV_i(3)) :: f_i
    real(kind=PRC), intent(out), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: rho

    integer :: x, y, z, vx, vy, vz
    real(kind=PRC) :: d3v_q_e, d3v_q_i

    d3v_q_e = dv_e(1)*q_e
    d3v_q_i = dv_i(1)*q_i
    if (dimensionality_v > 1) then
        d3v_q_e = d3v_q_e*dv_e(2)
        d3v_q_i = d3v_q_i*dv_i(2)
    end if
    if (dimensionality_v > 2) then
        d3v_q_e = d3v_q_e*dv_e(3)
        d3v_q_i = d3v_q_i*dv_i(3)
    end if

    !$acc data present(f_e,f_i,rho)

    !$acc kernels
    rho = 0.
    !$acc end kernels

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            !$acc loop vector
            do x = 0,dimX(1)
                !$acc loop seq
                do vz = 1,dimV_e(3)
                    !$acc loop seq
                    do vy = 1,dimV_e(2)
                        !$acc loop seq
                        do vx = 1,dimV_e(1)
                            rho(x,y,z) = rho(x,y,z) + f_e(x,y,z,vx,vy,vz)*d3v_q_e
                        enddo
                    enddo
                enddo

                !$acc loop seq
                do vz = 1,dimV_i(3)
                    !$acc loop seq
                    do vy = 1,dimV_i(2)
                        !$acc loop seq
                        do vx = 1,dimV_i(1)
                            rho(x,y,z) = rho(x,y,z) + f_i(x,y,z,vx,vy,vz)*d3v_q_i
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc end data

end subroutine convert_f_to_rho


subroutine convert_f_to_u(u,f,n,dimX,dimV,BD,dv,vb,dimensionality_v)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), dimensionality_v
    real(kind=PRC), intent(in) :: dv(3), vb(3)
    real(kind=PRC), intent(in) :: f(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))
    real(kind=PRC), intent(inout) :: u(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3)
    real(kind=PRC), intent(in) :: n(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3))

    integer :: x, y, z, vx, vy, vz
    real(kind=PRC) :: vxVal, vyVal, vzVal, d3v

    d3v = dv(1)
    if (dimensionality_v > 1) d3v = d3v*dv(2)
    if (dimensionality_v > 2) d3v = d3v*dv(3)

    !$acc data present(f,u)

    !$acc kernels
    u = 0._PRC
    !$acc end kernels

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            !$acc loop vector
            do x = 0,dimX(1)
                !$acc loop seq
                do vz = 1,dimV(3)
                    vzVal = (vb(3)+(vz-0.5_PRC)*dv(3))
                    !$acc loop seq
                    do vy = 1,dimV(2)
                        vyVal = (vb(2)+(vy-0.5_PRC)*dv(2))
                        !$acc loop seq
                        do vx = 1,dimV(1)
                            vxVal = (vb(1)+(vx-0.5_PRC)*dv(1))

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

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = -BD(3),dimX(3)+BD(3)
        do y = -BD(2),dimX(2)+BD(2)
            !$acc loop vector
            do x = -BD(1),dimX(1)+BD(1)
                u(x,y,z,1) = u(x,y,z,1) / n(x,y,z) * d3v
                u(x,y,z,2) = u(x,y,z,2) / n(x,y,z) * d3v
                u(x,y,z,3) = u(x,y,z,3) / n(x,y,z) * d3v
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc end data

end subroutine convert_f_to_u


subroutine convert_f_to_j(j,f_e,f_i,q_e,q_i,dimX,dimV_e,dimV_i,BD,dv_e,dv_i,vb_e,vb_i,dimensionality_v)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV_e(3), dimV_i(3), BD(3), dimensionality_v
    real(kind=PRC), intent(in) :: dv_e(3), dv_i(3), vb_e(3), vb_i(3), q_e, q_i
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),&
                                        -BD(3):dimX(3)+BD(3),dimV_e(1),dimV_e(2),dimV_e(3)) :: f_e
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),&
                                        -BD(3):dimX(3)+BD(3),dimV_i(1),dimV_i(2),dimV_i(3)) :: f_i
    real(kind=PRC), intent(inout) :: j(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3)

    integer :: x, y, z, vx, vy, vz
    real(kind=PRC) :: vxVal, vyVal, vzVal, d3v_q_e, d3v_q_i

    d3v_q_e = dv_e(1)*q_e
    d3v_q_i = dv_i(1)*q_i
    if (dimensionality_v > 1) then
        d3v_q_e = d3v_q_e*dv_e(2)
        d3v_q_i = d3v_q_i*dv_i(2)
    end if
    if (dimensionality_v > 2) then
        d3v_q_e = d3v_q_e*dv_e(3)
        d3v_q_i = d3v_q_i*dv_i(3)
    end if

    !$acc data present(f_e,f_i,j)

    !$acc kernels
    j = 0._PRC
    !$acc end kernels

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            !$acc loop vector
            do x = 0,dimX(1)
                !$acc loop seq
                do vz = 1,dimV_e(3)
                    vzVal = (vb_e(3)+(vz-0.5_PRC)*dv_e(3))
                    !$acc loop seq
                    do vy = 1,dimV_e(2)
                        vyVal = (vb_e(2)+(vy-0.5_PRC)*dv_e(2))
                        !$acc loop seq
                        do vx = 1,dimV_e(1)
                            vxVal = (vb_e(1)+(vx-0.5_PRC)*dv_e(1))

                            j(x,y,z,1) = j(x,y,z,1) + vxVal*d3v_q_e*f_e(x,y,z,vx,vy,vz);
                            j(x,y,z,2) = j(x,y,z,2) + vyVal*d3v_q_e*f_e(x,y,z,vx,vy,vz);
                            j(x,y,z,3) = j(x,y,z,3) + vzVal*d3v_q_e*f_e(x,y,z,vx,vy,vz);
                        enddo
                    enddo
                enddo

                !$acc loop seq
                do vz = 1,dimV_i(3)
                    vzVal = (vb_i(3)+(vz-0.5_PRC)*dv_i(3))
                    !$acc loop seq
                    do vy = 1,dimV_i(2)
                        vyVal = (vb_i(2)+(vy-0.5_PRC)*dv_i(2))
                        !$acc loop seq
                        do vx = 1,dimV_i(1)
                            vxVal = (vb_i(1)+(vx-0.5_PRC)*dv_i(1))

                            j(x,y,z,1) = j(x,y,z,1) + vxVal*d3v_q_i*f_i(x,y,z,vx,vy,vz);
                            j(x,y,z,2) = j(x,y,z,2) + vyVal*d3v_q_i*f_i(x,y,z,vx,vy,vz);
                            j(x,y,z,3) = j(x,y,z,3) + vzVal*d3v_q_i*f_i(x,y,z,vx,vy,vz);
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc end data

end subroutine convert_f_to_j


subroutine convert_f_to_un(un,f,dimX,dimV,BD,dv,vb,dimensionality_v)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), dimensionality_v
    real(kind=PRC), intent(in) :: dv(3), vb(3)
    real(kind=PRC), intent(in) :: f(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))
    real(kind=PRC), intent(inout) :: un(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3)

    integer :: x, y, z, vx, vy, vz
    real(kind=PRC) :: vxVal, vyVal, vzVal, d3v

    d3v = dv(1)
    if (dimensionality_v > 1) d3v = d3v*dv(2)
    if (dimensionality_v > 2) d3v = d3v*dv(3)

    !$acc data present(f,un)

    !$acc kernels
    un = 0._PRC
    !$acc end kernels

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            !$acc loop vector
            do x = 0,dimX(1)
                !$acc loop seq
                do vz = 1,dimV(3)
                    vzVal = (vb(3)+(vz-0.5_PRC)*dv(3))
                    !$acc loop seq
                    do vy = 1,dimV(2)
                        vyVal = (vb(2)+(vy-0.5_PRC)*dv(2))
                        !$acc loop seq
                        do vx = 1,dimV(1)
                            vxVal = (vb(1)+(vx-0.5_PRC)*dv(1))

                            un(x,y,z,1) = un(x,y,z,1) + vxVal * f(x,y,z,vx,vy,vz)
                            un(x,y,z,2) = un(x,y,z,2) + vyVal * f(x,y,z,vx,vy,vz)
                            un(x,y,z,3) = un(x,y,z,3) + vzVal * f(x,y,z,vx,vy,vz)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel
    
    !$acc kernels
    un = un * d3v
    !$acc end kernels

    !$acc end data

end subroutine convert_f_to_un


subroutine convert_f_to_eps(eps,f,dimX,dimV,BD,dv,vb,dimensionality_v)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), dimensionality_v
    real(kind=PRC), intent(in) :: dv(3), vb(3)
    real(kind=PRC), intent(in) :: f(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))
    real(kind=PRC), intent(inout) :: eps(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3))

    integer :: x, y, z, vx, vy, vz
    real(kind=PRC) :: vxValSq, vyValSq, vzValSq, d3v

    ! eps_ij = \int v_i v_j f dv (second moment; energy density w/o normalization)
    ! eps = tr(eps_ij)
    ! p_ij = m \int (v_i - u_i) (v_j - u_j) f dv (pressure tensor)
    ! p = tr(p_ij)/N (scalar pressure; N dimensionality)

    d3v = dv(1)
    if (dimensionality_v > 1) d3v = d3v*dv(2)
    if (dimensionality_v > 2) d3v = d3v*dv(3)

    !$acc data present(f,eps)

    !$acc kernels
    eps = 0._PRC
    !$acc end kernels

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            !$acc loop vector
            do x = 0,dimX(1)
                !$acc loop seq
                do vz = 1,dimV(3)
                    vzValSq = (vb(3)+(vz-0.5_PRC)*dv(3))
                    vzValSq = vzValSq*vzValSq
                    !$acc loop seq
                    do vy = 1,dimV(2)
                        vyValSq = (vb(2)+(vy-0.5_PRC)*dv(2))
                        vyValSq = vyValSq*vyValSq
                        !$acc loop seq
                        do vx = 1,dimV(1)
                            vxValSq = (vb(1)+(vx-0.5_PRC)*dv(1))
                            vxValSq = vxValSq*vxValSq

                            eps(x,y,z) = eps(x,y,z) + (vxValSq + vyValSq + vzValSq)*&
                                        f(x,y,z,vx,vy,vz)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc kernels
    eps = eps * d3v
    !$acc end kernels

    !$acc end data

end subroutine convert_f_to_eps


subroutine convert_f_to_heatflux(Q_h,f,u,m,dimX,dimV,BD,dv,vb,dimensionality_v)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), dimensionality_v
    real(kind=PRC), intent(in) :: m, dv(3), vb(3)
    real(kind=PRC), intent(in) :: f(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))
    real(kind=PRC), intent(in) :: u(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3)
    real(kind=PRC), intent(inout) :: Q_h(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),10)

    integer :: i, x, y, z, vx, vy, vz
    real(kind=PRC) :: vxVal, vyVal, vzVal, d3v

    d3v = dv(1)
    if (dimensionality_v > 1) d3v = d3v*dv(2)
    if (dimensionality_v > 2) d3v = d3v*dv(3)

    !$acc data present(f) present_or_copy(Q_h) present_or_copyin(u)

    !$acc kernels
    Q_h = 0._PRC
    !$acc end kernels

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            !$acc loop vector
            do x = 0,dimX(1)
                !$acc loop seq
                do vz = 1,dimV(3)
                    vzVal = vb(3) + (vz-0.5_PRC)*dv(3)
                    !$acc loop seq
                    do vy = 1,dimV(2)
                        vyVal = vb(2) + (vy-0.5_PRC)*dv(2)
                        !$acc loop seq
                        do vx = 1,dimV(1)
                            vxVal = vb(1) + (vx-0.5_PRC)*dv(1)

                            Q_h(x,y,z,1 ) = Q_h(x,y,z,1 ) + (vxVal-u(x,y,z,1))*(vxVal-u(x,y,z,1))*(vxVal-u(x,y,z,1)) * f(x,y,z,vx,vy,vz) ! xxx
                            Q_h(x,y,z,2 ) = Q_h(x,y,z,2 ) + (vxVal-u(x,y,z,1))*(vxVal-u(x,y,z,1))*(vyVal-u(x,y,z,2)) * f(x,y,z,vx,vy,vz) ! xxy
                            Q_h(x,y,z,3 ) = Q_h(x,y,z,3 ) + (vxVal-u(x,y,z,1))*(vxVal-u(x,y,z,1))*(vzVal-u(x,y,z,3)) * f(x,y,z,vx,vy,vz) ! xxz
                            Q_h(x,y,z,4 ) = Q_h(x,y,z,4 ) + (vxVal-u(x,y,z,1))*(vyVal-u(x,y,z,2))*(vyVal-u(x,y,z,2)) * f(x,y,z,vx,vy,vz) ! xyy
                            Q_h(x,y,z,5 ) = Q_h(x,y,z,5 ) + (vxVal-u(x,y,z,1))*(vyVal-u(x,y,z,2))*(vzVal-u(x,y,z,3)) * f(x,y,z,vx,vy,vz) ! xyz
                            Q_h(x,y,z,6 ) = Q_h(x,y,z,6 ) + (vxVal-u(x,y,z,1))*(vzVal-u(x,y,z,3))*(vzVal-u(x,y,z,3)) * f(x,y,z,vx,vy,vz) ! xzz
                            Q_h(x,y,z,7 ) = Q_h(x,y,z,7 ) + (vyVal-u(x,y,z,2))*(vyVal-u(x,y,z,2))*(vyVal-u(x,y,z,2)) * f(x,y,z,vx,vy,vz) ! yyy
                            Q_h(x,y,z,8 ) = Q_h(x,y,z,8 ) + (vyVal-u(x,y,z,2))*(vyVal-u(x,y,z,2))*(vzVal-u(x,y,z,3)) * f(x,y,z,vx,vy,vz) ! yyz
                            Q_h(x,y,z,9 ) = Q_h(x,y,z,9 ) + (vyVal-u(x,y,z,2))*(vzVal-u(x,y,z,3))*(vzVal-u(x,y,z,3)) * f(x,y,z,vx,vy,vz) ! yzz
                            Q_h(x,y,z,10) = Q_h(x,y,z,10) + (vzVal-u(x,y,z,3))*(vzVal-u(x,y,z,3))*(vzVal-u(x,y,z,3)) * f(x,y,z,vx,vy,vz) ! zzz
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc kernels
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            do x = 0,dimX(1)
                Q_h(x,y,z,1 ) = Q_h(x,y,z,1 )*m*d3v
                Q_h(x,y,z,2 ) = Q_h(x,y,z,2 )*m*d3v
                Q_h(x,y,z,3 ) = Q_h(x,y,z,3 )*m*d3v
                Q_h(x,y,z,4 ) = Q_h(x,y,z,4 )*m*d3v
                Q_h(x,y,z,5 ) = Q_h(x,y,z,5 )*m*d3v
                Q_h(x,y,z,6 ) = Q_h(x,y,z,6 )*m*d3v
                Q_h(x,y,z,7 ) = Q_h(x,y,z,7 )*m*d3v
                Q_h(x,y,z,8 ) = Q_h(x,y,z,8 )*m*d3v
                Q_h(x,y,z,9 ) = Q_h(x,y,z,9 )*m*d3v
                Q_h(x,y,z,10) = Q_h(x,y,z,10)*m*d3v
            enddo
        enddo
    enddo
    !$acc end kernels

    ! constant extrapolation into boundary cells for coupling
    ! with non-Vlasov schemes

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 1,BD(3)
        do y = -BD(2),dimX(2)+BD(2)
            !$acc loop vector collapse(2)
            do x = -BD(1),dimX(1)+BD(1)
                do i = 1,10
                    Q_h(x,y,-z,i) = Q_h(x,y,0,i)
                    Q_h(x,y,dimX(3)+z,i) = Q_h(x,y,dimX(3),i)
                end do
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = -BD(3),dimX(3)+BD(3)
        do y = 1,BD(2)
            !$acc loop vector collapse(2)
            do x = -BD(1),dimX(1)+BD(1)
                do i = 1,10
                    Q_h(x,-y,z,i) = Q_h(x,0,y,i)
                    Q_h(x,dimX(2)+y,z,i) = Q_h(x,dimX(2),z,i)
                end do
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = -BD(3),dimX(3)+BD(3)
        do y = -BD(2),dimX(2)+BD(2)
            !$acc loop vector collapse(2)
            do x = 1,BD(1)
                do i = 1,10
                    Q_h(-x,y,z,i) = Q_h(0,y,z,i)
                    Q_h(dimX(1)+x,y,z,i) = Q_h(dimX(1),y,z,i)
                end do
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc end data

end subroutine convert_f_to_heatflux


subroutine convert_f_to_raw_heatflux(Q_raw,f,dimX,dimV,BD,dv,vb,&
                                     dimensionality_v,exchange_in_this_direction)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), dimensionality_v
    real(kind=PRC), intent(in) :: dv(3), vb(3)
    real(kind=PRC), intent(in) :: f(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))
    real(kind=PRC), intent(inout) :: Q_raw(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),10)
    logical, intent(in) :: exchange_in_this_direction(6)

    integer :: i, x, y, z, vx, vy, vz, BD_dir(6)
    real(kind=PRC) :: vxVal, vyVal, vzVal, d3v

    ! raw 3rd moment

    d3v = dv(1)
    if (dimensionality_v > 1) d3v = d3v*dv(2)
    if (dimensionality_v > 2) d3v = d3v*dv(3)

    ! loop over boundary cells where boundary exchange of
    ! heat flux is not possible (fluid-Vlasov coupling border)
    do i = 1,6
        if (exchange_in_this_direction(i)) then
            BD_dir(i) = 0
        else
            BD_dir(i) = BD((i+1)/2)
        end if
    end do

    !$acc data present(f) present_or_copy(Q_raw)

    !$acc kernels
    Q_raw = 0._PRC
    !$acc end kernels

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = -BD_dir(5),dimX(3)+BD_dir(6)
        do y = -BD_dir(3),dimX(2)+BD_dir(4)
            !$acc loop vector
            do x = -BD_dir(1),dimX(1)+BD_dir(2)
                !$acc loop seq
                do vz = 1,dimV(3)
                    vzVal = vb(3) + (vz-0.5_PRC)*dv(3)
                    !$acc loop seq
                    do vy = 1,dimV(2)
                        vyVal = vb(2) + (vy-0.5_PRC)*dv(2)
                        !$acc loop seq
                        do vx = 1,dimV(1)
                            vxVal = vb(1) + (vx-0.5_PRC)*dv(1)

                            Q_raw(x,y,z,1 ) = Q_raw(x,y,z,1 ) + vxVal*vxVal*vxVal * f(x,y,z,vx,vy,vz) ! xxx
                            Q_raw(x,y,z,2 ) = Q_raw(x,y,z,2 ) + vxVal*vxVal*vyVal * f(x,y,z,vx,vy,vz) ! xxy
                            Q_raw(x,y,z,3 ) = Q_raw(x,y,z,3 ) + vxVal*vxVal*vzVal * f(x,y,z,vx,vy,vz) ! xxz
                            Q_raw(x,y,z,4 ) = Q_raw(x,y,z,4 ) + vxVal*vyVal*vyVal * f(x,y,z,vx,vy,vz) ! xyy
                            Q_raw(x,y,z,5 ) = Q_raw(x,y,z,5 ) + vxVal*vyVal*vzVal * f(x,y,z,vx,vy,vz) ! xyz
                            Q_raw(x,y,z,6 ) = Q_raw(x,y,z,6 ) + vxVal*vzVal*vzVal * f(x,y,z,vx,vy,vz) ! xzz
                            Q_raw(x,y,z,7 ) = Q_raw(x,y,z,7 ) + vyVal*vyVal*vyVal * f(x,y,z,vx,vy,vz) ! yyy
                            Q_raw(x,y,z,8 ) = Q_raw(x,y,z,8 ) + vyVal*vyVal*vzVal * f(x,y,z,vx,vy,vz) ! yyz
                            Q_raw(x,y,z,9 ) = Q_raw(x,y,z,9 ) + vyVal*vzVal*vzVal * f(x,y,z,vx,vy,vz) ! yzz
                            Q_raw(x,y,z,10) = Q_raw(x,y,z,10) + vzVal*vzVal*vzVal * f(x,y,z,vx,vy,vz) ! zzz
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc kernels
    do z = -BD_dir(5),dimX(3)+BD_dir(6)
        do y = -BD_dir(3),dimX(2)+BD_dir(4)
            do x = -BD_dir(1),dimX(1)+BD_dir(2)
                Q_raw(x,y,z,1 ) = Q_raw(x,y,z,1 )*d3v
                Q_raw(x,y,z,2 ) = Q_raw(x,y,z,2 )*d3v
                Q_raw(x,y,z,3 ) = Q_raw(x,y,z,3 )*d3v
                Q_raw(x,y,z,4 ) = Q_raw(x,y,z,4 )*d3v
                Q_raw(x,y,z,5 ) = Q_raw(x,y,z,5 )*d3v
                Q_raw(x,y,z,6 ) = Q_raw(x,y,z,6 )*d3v
                Q_raw(x,y,z,7 ) = Q_raw(x,y,z,7 )*d3v
                Q_raw(x,y,z,8 ) = Q_raw(x,y,z,8 )*d3v
                Q_raw(x,y,z,9 ) = Q_raw(x,y,z,9 )*d3v
                Q_raw(x,y,z,10) = Q_raw(x,y,z,10)*d3v
            enddo
        enddo
    enddo
    !$acc end kernels

    !$acc end data

end subroutine convert_f_to_raw_heatflux


subroutine convert_u_to_j(j,u_e,u_i,n_e,n_i,dimX,BD,q_e,q_i)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in) :: q_e, q_i
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u_e, u_i
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n_e, n_i
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: j

    integer :: x, y, z

    !$acc data present_or_copy(j, u_e, u_i, n_e, n_i)

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            !$acc loop vector
            do x = 0,dimX(1)
                j(x,y,z,1) = q_e*u_e(x,y,z,1)*n_e(x,y,z) + q_i*u_i(x,y,z,1)*n_i(x,y,z)
                j(x,y,z,2) = q_e*u_e(x,y,z,2)*n_e(x,y,z) + q_i*u_i(x,y,z,2)*n_i(x,y,z)
                j(x,y,z,3) = q_e*u_e(x,y,z,3)*n_e(x,y,z) + q_i*u_i(x,y,z,3)*n_i(x,y,z)
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc end data

end subroutine convert_u_to_j


subroutine convert_un_to_j(j,un_e,un_i,dimX,BD,q_e,q_i)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in) :: q_e, q_i
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: un_e, un_i
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: j

    !$acc data present(j, un_e, un_i)
    !$acc kernels
    j = q_e*un_e + q_i*un_i
    !$acc end kernels
    !$acc end data

end subroutine convert_un_to_j


subroutine convert_un_to_j_reflecting_wall(j_face,un_e_face,un_i_face,un_e_center,un_i_center, &
                                              dimX,BD,q_e,q_i,reflecting_wall)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    logical, intent(in) :: reflecting_wall(6)
    real(kind=PRC), intent(in) :: q_e, q_i
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: &
                                  un_e_face, un_i_face, un_e_center, un_i_center
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: j_face

    integer :: x, y, z, i, offset(6)

    ! calculate j from face values in general, but at reflecting walls from cell
    ! center interpolation in order to avoid boundary effects

    do i = 1,6
        if (reflecting_wall(i)) then
            offset(i) = 3
        else
            offset(i) = 0
        end if
    enddo

    !$acc data present(j_face, un_e_face, un_i_face, un_e_center, un_i_center)

    !$acc kernels
    j_face = q_e*un_e_face + q_i*un_i_face
    j_face = q_e*un_e_face + q_i*un_i_face
    j_face = q_e*un_e_face + q_i*un_i_face
    !$acc end kernels

    if (reflecting_wall(1)) then ! x lower
      !$acc parallel
      !$acc loop gang
      do z = 0,dimX(3)
          !$acc loop vector collapse(2)
          do y = 0,dimX(2)
              do x = 0,offset(1)-1
                  j_face(x,y,z,1) = 0.5_PRC*(q_e*(un_e_center(x,y,z,1) + un_e_center(x-1,y,z,1)) + &
                                             q_i*(un_i_center(x,y,z,1) + un_i_center(x-1,y,z,1)))
                  j_face(x,y,z,2) = 0.5_PRC*(q_e*(un_e_center(x,y,z,2) + un_e_center(x,y-1,z,2)) + &
                                             q_i*(un_i_center(x,y,z,2) + un_i_center(x,y-1,z,2)))
                  j_face(x,y,z,3) = 0.5_PRC*(q_e*(un_e_center(x,y,z,3) + un_e_center(x,y,z-1,3)) + &
                                             q_i*(un_i_center(x,y,z,3) + un_i_center(x,y,z-1,3)))
              enddo
          enddo
      enddo
      !$acc end parallel
    end if

    if (reflecting_wall(2)) then ! x upper
      !$acc parallel
      !$acc loop gang
      do z = 0,dimX(3)
          !$acc loop vector collapse(2)
          do y = 0,dimX(2)
              do x = dimX(1)-offset(2)+1,dimX(1)
                  j_face(x,y,z,1) = 0.5_PRC*(q_e*(un_e_center(x,y,z,1) + un_e_center(x-1,y,z,1)) + &
                                             q_i*(un_i_center(x,y,z,1) + un_i_center(x-1,y,z,1)))
                  j_face(x,y,z,2) = 0.5_PRC*(q_e*(un_e_center(x,y,z,2) + un_e_center(x,y-1,z,2)) + &
                                             q_i*(un_i_center(x,y,z,2) + un_i_center(x,y-1,z,2)))
                  j_face(x,y,z,3) = 0.5_PRC*(q_e*(un_e_center(x,y,z,3) + un_e_center(x,y,z-1,3)) + &
                                             q_i*(un_i_center(x,y,z,3) + un_i_center(x,y,z-1,3)))
              enddo
          enddo
      enddo
      !$acc end parallel
    end if

    if (reflecting_wall(3)) then ! y lower
      !$acc parallel
      !$acc loop gang collapse(2)
      do z = 0,dimX(3)
          do y = 0,offset(3)-1
              !$acc loop vector
              do x = 0,dimX(1)
                  j_face(x,y,z,1) = 0.5_PRC*(q_e*(un_e_center(x,y,z,1) + un_e_center(x-1,y,z,1)) + &
                                             q_i*(un_i_center(x,y,z,1) + un_i_center(x-1,y,z,1)))
                  j_face(x,y,z,2) = 0.5_PRC*(q_e*(un_e_center(x,y,z,2) + un_e_center(x,y-1,z,2)) + &
                                             q_i*(un_i_center(x,y,z,2) + un_i_center(x,y-1,z,2)))
                  j_face(x,y,z,3) = 0.5_PRC*(q_e*(un_e_center(x,y,z,3) + un_e_center(x,y,z-1,3)) + &
                                             q_i*(un_i_center(x,y,z,3) + un_i_center(x,y,z-1,3)))
              enddo
          enddo
      enddo
      !$acc end parallel
    end if

    if (reflecting_wall(4)) then ! y upper
      !$acc parallel
      !$acc loop gang collapse(2)
      do z = 0,dimX(3)
          do y = dimX(2)-offset(4)+1,dimX(2)
              !$acc loop vector
              do x = 0,dimX(1)
                  j_face(x,y,z,1) = 0.5_PRC*(q_e*(un_e_center(x,y,z,1) + un_e_center(x-1,y,z,1)) + &
                                             q_i*(un_i_center(x,y,z,1) + un_i_center(x-1,y,z,1)))
                  j_face(x,y,z,2) = 0.5_PRC*(q_e*(un_e_center(x,y,z,2) + un_e_center(x,y-1,z,2)) + &
                                             q_i*(un_i_center(x,y,z,2) + un_i_center(x,y-1,z,2)))
                  j_face(x,y,z,3) = 0.5_PRC*(q_e*(un_e_center(x,y,z,3) + un_e_center(x,y,z-1,3)) + &
                                             q_i*(un_i_center(x,y,z,3) + un_i_center(x,y,z-1,3)))
              enddo
          enddo
      enddo
      !$acc end parallel
    end if

    if (reflecting_wall(5)) then ! z lower
      !$acc parallel
      !$acc loop gang collapse(2)
      do z = 0,offset(5)-1
          do y = 0,dimX(2)
              !$acc loop vector
              do x = 0,dimX(1)
                  j_face(x,y,z,1) = 0.5_PRC*(q_e*(un_e_center(x,y,z,1) + un_e_center(x-1,y,z,1)) + &
                                             q_i*(un_i_center(x,y,z,1) + un_i_center(x-1,y,z,1)))
                  j_face(x,y,z,2) = 0.5_PRC*(q_e*(un_e_center(x,y,z,2) + un_e_center(x,y-1,z,2)) + &
                                             q_i*(un_i_center(x,y,z,2) + un_i_center(x,y-1,z,2)))
                  j_face(x,y,z,3) = 0.5_PRC*(q_e*(un_e_center(x,y,z,3) + un_e_center(x,y,z-1,3)) + &
                                             q_i*(un_i_center(x,y,z,3) + un_i_center(x,y,z-1,3)))
              enddo
          enddo
      enddo
      !$acc end parallel
    end if

    if (reflecting_wall(6)) then ! z upper
      !$acc parallel
      !$acc loop gang collapse(2)
      do z = dimX(3)-offset(6)+1,dimX(3)
          do y = 0,dimX(2)
              !$acc loop vector
              do x = 0,dimX(1)
                  j_face(x,y,z,1) = 0.5_PRC*(q_e*(un_e_center(x,y,z,1) + un_e_center(x-1,y,z,1)) + &
                                             q_i*(un_i_center(x,y,z,1) + un_i_center(x-1,y,z,1)))
                  j_face(x,y,z,2) = 0.5_PRC*(q_e*(un_e_center(x,y,z,2) + un_e_center(x,y-1,z,2)) + &
                                             q_i*(un_i_center(x,y,z,2) + un_i_center(x,y-1,z,2)))
                  j_face(x,y,z,3) = 0.5_PRC*(q_e*(un_e_center(x,y,z,3) + un_e_center(x,y,z-1,3)) + &
                                             q_i*(un_i_center(x,y,z,3) + un_i_center(x,y,z-1,3)))
              enddo
          enddo
      enddo
      !$acc end parallel
    end if

    !$acc end data

end subroutine convert_un_to_j_reflecting_wall


subroutine convert_un_to_u(u,un,n,dimX,BD)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(out), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u
    real(kind=PRC), intent(in),  dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3))   :: n
    real(kind=PRC), intent(in),  dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: un

    integer :: x, y, z

    !TODO: limiter for 1/n?

    !$acc data present(u, un, n)

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            !$acc loop vector
            do x = 0,dimX(1)
                u(x,y,z,1) = un(x,y,z,1)/n(x,y,z)
                u(x,y,z,2) = un(x,y,z,2)/n(x,y,z)
                u(x,y,z,3) = un(x,y,z,3)/n(x,y,z)
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc end data

end subroutine convert_un_to_u


subroutine convert_heatflux_to_divheatflux(divQ_h,Q_h,dx,dimX,BD)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in) :: dx(3)
    real(kind=PRC), intent(inout) :: divQ_h(0:dimX(1),0:dimX(2),0:dimX(3),6)
    real(kind=PRC), intent(in) :: Q_h(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),10)

    integer :: x, y, z

    ! Q_h components: xxx 1, xxy 2, xxz 3, xyy 4, xyz 5, xzz 6, yyy 7, yyz 8, yzz 9, zzz 10

    !$acc data present(divQ_h,Q_h)

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
                ! second order
                ! (divQ)xx = d/dx Qxxx + d/dy Qxxy + d/dz Qxxz
                divQ_h(x,y,z,1) = (Q_h(x+1,y,z,1) - Q_h(x-1,y,z,1))/(2._PRC*dx(1)) + &
                                  (Q_h(x,y+1,z,2) - Q_h(x,y-1,z,2))/(2._PRC*dx(2)) + &
                                  (Q_h(x,y,z+1,3) - Q_h(x,y,z-1,3))/(2._PRC*dx(3))
                ! (divQ)xy = d/dx Qxyx + d/dy Qxyy + d/dz Qxyz
                divQ_h(x,y,z,2) = (Q_h(x+1,y,z,2) - Q_h(x-1,y,z,2))/(2._PRC*dx(1)) + &
                                  (Q_h(x,y+1,z,4) - Q_h(x,y-1,z,4))/(2._PRC*dx(2)) + &
                                  (Q_h(x,y,z+1,5) - Q_h(x,y,z-1,5))/(2._PRC*dx(3))
                ! (divQ)xz = d/dx Qxzx + d/dy Qxzy + d/dz Qxzz
                divQ_h(x,y,z,3) = (Q_h(x+1,y,z,3) - Q_h(x-1,y,z,3))/(2._PRC*dx(1)) + &
                                  (Q_h(x,y+1,z,5) - Q_h(x,y-1,z,5))/(2._PRC*dx(2)) + &
                                  (Q_h(x,y,z+1,6) - Q_h(x,y,z-1,6))/(2._PRC*dx(3))
                ! (divQ)yy = d/dx Qyyx + d/dy Qyyy + d/dz Qyyz
                divQ_h(x,y,z,4) = (Q_h(x+1,y,z,4) - Q_h(x-1,y,z,4))/(2._PRC*dx(1)) + &
                                  (Q_h(x,y+1,z,7) - Q_h(x,y-1,z,7))/(2._PRC*dx(2)) + &
                                  (Q_h(x,y,z+1,8) - Q_h(x,y,z-1,8))/(2._PRC*dx(3))
                ! (divQ)yz = d/dx Qyzx + d/dy Qyzy + d/dz Qyzz
                divQ_h(x,y,z,5) = (Q_h(x+1,y,z,5) - Q_h(x-1,y,z,5))/(2._PRC*dx(1)) + &
                                  (Q_h(x,y+1,z,8) - Q_h(x,y-1,z,8))/(2._PRC*dx(2)) + &
                                  (Q_h(x,y,z+1,9) - Q_h(x,y,z-1,9))/(2._PRC*dx(3))
                ! (divQ)zz = d/dx Qzzx + d/dy Qzzy + d/dz Qzzz
                divQ_h(x,y,z,6) = (Q_h(x+1,y,z,6) - Q_h(x-1,y,z,6))/(2._PRC*dx(1)) + &
                                  (Q_h(x,y+1,z,9) - Q_h(x,y-1,z,9))/(2._PRC*dx(2)) + &
                                  (Q_h(x,y,z+1,10) - Q_h(x,y,z-1,10))/(2._PRC*dx(3))

            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc end data

end subroutine convert_heatflux_to_divheatflux


subroutine convert_mhd_1temperature_to_electron_pressure(P_e,five_moments,m_e,T_e,T_i,dimX,BD)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in) :: m_e, T_e, T_i
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),5) :: five_moments
    real(kind=PRC), intent(out), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: P_e

    integer :: x,y,z

    do z = -BD(3),dimX(3)+BD(3)
        do y = -BD(2),dimX(2)+BD(2)
            do x = -BD(1),dimX(1)+BD(1)
                P_e(x,y,z) = m_e * (five_moments(x,y,z,5)*T_e/(T_e+T_i) - &
                    (five_moments(x,y,z,2)**2 + five_moments(x,y,z,3)**2 + &
                     five_moments(x,y,z,4)**2) / five_moments(x,y,z,1)) / 3._PRC
            enddo
        enddo
    enddo

end subroutine convert_mhd_1temperature_to_electron_pressure

