
! This file is part of the muphyII multiphysics plasma simulation project.
! 
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at https://mozilla.org/MPL/2.0/.


! Interpolation using five point Lagrange polynomials

! The order step_x->step_y->step_z is important because boundary cells
! are only filled, if they are needed by the following steps. This is
! also the reason why step_z does not need the dimensionality_x argument.

subroutine step_x_vlasov_lagrange5(f,dimX,dimV,BD,dx,dv,v_b,dt,dimensionality_x)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3),dimV(3),BD(3),dimensionality_x
    real(kind=PRC), intent(in) :: dt,dx(3),v_b(3),dv(3)
    real(kind=PRC), intent(inout) :: f(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))

    integer :: x,y,z,vx,vy,vz,BDY,BDZ
    real(kind=PRC) :: a,dt_dx,E1,E2,E4,E5,W1,W2,W3,W4,W5,f_minus_2,f_minus_1,f_current

    dt_dx = -dt/dx(1)

    BDY = 0
    BDZ = 0
    if (dimensionality_x > 1) BDY = 2
    if (dimensionality_x > 2) BDZ = 2

    !$acc data present(f)
    !$acc parallel
    !$acc loop gang collapse(3)
    do vz = 1,dimV(3)
        do vy = 1,dimV(2)
            do vx = 1,dimV(1)
                !$acc loop vector collapse(2)
                do z = -BDZ,dimX(3)+BDZ
                    do y = -BDY,dimX(2)+BDY
                        a = (v_b(1) + (vx-0.5_PRC)*dv(1)) * dt_dx
                        ! Lagrange polynomial
                        E1 = a+2;
                        E2 = a+1;
                        E4 = a-1;
                        E5 = a-2;
                        W1 = E2*a*E4*E5/24._PRC
                        W2 = -E1*a*E4*E5/6._PRC
                        W3 = E1*E2*E4*E5/4._PRC
                        W4 = -E1*E2*a*E5/6._PRC
                        W5 = E1*E2*a*E4/24._PRC

                        f_minus_1 = f(-2,y,z,vx,vy,vz)
                        f_current = f(-1,y,z,vx,vy,vz)
                        !$acc loop seq
                        do x = 0,dimX(1)
                            f_minus_2 = f_minus_1
                            f_minus_1 = f_current
                            f_current = f(x,y,z,vx,vy,vz)
                            f(x,y,z,vx,vy,vz) = W1*f_minus_2 + &
                                                W2*f_minus_1 + &
                                                W3*f_current + &
                                                W4*f(x+1,y,z,vx,vy,vz) + &
                                                W5*f(x+2,y,z,vx,vy,vz)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel
    !$acc end data

end subroutine step_x_vlasov_lagrange5


subroutine step_y_vlasov_lagrange5(f,dimX,dimV,BD,dx,dv,v_b,dt,dimensionality_x)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3),dimV(3),BD(3),dimensionality_x
    real(kind=PRC), intent(in) :: dt,dx(3),v_b(3),dv(3)
    real(kind=PRC), intent(inout) :: f(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))

    integer :: x,y,z,vx,vy,vz,BDZ
    real(kind=PRC) :: a,dt_dx,E1,E2,E4,E5,W1,W2,W3,W4,W5,f_minus_2,f_minus_1,f_current

    dt_dx = -dt/dx(2)

    BDZ = 0
    if (dimensionality_x > 2) BDZ = 2

    !$acc data present(f)
    !$acc parallel
    !$acc loop gang collapse(3)
    do vz = 1,dimV(3)
        do vy = 1,dimV(2)
            do vx = 1,dimV(1)
                !$acc loop vector collapse(2)
                do z = -BDZ,dimX(3)+BDZ
                    do x = 0,dimX(1)
                        a = (v_b(2) + (vy-0.5_PRC)*dv(2)) * dt_dx
                        ! Lagrange polynomial
                        E1 = a+2;
                        E2 = a+1;
                        E4 = a-1;
                        E5 = a-2;
                        W1 = E2*a*E4*E5/24._PRC
                        W2 = -E1*a*E4*E5/6._PRC
                        W3 = E1*E2*E4*E5/4._PRC
                        W4 = -E1*E2*a*E5/6._PRC
                        W5 = E1*E2*a*E4/24._PRC

                        f_minus_1 = f(x,-2,z,vx,vy,vz)
                        f_current = f(x,-1,z,vx,vy,vz)
                        !$acc loop seq
                        do y = 0,dimX(2)
                            f_minus_2 = f_minus_1
                            f_minus_1 = f_current
                            f_current = f(x,y,z,vx,vy,vz)
                            f(x,y,z,vx,vy,vz) = W1*f_minus_2 + &
                                                W2*f_minus_1 + &
                                                W3*f_current + &
                                                W4*f(x,y+1,z,vx,vy,vz) + &
                                                W5*f(x,y+2,z,vx,vy,vz)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel
    !$acc end data

end subroutine step_y_vlasov_lagrange5


subroutine step_z_vlasov_lagrange5(f,dimX,dimV,BD,dx,dv,v_b,dt)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3),dimV(3),BD(3)
    real(kind=PRC), intent(in) :: dt,dx(3),v_b(3),dv(3)
    real(kind=PRC), intent(inout) :: f(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))

    integer :: x,y,z,vx,vy,vz
    real(kind=PRC) :: a,dt_dx,E1,E2,E4,E5,W1,W2,W3,W4,W5,f_minus_2,f_minus_1,f_current

    dt_dx = -dt/dx(3)

    !$acc data present(f)
    !$acc parallel
    !$acc loop gang collapse(3)
    do vz = 1,dimV(3)
        do vy = 1,dimV(2)
            do vx = 1,dimV(1)
                !$acc loop vector collapse(2)
                do y = 0,dimX(2)
                    do x = 0,dimX(1)
                        a = (v_b(3) + (vz-0.5_PRC)*dv(3)) * dt_dx
                        ! Lagrange polynomial
                        E1 = a+2;
                        E2 = a+1;
                        E4 = a-1;
                        E5 = a-2;
                        W1 = E2*a*E4*E5/24._PRC
                        W2 = -E1*a*E4*E5/6._PRC
                        W3 = E1*E2*E4*E5/4._PRC
                        W4 = -E1*E2*a*E5/6._PRC
                        W5 = E1*E2*a*E4/24._PRC

                        f_minus_1 = f(x,y,-2,vx,vy,vz)
                        f_current = f(x,y,-1,vx,vy,vz)
                        !$acc loop seq
                        do z = 0,dimX(3)
                            f_minus_2 = f_minus_1
                            f_minus_1 = f_current
                            f_current = f(x,y,z,vx,vy,vz)
                            f(x,y,z,vx,vy,vz) = W1*f_minus_2 + &
                                                W2*f_minus_1 + &
                                                W3*f_current + &
                                                W4*f(x,y,z+1,vx,vy,vz) + &
                                                W5*f(x,y,z+2,vx,vy,vz)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel
    !$acc end data

end subroutine step_z_vlasov_lagrange5


subroutine step_vx_vlasov_lagrange5(f,E,B,dimX,dimV,BD,dv,v_b,q,m,dt,dimensionality_x,include_boundary_cells)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), dimensionality_x
    logical, intent(in) :: include_boundary_cells
    real(kind=PRC), intent(in) :: dt, dv(3), v_b(3), q, m
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: E, B
    real(kind=PRC), intent(inout) :: f(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))

    integer :: x,y,z,vx,vy,vz,BDX,BDY,BDZ
    real(kind=PRC) :: vyVal,vzVal,a,qmd,E1,E2,E4,E5,W1,W2,W3,W4,W5,f_minus_2,f_minus_1,f_current

    qmd = -q/m*dt/dv(1)

    BDX = 0
    BDY = 0
    BDZ = 0
    if (include_boundary_cells) then
        BDX = BD(1)
        if (dimensionality_x > 1) BDY = BD(2)
        if (dimensionality_x > 2) BDZ = BD(3)
    end if

    !$acc data present(f,E,B)
    !$acc parallel
    !$acc loop gang collapse(3)
    do vz = 1,dimV(3)
        do vy = 1,dimV(2)
            do z = -BDZ,dimX(3)+BDZ
                !$acc loop vector collapse(2)
                do y = -BDY,dimX(2)+BDY
                    do x = -BDX,dimX(1)+BDX
                        vyVal = v_b(2) + (vy-0.5_PRC)*dv(2)
                        vzVal = v_b(3) + (vz-0.5_PRC)*dv(3)
                        a = qmd * (E(x,y,z,1) + vyVal*B(x,y,z,3) - vzVal*B(x,y,z,2))
                        ! Lagrange polynomial
                        E1 = a+2;
                        E2 = a+1;
                        E4 = a-1;
                        E5 = a-2;
                        W1 = E2*a*E4*E5/24._PRC
                        W2 = -E1*a*E4*E5/6._PRC
                        W3 = E1*E2*E4*E5/4._PRC
                        W4 = -E1*E2*a*E5/6._PRC
                        W5 = E1*E2*a*E4/24._PRC

                        f_minus_1 = f(x,y,z,1,vy,vz)
                        f_current = f(x,y,z,2,vy,vz)
                        !$acc loop seq
                        do vx = 3,dimV(1)-2
                            f_minus_2 = f_minus_1
                            f_minus_1 = f_current
                            f_current = f(x,y,z,vx,vy,vz)
                            f(x,y,z,vx,vy,vz) = W1*f_minus_2 + &
                                                W2*f_minus_1 + &
                                                W3*f_current + &
                                                W4*f(x,y,z,vx+1,vy,vz) + &
                                                W5*f(x,y,z,vx+2,vy,vz)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel
    !$acc end data

end subroutine step_vx_vlasov_lagrange5


subroutine step_vy_vlasov_lagrange5(f,E,B,dimX,dimV,BD,dv,v_b,q,m,dt,dimensionality_x,include_boundary_cells)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), dimensionality_x
    logical, intent(in) :: include_boundary_cells
    real(kind=PRC), intent(in) :: dt, dv(3), v_b(3), q, m
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: E, B
    real(kind=PRC), intent(inout) :: f(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))

    integer :: x,y,z,vx,vy,vz,BDX,BDY,BDZ
    real(kind=PRC) :: vxVal,vzVal,a,qmd,E1,E2,E4,E5,W1,W2,W3,W4,W5,f_minus_2,f_minus_1,f_current

    qmd = -q/m*dt/dv(2)

    BDX = 0
    BDY = 0
    BDZ = 0
    if (include_boundary_cells) then
        BDX = BD(1)
        if (dimensionality_x > 1) BDY = BD(2)
        if (dimensionality_x > 2) BDZ = BD(3)
    end if

    !$acc data present(f,E,B)
    !$acc parallel
    !$acc loop gang collapse(3)
    do vz = 1,dimV(3)
        do vx = 1,dimV(1)
            do z = -BDZ,dimX(3)+BDZ
                !$acc loop vector collapse(2)
                do y = -BDY,dimX(2)+BDY
                    do x = -BDX,dimX(1)+BDX
                        vxVal = v_b(1) + (vx-0.5_PRC)*dv(1)
                        vzVal = v_b(3) + (vz-0.5_PRC)*dv(3)
                        a = qmd * (E(x,y,z,2) + vzVal*B(x,y,z,1) - vxVal*B(x,y,z,3))
                        ! Lagrange polynomial
                        E1 = a+2;
                        E2 = a+1;
                        E4 = a-1;
                        E5 = a-2;
                        W1 = E2*a*E4*E5/24._PRC
                        W2 = -E1*a*E4*E5/6._PRC
                        W3 = E1*E2*E4*E5/4._PRC
                        W4 = -E1*E2*a*E5/6._PRC
                        W5 = E1*E2*a*E4/24._PRC

                        f_minus_1 = f(x,y,z,vx,1,vz)
                        f_current = f(x,y,z,vx,2,vz)
                        !$acc loop seq
                        do vy = 3,dimV(2)-2
                            f_minus_2 = f_minus_1
                            f_minus_1 = f_current
                            f_current = f(x,y,z,vx,vy,vz)
                            f(x,y,z,vx,vy,vz) = W1*f_minus_2 + &
                                                W2*f_minus_1 + &
                                                W3*f_current + &
                                                W4*f(x,y,z,vx,vy+1,vz) + &
                                                W5*f(x,y,z,vx,vy+2,vz)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel
    !$acc end data

end subroutine step_vy_vlasov_lagrange5


subroutine step_vz_vlasov_lagrange5(f,E,B,dimX,dimV,BD,dv,v_b,q,m,dt,dimensionality_x,include_boundary_cells)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), dimensionality_x
    logical, intent(in) :: include_boundary_cells
    real(kind=PRC), intent(in) :: dt, dv(3), v_b(3), q, m
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: E, B
    real(kind=PRC), intent(inout) :: f(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))

    integer :: x,y,z,vx,vy,vz,BDX,BDY,BDZ
    real(kind=PRC) :: vxVal,vyVal,a,qmd,E1,E2,E4,E5,W1,W2,W3,W4,W5,f_minus_2,f_minus_1,f_current

    qmd = -q/m*dt/dv(3)

    BDX = 0
    BDY = 0
    BDZ = 0
    if (include_boundary_cells) then
        BDX = BD(1)
        if (dimensionality_x > 1) BDY = BD(2)
        if (dimensionality_x > 2) BDZ = BD(3)
    end if

    !$acc data present(f,E,B)
    !$acc parallel
    !$acc loop gang collapse(3)
    do vy = 1,dimV(2)
        do vx = 1,dimV(1)
            do z = -BDZ,dimX(3)+BDZ
                !$acc loop vector collapse(2)
                do y = -BDY,dimX(2)+BDY
                    do x = -BDX,dimX(1)+BDX
                        vxVal = v_b(1) + (vx-0.5_PRC)*dv(1)
                        vyVal = v_b(2) + (vy-0.5_PRC)*dv(2)
                        a = qmd * (E(x,y,z,3) + vxVal*B(x,y,z,2) - vyVal*B(x,y,z,1))
                        ! Lagrange polynomial
                        E1 = a+2;
                        E2 = a+1;
                        E4 = a-1;
                        E5 = a-2;
                        W1 = E2*a*E4*E5/24._PRC
                        W2 = -E1*a*E4*E5/6._PRC
                        W3 = E1*E2*E4*E5/4._PRC
                        W4 = -E1*E2*a*E5/6._PRC
                        W5 = E1*E2*a*E4/24._PRC

                        f_minus_1 = f(x,y,z,vx,vy,1)
                        f_current = f(x,y,z,vx,vy,2)
                        !$acc loop seq
                        do vz = 3,dimV(3)-2
                            f_minus_2 = f_minus_1
                            f_minus_1 = f_current
                            f_current = f(x,y,z,vx,vy,vz)
                            f(x,y,z,vx,vy,vz) = W1*f_minus_2 + &
                                                W2*f_minus_1 + &
                                                W3*f_current + &
                                                W4*f(x,y,z,vx,vy,vz+1) + &
                                                W5*f(x,y,z,vx,vy,vz+2)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel
    !$acc end data

end subroutine step_vz_vlasov_lagrange5

