
! This file is part of the muphyII multiphysics plasma simulation project.
! 
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at https://mozilla.org/MPL/2.0/.


!_______________ setup converters _______________

subroutine convert_fluid_quantities_to_f(f,n,u,T,m,dimX,dimV,BD,dv,vb,dimensionality_v)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), dimensionality_v
    real(kind=PRC), intent(in) :: dv(3), vb(3), m
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3)) :: f
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n, T
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u

    integer :: x,y,z,vx,vy,vz
    real(kind=PRC) :: vxVal, vyVal, vzVal

    !$acc data present_or_copyout(f) copyin(n,u,T)
    !$acc parallel
    !$acc loop gang collapse(3)
    do z = -BD(3),dimX(3)+BD(3)
        do y = -BD(2),dimX(2)+BD(2)
            do x = -BD(1),dimX(1)+BD(1)
                !$acc loop vector collapse(3)
                do vz = 1,dimV(3)
                    do vy = 1,dimV(2)
                        do vx = 1,dimV(1)
                            vyVal = vb(2) + (vy-0.5_PRC)*dv(2)
                            vzVal = vb(3) + (vz-0.5_PRC)*dv(3)
                            vxVal = vb(1) + (vx-0.5_PRC)*dv(1)
                            f(x,y,z,vx,vy,vz) = sqrt(m/(T(x,y,z)*2._PRC*M_PI))**dimensionality_v * n(x,y,z) * &
                                    exp(-m*((vxVal-u(x,y,z,1))**2+(vyVal-u(x,y,z,2))**2 + &
                                    (vzVal-u(x,y,z,3))**2)/(2._PRC*T(x,y,z)))
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel
    !$acc end data

end subroutine convert_fluid_quantities_to_f


subroutine convert_fluid_quantities_to_f_separate_background(f,n,u,T,n_bg,u_bg,T_bg,m,dimX,dimV,BD,dv,vb,dimensionality_v)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), dimensionality_v
    real(kind=PRC), intent(in) :: dv(3), vb(3), m
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3)) :: f
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n, T, n_bg, T_bg
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u, u_bg

    integer :: x,y,z,vx,vy,vz
    real(kind=PRC) :: vxVal, vyVal, vzVal

    !$acc data present_or_copyout(f) copyin(n,u,T,n_bg,u_bg,T_bg)
    !$acc parallel
    !$acc loop gang collapse(3)
    do z = -BD(3),dimX(3)+BD(3)
        do y = -BD(2),dimX(2)+BD(2)
            do x = -BD(1),dimX(1)+BD(1)
                !$acc loop vector collapse(3)
                do vz = 1,dimV(3)
                    do vy = 1,dimV(2)
                        do vx = 1,dimV(1)
                            vyVal = vb(2) + (vy-0.5_PRC)*dv(2)
                            vzVal = vb(3) + (vz-0.5_PRC)*dv(3)
                            vxVal = vb(1) + (vx-0.5_PRC)*dv(1)
                            f(x,y,z,vx,vy,vz) = sqrt(m/(T(x,y,z)*2._PRC*M_PI))**dimensionality_v * n(x,y,z) * &
                                    exp(-m*((vxVal-u(x,y,z,1))**2+(vyVal-u(x,y,z,2))**2 + &
                                    (vzVal-u(x,y,z,3))**2)/(2._PRC*T(x,y,z))) + &
                                    sqrt(m/(T_bg(x,y,z)*2._PRC*M_PI))**dimensionality_v * n_bg(x,y,z) * &
                                    exp(-m*((vxVal-u_bg(x,y,z,1))**2+(vyVal-u_bg(x,y,z,2))**2 + &
                                    (vzVal-u_bg(x,y,z,3))**2)/(2._PRC*T_bg(x,y,z)))
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel
    !$acc end data

end subroutine convert_fluid_quantities_to_f_separate_background


subroutine convert_fluid_quantities_to_mhd(five_moments,&
        n_e,n_i,u_e,u_i,T_e,T_i,m_e,m_i,dimX,BD)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in) :: m_e, m_i
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),5) :: five_moments
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n_e, n_i, T_e, T_i
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u_e, u_i

    five_moments(:,:,:,1) = ( m_e*n_e + m_i*n_i ) / ( m_e + m_i )    ! rho[n_e,n_i]

    five_moments(:,:,:,2) = ( m_e*n_e*u_e(:,:,:,1) + m_i*n_i*u_i(:,:,:,1) ) ! (u_mhd*rho)[u_e,u_i]
    five_moments(:,:,:,3) = ( m_e*n_e*u_e(:,:,:,2) + m_i*n_i*u_i(:,:,:,2) )
    five_moments(:,:,:,4) = ( m_e*n_e*u_e(:,:,:,3) + m_i*n_i*u_i(:,:,:,3) )

    ! @TODO: use rho in place of n_e and n_i?
    five_moments(:,:,:,5) = 3._PRC*(n_e*T_e+n_i*T_i) + ( five_moments(:,:,:,2)**2 + &
            five_moments(:,:,:,3)**2 + five_moments(:,:,:,4)**2 )/five_moments(:,:,:,1) ! Eps_mhd[T_e,T_i,...]
    ! @TODO: Limiter in 1/rho?

end subroutine convert_fluid_quantities_to_mhd


subroutine convert_fluid_quantities_to_five_moments(five_moments,n,u,T,m,dimX,BD,dimensionality_v)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), dimensionality_v
    real(kind=PRC), intent(in) :: m
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),5) :: five_moments
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n, T
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u

    integer :: x,y,z

    !$acc data present_or_copyout(five_moments) copyin(n, u, T)

    !$acc parallel
    !$acc loop gang
    do z = -BD(3),dimX(3)+BD(3)
        !$acc loop vector collapse(2)
        do y = -BD(2),dimX(2)+BD(2)
            do x = -BD(1),dimX(1)+BD(1)
                five_moments(x,y,z,1) = n(x,y,z)
                five_moments(x,y,z,2) = u(x,y,z,1) * n(x,y,z)
                five_moments(x,y,z,3) = u(x,y,z,2) * n(x,y,z)
                five_moments(x,y,z,4) = u(x,y,z,3) * n(x,y,z)
                five_moments(x,y,z,5) = n(x,y,z) * (dimensionality_v*T(x,y,z)/m + &
                    (u(x,y,z,1)*u(x,y,z,1) + u(x,y,z,2)*u(x,y,z,2) + u(x,y,z,3)*u(x,y,z,3)))
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc end data

end subroutine convert_fluid_quantities_to_five_moments


subroutine convert_fluid_quantities_to_ten_moments(ten_moments,n,u,T,m,dimX,BD,dimensionality_v)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), dimensionality_v
    real(kind=PRC), intent(in) :: m
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),10) :: ten_moments
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n, T
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u

    integer :: x,y,z

    !$acc data present_or_copyout(ten_moments) copyin(n, u, T)

    !$acc parallel
    !$acc loop gang
    do z = -BD(3),dimX(3)+BD(3)
        !$acc loop vector collapse(2)
        do y = -BD(2),dimX(2)+BD(2)
            do x = -BD(1),dimX(1)+BD(1)

                ten_moments(x,y,z,1) = n(x,y,z)
                ten_moments(x,y,z,2) = u(x,y,z,1) * n(x,y,z)
                ten_moments(x,y,z,3) = u(x,y,z,2) * n(x,y,z)
                ten_moments(x,y,z,4) = u(x,y,z,3) * n(x,y,z)
                ten_moments(x,y,z,5) = n(x,y,z) * (T(x,y,z)/m + u(x,y,z,1)*u(x,y,z,1))  ! Eps_xx
                ten_moments(x,y,z,6) = 0._PRC                                           ! Eps_xy
                if (dimensionality_v > 1) then
                  ten_moments(x,y,z,8) = n(x,y,z) * (T(x,y,z)/m + u(x,y,z,2)*u(x,y,z,2))  ! Eps_yy
                else
                  ten_moments(x,y,z,8) = 0._PRC
                end if
                ten_moments(x,y,z,7) = 0._PRC                                           ! Eps_xz
                ten_moments(x,y,z,9) = 0._PRC                                           ! Eps_yz
                if (dimensionality_v > 2) then
                  ten_moments(x,y,z,10) = n(x,y,z) * (T(x,y,z)/m + u(x,y,z,3)*u(x,y,z,3)) ! Eps_zz
                else
                  ten_moments(x,y,z,10) = 0._PRC
                end if
 
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc end data

end subroutine convert_fluid_quantities_to_ten_moments


