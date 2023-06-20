
! This file is part of the muphyII multiphysics plasma simulation project.
! 
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at https://mozilla.org/MPL/2.0/.


subroutine setup_landau_damping(f_e,f_i,E,B,m_e,T0_e,eps0,alpha,dimX,dimV_e,dimV_i,BD,dx,&
                                vb_e,dv_e,dv_i,dimensionality_x,dimensionality_v,xb_loc)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV_e(3), dimV_i(3), BD(3), dimensionality_x, dimensionality_v
    real(kind=PRC), intent(in) :: vb_e(3), dx(3), dv_e(3), dv_i(3), xb_loc(3), m_e, T0_e, eps0, alpha
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),&
                                              dimV_e(1),dimV_e(2),dimV_e(3)) :: f_e
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),&
                                              dimV_i(1),dimV_i(2),dimV_i(3)) :: f_i
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: E, B

    integer :: x, y, z, vx, vy, vz
    real(kind=PRC) :: xVal, yVal, zVal, vxVal, vyVal, vzVal, v2, alpha_y, alpha_z

    alpha_y = 0._PRC
    alpha_z = 0._PRC
    if (dimensionality_x > 1) alpha_y = alpha
    if (dimensionality_x > 2) alpha_z = alpha

    do vz=1,dimV_e(3)
        vzVal = vb_e(3) + (vz-0.5_PRC)*dv_e(3)
        do vy=1,dimV_e(2)
            vyVal = vb_e(2) + (vy-0.5_PRC)*dv_e(2)
            do vx=1,dimV_e(1)
                vxVal = vb_e(1) + (vx-0.5_PRC)*dv_e(1)
                v2 = vxVal**2+vyVal**2+vzVal**2
                do z=-BD(3),dimX(3)+BD(3)
                    zVal = xb_loc(3) + (z+0.5_PRC)*dx(3)
                    do y=-BD(2),dimX(2)+BD(2)
                        yVal = xb_loc(2) + (y+0.5_PRC)*dx(2)
                        do x=-BD(1),dimX(1)+BD(1)
                            xVal = xb_loc(1) + (x+0.5_PRC)*dx(1)

                            f_e(x,y,z,vx,vy,vz) = sqrt(m_e/(M_PI*2._PRC*T0_e))**dimensionality_v*exp(-m_e*v2/(2._PRC*T0_e)) * &
                                  (1._PRC + alpha*cos(0.5_PRC*xVal)+alpha_y*cos(0.5_PRC*yVal)+alpha_z*cos(0.5_PRC*zVal))
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

    if (dimensionality_v == 1) then
        f_i = 1._PRC/(dimV_i(1)*dimV_i(2)*dimV_i(3))/dv_i(1)
    else if (dimensionality_v == 2) then
        f_i = 1._PRC/(dimV_i(1)*dimV_i(2)*dimV_i(3))/(dv_i(1)*dv_i(2))
    else if (dimensionality_v == 3) then
        f_i = 1._PRC/(dimV_i(1)*dimV_i(2)*dimV_i(3))/(dv_i(1)*dv_i(2)*dv_i(3))
    end if

    do z=-BD(3),dimX(3)+BD(3)
        zVal = xb_loc(3) + (z+0.5d0)*dx(3)
        do y=-BD(2),dimX(2)+BD(2)
            yVal = xb_loc(2) + (y+0.5d0)*dx(2)
            do x=-BD(1),dimX(1)+BD(1)
                xVal = xb_loc(1) + (x+0.5d0)*dx(1)

                ! Yee cube
                ! y,z center; x face
                E(x,y,z,1) = -2._PRC*alpha/eps0*sin(0.5*(xVal-0.5_PRC*dx(1)))
                ! x,z center; y face
                if (dimensionality_x > 1) then
                    E(x,y,z,2) = -2._PRC*alpha_y/eps0*sin(0.5*(yVal-0.5_PRC*dx(2)))
                else
                    E(x,y,z,2) = 0._PRC
                end if
                ! x,y center; z face
                if (dimensionality_x > 2) then
                    E(x,y,z,3) = -2._PRC*alpha_z/eps0*sin(0.5*(zVal-0.5_PRC*dx(3)))
                else
                    E(x,y,z,3) = 0._PRC
                end if
                ! y,z face; x center
                B(x,y,z,1) = 0._PRC
                ! x,z face; y center
                B(x,y,z,2) = 0._PRC
                ! x,y face; z center
                B(x,y,z,3) = 0._PRC

            end do
        end do
    end do

end subroutine setup_landau_damping


subroutine setup_two_stream_instability(f_e,f_i,E,B,m_e,T0_e,alpha,dimX,dimV_e,dimV_i,BD,dx,&
                                        vb_e,dv_e,dv_i,dimensionality_v,xb_loc)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV_e(3), dimV_i(3), BD(3), dimensionality_v
    real(kind=PRC), intent(in) :: vb_e(3), dx(3), dv_e(3), dv_i(3), xb_loc(3), m_e, T0_e, alpha
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),&
                                              dimV_e(1),dimV_e(2),dimV_e(3)) :: f_e
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),&
                                              dimV_i(1),dimV_i(2),dimV_i(3)) :: f_i
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: E, B

    integer :: x, y, z, vx, vy, vz
    real(kind=PRC) :: xVal, yVal, zVal, vxVal, vyVal, vzVal, v2

    do vz=1,dimV_e(3)
        vzVal = vb_e(3) + (vz-0.5_PRC)*dv_e(3)
        do vy=1,dimV_e(2)
            vyVal = vb_e(2) + (vy-0.5_PRC)*dv_e(2)
            do vx=1,dimV_e(1)
                vxVal = vb_e(1) + (vx-0.5_PRC)*dv_e(1)
                v2 = vxVal**2+vyVal**2+vzVal**2
                do z=-BD(3),dimX(3)+BD(3)
                    zVal = xb_loc(3) + (z+0.5_PRC)*dx(3)
                    do y=-BD(2),dimX(2)+BD(2)
                        yVal = xb_loc(2) + (y+0.5_PRC)*dx(2)
                        do x=-BD(1),dimX(1)+BD(1)
                            xVal = xb_loc(1) + (x+0.5_PRC)*dx(1)

                            f_e(x,y,z,vx,vy,vz) = sqrt(m_e/(M_PI*2._PRC*T0_e))**dimensionality_v*exp(-m_e*v2/(2._PRC*T0_e)) * &
                                  v2 * (1._PRC + alpha*cos(0.5_PRC*xVal))
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

    if (dimensionality_v == 1) then
        f_i = 1._PRC/(dimV_i(1)*dimV_i(2)*dimV_i(3))/dv_i(1)
    else if (dimensionality_v == 2) then
        f_i = 1._PRC/(dimV_i(1)*dimV_i(2)*dimV_i(3))/(dv_i(1)*dv_i(2))
    else if (dimensionality_v == 3) then
        f_i = 1._PRC/(dimV_i(1)*dimV_i(2)*dimV_i(3))/(dv_i(1)*dv_i(2)*dv_i(3))
    end if

    E = 0._PRC
    B = 0._PRC

end subroutine setup_two_stream_instability


subroutine setup_harris_sheet(n_e,n_i,u_e,u_i,n_e_bg,n_i_bg,u_e_bg,u_i_bg,T0_e,T0_i,E,B,lambda,psi,&
                              n_bg,guide_field,noise_level,drifting_background,sine_perturbation,&
                              dimX,BD,xb,xe,dx,xb_loc)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in) :: xb(3), xe(3), dx(3), xb_loc(3), T0_e, T0_i, lambda, psi, n_bg, guide_field, noise_level
    logical, intent(in) :: drifting_background, sine_perturbation
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n_e, n_i, n_e_bg, n_i_bg
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u_e, u_i, u_e_bg, u_i_bg, E, B

    integer :: x, y, z, seed_size
    integer, allocatable :: my_seed(:)
    real(kind=PRC) :: xVal, yVal, zVal, xsize, ysize, a, rnd

    ! in Birn et al. 2001 (GEM)
    ! https://doi.org/10.1029/1999JA900449

    ! in Wang, Hakim, Bhattacharjee, Germaschewski 2015 (WHBG)
    ! https://arxiv.org/abs/1409.0262
    ! https://doi.org/10.1063/1.4906063

    ! rng different on each process
    call random_seed(size=seed_size)
    allocate(my_seed(seed_size))
    call random_seed(get=my_seed)
    my_seed(1:min(3,seed_size)) = nint(xb_loc(1:min(3,seed_size))/dx(1:min(3,seed_size)))
    call random_seed(put=my_seed)

    ! initialize quarter domain if it starts/ends at zero
    if ((abs(xb(1)) < 1.e-5_PRC .or. abs(xe(1)) < 1.e-5_PRC) .and. &
        (abs(xb(2)) < 1.e-5_PRC .or. abs(xe(2)) < 1.e-5_PRC)) then
        xsize = (xe(1) - xb(1))*2._PRC
        ysize = (xe(2) - xb(2))*2._PRC
    else ! full domain
        xsize = xe(1) - xb(1)
        ysize = xe(2) - xb(2)
    end if

    do z = -BD(3),dimX(3)+BD(3)
        zVal = xb_loc(3) + (z+0.5_PRC)*dx(3)
        do y = -BD(2),dimX(2)+BD(2)
            yVal = xb_loc(2) + (y+0.5_PRC)*dx(2)
            do x = -BD(1),dimX(1)+BD(1)
                xVal = xb_loc(1) + (x+0.5_PRC)*dx(1)

                call random_number(rnd)

                ! cell center
                n_e(x,y,z) = (1._PRC/cosh(yVal/lambda))**2
                n_i(x,y,z) = (1._PRC/cosh(yVal/lambda))**2
                n_e_bg(x,y,z) = n_bg
                n_i_bg(x,y,z) = n_bg

                u_e(x,y,z,1) = 0._PRC
                u_e(x,y,z,2) = 0._PRC
                u_e_bg(x,y,z,1) = 0._PRC
                u_e_bg(x,y,z,2) = 0._PRC

                u_i(x,y,z,1) = 0._PRC
                u_i(x,y,z,2) = 0._PRC
                u_i_bg(x,y,z,1) = 0._PRC
                u_i_bg(x,y,z,2) = 0._PRC

                if (drifting_background) then
                    u_e(x,y,z,3) =  2._PRC * T0_e * (1._PRC/lambda*(1._PRC/cosh(yVal/lambda))**2) / (n_e(x,y,z)+n_e_bg(x,y,z))
                    u_i(x,y,z,3) = -2._PRC * T0_i * (1._PRC/lambda*(1._PRC/cosh(yVal/lambda))**2) / (n_i(x,y,z)+n_i_bg(x,y,z))
                    u_e_bg(x,y,z,3) = u_e(x,y,z,3)
                    u_i_bg(x,y,z,3) = u_i(x,y,z,3)
                else
                    u_e(x,y,z,3) =  2._PRC * T0_e / lambda
                    u_i(x,y,z,3) = -2._PRC * T0_i / lambda
                    u_e_bg(x,y,z,3) = 0._PRC
                    u_i_bg(x,y,z,3) = 0._PRC
                end if

                ! Yee cube
                ! y,z center; x face
                E(x,y,z,1) = 0._PRC
                ! x,z center; y face
                E(x,y,z,2) = 0._PRC
                ! x,y center; z face
                E(x,y,z,3) = 0._PRC
                ! y,z face; x center
                B(x,y,z,1) = tanh((yVal-0.5_PRC*dx(2))/lambda)
                ! x,z face; y center
                B(x,y,z,2) = 0._PRC
                ! x,y face; z center
                B(x,y,z,3) = guide_field

                ! add perturbation
                if (sine_perturbation) then ! GEM-like
                    ! y,z face; x center
                    B(x,y,z,1) = B(x,y,z,1) - psi * M_PI/ysize * &
                        cos(2._PRC*M_PI/xsize*xVal) * sin(M_PI/ysize*(yVal-0.5_PRC*dx(2)))
                    ! x,z face; y center
                    B(x,y,z,2) = B(x,y,z,2) + psi * 2._PRC*M_PI/xsize * &
                        sin(2._PRC*M_PI/xsize*(xVal-0.5_PRC*dx(1))) * cos(M_PI/ysize*yVal)

                    if (drifting_background) then ! neglected in other cases for simplicity
                        u_e(x,y,z,3) = u_e(x,y,z,3) - 2._PRC * T0_e * ((2._PRC*M_PI/xsize)**2 + (M_PI/ysize)**2) * psi * &
                            cos(2._PRC*M_PI*xVal/xsize)*cos(M_PI*yVal/ysize) / (n_e(x,y,z)+n_e_bg(x,y,z))
                        u_i(x,y,z,3) = u_i(x,y,z,3) + 2._PRC * T0_i * ((2._PRC*M_PI/xsize)**2 + (M_PI/ysize)**2) * psi * &
                            cos(2._PRC*M_PI*xVal/xsize)*cos(M_PI*yVal/ysize) / (n_i(x,y,z)+n_i_bg(x,y,z))
                    end if
                else ! gauss perturbation
                    a = xsize/ysize
                    ! y,z face; x center
                    B(x,y,z,1) = B(x,y,z,1)  - psi * 2._PRC*(yVal-0.5_PRC*dx(2))/(lambda)**2 * &
                        exp(-(xVal/(a*lambda))**2) * exp(-((yVal-0.5_PRC*dx(2))/lambda)**2)
                    ! x,z face; y center
                    B(x,y,z,2) = B(x,y,z,2) + psi * 2._PRC*(xVal-0.5_PRC*dx(1))/(a*lambda)**2 * &
                        exp(-((xVal-0.5_PRC*dx(1))/(a*lambda))**2) * exp(-(yVal/lambda)**2)
                end if

                if (abs(yVal) < xe(2)*.8_PRC) then
                    B(x,y,z,1) = B(x,y,z,1) + noise_level*(rnd-0.5_PRC)
                end if
            enddo
        enddo
    enddo

end subroutine setup_harris_sheet


subroutine setup_double_harris_sheet(n_e,n_i,u_e,u_i,T0_e,T0_i,E,B,lambda,psi,n_bg,guide_field,noise_level,dimX,BD,xb,xe,dx,xb_loc)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in) :: xb(3), xe(3), dx(3), xb_loc(3), T0_e, T0_i, lambda, psi, n_bg, guide_field, noise_level
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n_e, n_i
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u_e, u_i, E, B

    integer :: x, y, z, seed_size
    integer, allocatable :: my_seed(:)
    real(kind=PRC) :: xVal, yVal, zVal, xsize, ysize, rnd

    ! in doi:10.1029/2006GL025957

    ! always use xb = yb = 0 for this setup

    ! rng different on each process
    call random_seed(size=seed_size)
    allocate(my_seed(seed_size))
    call random_seed(get=my_seed)
    my_seed(1:min(3,seed_size)) = nint(xb_loc(1:min(3,seed_size))/dx(1:min(3,seed_size)))
    call random_seed(put=my_seed)

    xsize = xe(1) - xb(1)
    ysize = xe(2) - xb(2)

    do z = -BD(3),dimX(3)+BD(3)
        zVal = xb_loc(3) + (z+0.5_PRC)*dx(3)
        do y = -BD(2),dimX(2)+BD(2)
            yVal = xb_loc(2) + (y+0.5_PRC)*dx(2)
            do x = -BD(1),dimX(1)+BD(1)
                xVal = xb_loc(1) + (x+0.5_PRC)*dx(1)

                call random_number(rnd)

                ! cell center
                n_e(x,y,z) = n_bg + (1._PRC/cosh((yVal-0.25*ysize)/lambda))**2 + &
                                    (1._PRC/cosh((yVal-0.75*ysize)/lambda))**2
                n_i(x,y,z) = n_bg + (1._PRC/cosh((yVal-0.25*ysize)/lambda))**2 + &
                                    (1._PRC/cosh((yVal-0.75*ysize)/lambda))**2

                u_e(x,y,z,1) = 0._PRC
                u_e(x,y,z,2) = 0._PRC
                u_e(x,y,z,3) = 2._PRC * T0_e * (1._PRC/lambda*((1._PRC/cosh((yVal-0.25*ysize)/lambda))**2 - &
                                                               (1._PRC/cosh((yVal-0.75*ysize)/lambda))**2)) / n_e(x,y,z)

                u_i(x,y,z,1) = 0._PRC
                u_i(x,y,z,2) = 0._PRC
                u_i(x,y,z,3) = -2._PRC * T0_i * (1._PRC/lambda*((1._PRC/cosh((yVal-0.25*ysize)/lambda))**2 - &
                                                                (1._PRC/cosh((yVal-0.75*ysize)/lambda))**2)) / n_i(x,y,z)

                ! Yee cube
                ! y,z center; x face
                E(x,y,z,1) = 0._PRC
                ! x,z center; y face
                E(x,y,z,2) = 0._PRC
                ! x,y center; z face
                E(x,y,z,3) = 0._PRC
                ! y,z face; x center
                B(x,y,z,1) = tanh((yVal-0.5_PRC*dx(2)-0.25_PRC*ysize)/lambda) - &
                             tanh((yVal-0.5_PRC*dx(2)-0.75_PRC*ysize)/lambda) - 1._PRC &
                             - psi * M_PI/ysize * ( &
                             cos(2._PRC*M_PI/xsize*xVal) * sin(M_PI/ysize*(yVal-0.5_PRC*dx(2)-0.25_PRC*ysize)) + &
                             cos(2._PRC*M_PI/xsize*xVal) * sin(M_PI/ysize*(yVal-0.5_PRC*dx(2)-0.75_PRC*ysize)) )
                ! x,z face; y center
                B(x,y,z,2) = psi * 2._PRC*M_PI/xsize * ( &
                             sin(2._PRC*M_PI/xsize*(xVal-0.5_PRC*dx(1))) * cos(M_PI/ysize*(yVal-0.25_PRC*ysize)) + &
                             sin(2._PRC*M_PI/xsize*(xVal-0.5_PRC*dx(1))) * cos(M_PI/ysize*(yVal-0.75_PRC*ysize)) )
                ! x,y face; z center
                B(x,y,z,3) = guide_field

                if (yVal > xe(2)*.2_PRC .and. yVal < xe(2)*.8_PRC) then
                    B(x,y,z,1) = B(x,y,z,1) + noise_level*(rnd-0.5)
                end if
            enddo
        enddo
    enddo

end subroutine setup_double_harris_sheet


subroutine setup_double_harris_sheet_gaussian_hump(n_e,n_i,u_e,u_i,n_e_bg,n_i_bg,u_e_bg,u_i_bg,T0_e,T0_i,&
                                                   E,B,n_bg,lambda,guide_field,noise_level,dimX,BD,xb,xe,dx,xb_loc)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in) :: xb(3), xe(3), dx(3), xb_loc(3), T0_e, T0_i, n_bg, lambda, guide_field, noise_level
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: &
        n_e, n_i, n_e_bg, n_i_bg
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: &
        u_e, u_i, u_e_bg, u_i_bg, E, B

    integer :: x, y, z, seed_size
    integer, allocatable :: my_seed(:)
    real(kind=PRC) :: xVal, yVal, zVal, xsize, ysize, rnd, p_x, lambda_x, lambda_y

    ! EMfields3D::initDoublePeriodicHarrisWithGaussianHumpPerturbation in iPic3D GitHub
    ! TODO: non-drifting background
    ! xsize = 30 d_i0, ysize = 40 d_i0

    ! rng different on each process
    call random_seed(size=seed_size)
    allocate(my_seed(seed_size))
    call random_seed(get=my_seed)
    my_seed(1:min(3,seed_size)) = nint(xb_loc(1:min(3,seed_size))/dx(1:min(3,seed_size)))
    call random_seed(put=my_seed)

    xsize = xe(1) - xb(1)
    ysize = xe(2) - xb(2)
    p_x = 0.4_PRC
    lambda_x = 8._PRC*lambda
    lambda_y = 4._PRC*lambda

    do z = -BD(3),dimX(3)+BD(3)
        zVal = xb_loc(3) + (z+0.5_PRC)*dx(3)
        do y = -BD(2),dimX(2)+BD(2)
            yVal = xb_loc(2) + (y+0.5_PRC)*dx(2)
            do x = -BD(1),dimX(1)+BD(1)
                xVal = xb_loc(1) + (x+0.5_PRC)*dx(1)

                call random_number(rnd)

                ! cell center
                n_e(x,y,z) = (1._PRC/cosh((yVal-0.25_PRC*ysize)/lambda))**2 + &
                             (1._PRC/cosh((yVal-0.75_PRC*ysize)/lambda))**2
                n_i(x,y,z) = (1._PRC/cosh((yVal-0.25_PRC*ysize)/lambda))**2 + &
                             (1._PRC/cosh((yVal-0.75_PRC*ysize)/lambda))**2
                n_e_bg(x,y,z) = n_bg
                n_i_bg(x,y,z) = n_bg

                u_e(x,y,z,1) = 0._PRC
                u_e(x,y,z,2) = 0._PRC
                u_e(x,y,z,3) = 2._PRC*T0_e/lambda
                u_i(x,y,z,1) = 0._PRC
                u_i(x,y,z,2) = 0._PRC
                u_i(x,y,z,3) = -2._PRC*T0_i/lambda

                u_e_bg(x,y,z,1) = 0._PRC
                u_e_bg(x,y,z,2) = 0._PRC
                u_e_bg(x,y,z,3) = 0._PRC
                u_i_bg(x,y,z,1) = 0._PRC
                u_i_bg(x,y,z,2) = 0._PRC
                u_i_bg(x,y,z,3) = 0._PRC

                ! Yee cube
                ! y,z center; x face
                E(x,y,z,1) = 0._PRC
                ! x,z center; y face
                E(x,y,z,2) = 0._PRC
                ! x,y center; z face
                E(x,y,z,3) = 0._PRC
                ! y,z face; x center
                B(x,y,z,1) = tanh((yVal-0.5_PRC*dx(2)-0.25_PRC*ysize)/lambda) - &
                             tanh((yVal-0.5_PRC*dx(2)-0.75_PRC*ysize)/lambda) - 1._PRC + &
                             2._PRC*p_x*((-(yVal-0.5_PRC*dx(2)-0.25_PRC*ysize)/lambda_y) * &
                             exp(-((xVal-0.25_PRC*xsize)/lambda_x)**2 - &
                                  ((yVal-0.5_PRC*dx(2)-0.25_PRC*ysize)/lambda_y)**2) + &
                             (yVal-0.5_PRC*dx(2)-0.75_PRC*ysize)/lambda_y * &
                             exp(-((xVal-0.75_PRC*xsize)/lambda_x)**2 - &
                                  ((yVal-0.5_PRC*dx(2)-0.75_PRC*ysize)/lambda_y)**2))
                ! x,z face; y center
                B(x,y,z,2) = 2._PRC*p_x*((xVal-0.5_PRC*dx(1)-0.25_PRC*xsize)/lambda_x * &
                             exp(-((xVal-0.25_PRC*xsize)/lambda_x)**2 - &
                                  ((yVal-0.5_PRC*dx(2)-0.25_PRC*ysize)/lambda_y)**2) - &
                             (xVal-0.5_PRC*dx(1)-0.75_PRC*xsize)/lambda_x * &
                             exp(-((xVal-0.75_PRC*xsize)/lambda_x)**2 - &
                                  ((yVal-0.5_PRC*dx(2)-0.75_PRC*ysize)/lambda_y)**2))
                ! x,y face; z center
                B(x,y,z,3) = guide_field

                if (yVal > xe(2)*.2_PRC .and. yVal < xe(2)*.8_PRC) then
                    B(x,y,z,1) = B(x,y,z,1) + noise_level*(rnd-0.5)
                end if
            enddo
        enddo
    enddo

end subroutine setup_double_harris_sheet_gaussian_hump


subroutine setup_island_coalescence(n_e,n_i,u_e,u_i,T0_e,T0_i,E,B,lambda,eta,psi,n_bg,guide_field,noise_level,dimX,BD,xe,dx,xb_loc)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in) :: xe(3), dx(3), xb_loc(3), T0_e, T0_i, lambda, eta, psi, n_bg, guide_field, noise_level
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n_e, n_i
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u_e, u_i, E, B

    integer :: x, y, z, seed_size
    integer, allocatable :: my_seed(:)
    real(kind=PRC) :: xVal, yVal, zVal, rnd

    ! in Ng et al. 2015
    ! https://doi.org/10.1063/1.4935302

    ! rng different on each process
    call random_seed(size=seed_size)
    allocate(my_seed(seed_size))
    call random_seed(get=my_seed)
    my_seed(1:min(3,seed_size)) = nint(xb_loc(1:min(3,seed_size))/dx(1:min(3,seed_size)))
    call random_seed(put=my_seed)

    do z = -BD(3),dimX(3)+BD(3)
        zVal = xb_loc(3) + (z+0.5_PRC)*dx(3)
        do y = -BD(2),dimX(2)+BD(2)
            yVal = xb_loc(2) + (y+0.5_PRC)*dx(2)
            do x = -BD(1),dimX(1)+BD(1)
                xVal = xb_loc(1) + (x+0.5_PRC)*dx(1)

                call random_number(rnd)

                ! cell center
                n_e(x,y,z) = n_bg + (1._PRC - eta*eta) / (cosh(xVal/lambda) + eta*cos(yVal/lambda))**2
                n_i(x,y,z) = n_bg + (1._PRC - eta*eta) / (cosh(xVal/lambda) + eta*cos(yVal/lambda))**2

                if (abs(xVal) < xe(1)*.8_PRC) then
                    n_e(x,y,z) = n_e(x,y,z) + noise_level*(rnd-0.5)
                    n_i(x,y,z) = n_i(x,y,z) + noise_level*(rnd-0.5)
                end if

                u_e(x,y,z,1) = 0._PRC
                u_e(x,y,z,2) = 0._PRC
                u_e(x,y,z,3) = - 2._PRC*T0_e / lambda * ((1._PRC - eta*eta) / (cosh(xVal/lambda) + &
                            eta*cos(yVal/lambda))**2) / n_e(x,y,z) ! vlasov equilibrium
                u_i(x,y,z,1) = 0._PRC
                u_i(x,y,z,2) = 0._PRC
                u_i(x,y,z,3) = 2._PRC*T0_i / lambda * ((1._PRC - eta*eta) / (cosh(xVal/lambda) + &
                            eta*cos(yVal/lambda))**2) / n_i(x,y,z) ! vlasov equilibrium

                ! Yee cube
                ! y,z center; x face
                E(x,y,z,1) = 0._PRC
                ! x,z center; y face
                E(x,y,z,2) = 0._PRC
                ! x,y center; z face
                E(x,y,z,3) = 0._PRC
                ! y,z face; x center
                B(x,y,z,1) = eta*sin((yVal-0.5_PRC*dx(2))/lambda) / (cosh(xVal/lambda) + &
                            eta*cos((yVal-0.5_PRC*dx(2))/lambda)) + &
                            psi * sin((yVal-0.5_PRC*dx(2))/(2._PRC*lambda) - M_PI) * cos(xVal/(2._PRC*lambda))
                ! x,z face; y center
                B(x,y,z,2) = (exp((xVal-0.5_PRC*dx(1))/lambda)-cosh((xVal-0.5_PRC*dx(1))/lambda)) / &
                            ( cosh((xVal-0.5_PRC*dx(1))/lambda) + eta*cos(yVal/lambda)) &
                            - psi * cos(yVal/(2._PRC*lambda) - M_PI) * sin((xVal-0.5_PRC*dx(1))/(2._PRC*lambda))
                ! x,y face; z center
                B(x,y,z,3) = guide_field
            enddo
        enddo
    enddo

end subroutine setup_island_coalescence


subroutine setup_orszag_tang(n_e,n_i,u_e,u_i,E,B,c,L,n_bg,delta_u,dimX,BD,dx,xb_loc)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in) :: dx(3), xb_loc(3), c, L, n_bg, delta_u
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n_e, n_i
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u_e, u_i, E, B

    real(kind=PRC) :: delta_b
    integer :: x,y,z
    real(kind=PRC) :: xVal,yVal,zVal

    ! in Groselj et al. 2017
    ! https://doi.org/10.3847%2F1538-4357%2Faa894d

    delta_b = delta_u

    do z = -BD(3),dimX(3)+BD(3)
        zVal = xb_loc(3) + (z+0.5_PRC)*dx(3)
        do y = -BD(2),dimX(2)+BD(2)
            yVal = xb_loc(2) + (y+0.5_PRC)*dx(2)
            do x = -BD(1),dimX(1)+BD(1)
                xVal = xb_loc(1) + (x+0.5_PRC)*dx(1)

                ! cell center
                n_e(x,y,z) = n_bg + 2._PRC*M_PI/L * delta_u * &
                    (cos(xVal*2._PRC*M_PI/L) + cos(yVal*2._PRC*M_PI/L)) / (c*c)
                n_i(x,y,z) = n_bg

                u_e(x,y,z,1) = -delta_u * sin(yVal*2._PRC*M_PI/L)
                u_e(x,y,z,2) =  delta_u * sin(xVal*2._PRC*M_PI/L)
                u_e(x,y,z,3) = - 2._PRC*M_PI/L * delta_b * &
                    (2._PRC*cos(xVal*4._PRC*M_PI/L) + cos(yVal*2._PRC*M_PI/L)) ! aforementioned paper
                !u_e(x,y,z,3) = -m_i/(m_e+m_i) * 2._PRC*M_PI/L * delta_b * &
                !    (2._PRC*cos(xVal*4._PRC*M_PI/L) + cos(yVal*2._PRC*M_PI/L)) / n_e(x,y,z) ! alternative

                u_i(x,y,z,1) = -delta_u * sin(yVal*2._PRC*M_PI/L)
                u_i(x,y,z,2) =  delta_u * sin(xVal*2._PRC*M_PI/L)
                u_i(x,y,z,3) =  0._PRC ! aforementioned paper
                !u_i(x,y,z,3) =  m_e/(m_e+m_i) * 2._PRC*M_PI/L * delta_b * &
                !   (2._PRC*cos(xVal*4._PRC*M_PI/L) + cos(yVal*2._PRC*M_PI/L)) / n_i(x,y,z) ! alternative

                ! Yee cube
                ! y,z center; x face
                E(x,y,z,1) = -delta_u * sin((xVal-0.5_PRC*dx(1))*2._PRC*M_PI/L)
                ! x,z center; y face
                E(x,y,z,2) = -delta_u * sin((yVal-0.5_PRC*dx(2))*2._PRC*M_PI/L)
                ! x,y center; z face
                E(x,y,z,3) = 0._PRC
                ! y,z face; x center
                B(x,y,z,1) = -delta_b * sin((yVal-0.5_PRC*dx(2))*2._PRC*M_PI/L)
                ! x,z face; y center
                B(x,y,z,2) =  delta_b * sin((xVal-0.5_PRC*dx(1))*4._PRC*M_PI/L)
                ! x,y face; z center
                B(x,y,z,3) = 1._PRC

            enddo
        enddo
    enddo

end subroutine setup_orszag_tang


subroutine setup_orszag_tang_3d(n_e,n_i,u_e,u_i,E,B,L,L_z,n_bg,delta_u,cross_helicity,dimX,BD,dx,xb_loc)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in) :: dx(3), xb_loc(3), L, L_z, n_bg, delta_u, cross_helicity
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n_e, n_i
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u_e, u_i, E, B

    real(kind=PRC) :: delta_b, k_perp, k_par, phases(4), sigma
    integer :: x,y,z
    real(kind=PRC) :: xVal,yVal,zVal

    ! in Li et al. 2016
    ! https://arxiv.org/abs/1510.02842

    ! in dissertation D. Groselj 2018
    ! Fully Kinetic Simulations of Microscale Turbulence in Space and Astrophysical Plasmas
    ! and PRL 2018

    delta_b = delta_u

    if (abs(cross_helicity) > 0.5) then
        write(*,*) "Warning: abs(cross_helicity)>0.5 in Orszag-Tang 3D setup not allowed. Falling back to abs(cross_helicity)=0.5."
    end if
    sigma = cross_helicity
    if (sigma > 0.5_PRC) then
        sigma = 0.5_PRC
    else if (sigma < -0.5_PRC) then
        sigma = -0.5_PRC
    end if

    phases(1) = 0.1_PRC
    phases(2) = phases(1) - acos(2._PRC*sigma)
    phases(3) = 0.7_PRC
    phases(4) = 1.1_PRC

    k_perp = 2._PRC*M_PI/L
    k_par = 2._PRC*M_PI/L_z

    do z = -BD(3),dimX(3)+BD(3)
        zVal = xb_loc(3) + (z+0.5_PRC)*dx(3)
        do y = -BD(2),dimX(2)+BD(2)
            yVal = xb_loc(2) + (y+0.5_PRC)*dx(2)
            do x = -BD(1),dimX(1)+BD(1)
                xVal = xb_loc(1) + (x+0.5_PRC)*dx(1)

                ! cell center
                n_e(x,y,z) = n_bg
                n_i(x,y,z) = n_bg

                u_e(x,y,z,1) = delta_u * (cos(k_perp*yVal + phases(1))*sin(k_par*zVal) - &
                                          sin(k_perp*yVal + phases(2))*cos(k_par*zVal)) - &
                               k_par * delta_b * (cos(k_perp*xVal + phases(3))*cos(k_par*zVal) + &
                                                  sin(2._PRC*k_perp*xVal + phases(4))*sin(k_par*zVal))
                u_e(x,y,z,2) = delta_u * (sin(k_perp*xVal + phases(3))*cos(k_par*zVal) - &
                                          cos(2._PRC*k_perp*xVal + phases(4))*sin(k_par*zVal)) - &
                               k_par * delta_b * (cos(k_perp*yVal + phases(2))*cos(k_par*zVal) + &
                                                  sin(k_perp*yVal + phases(1))*sin(k_par*zVal))

                u_e(x,y,z,3) = -k_perp * delta_b * (&
                    (sin(k_perp*xVal + phases(3)) + sin(k_perp*yVal + phases(2)))*sin(k_par*zVal) + &
                    (2._PRC*cos(2._PRC*k_perp*xVal + phases(4)) + cos(k_perp*yVal) + phases(1))*cos(k_par*zVal))

                u_i(x,y,z,1) = delta_u * (cos(k_perp*yVal + phases(1))*sin(k_par*zVal) - &
                                          sin(k_perp*yVal + phases(2))*cos(k_par*zVal))
                u_i(x,y,z,2) = delta_u * (sin(k_perp*xVal + phases(3))*cos(k_par*zVal) - &
                                          cos(2._PRC*k_perp*xVal + phases(4))*sin(k_par*zVal))
                u_i(x,y,z,3) =  0._PRC

                ! Yee cube
                ! y,z center; x face
                E(x,y,z,1) = 0._PRC
                ! x,z center; y face
                E(x,y,z,2) = 0._PRC
                ! x,y center; z face
                E(x,y,z,3) = 0._PRC
                ! y,z face; x center
                B(x,y,z,1) = delta_b * (cos(k_perp*(yVal-0.5_PRC*dx(2)) + phases(1))*sin(k_par*(zVal-0.5_PRC*dx(3))) - &
                                        sin(k_perp*(yVal-0.5_PRC*dx(2)) + phases(2))*cos(k_par*(zVal-0.5_PRC*dx(3))))
                ! x,z face; y center
                B(x,y,z,2) = delta_b * (sin(k_perp*(xVal-0.5_PRC*dx(1)) + phases(3))*cos(k_par*(zVal-0.5_PRC*dx(3))) - &
                                        cos(2._PRC*k_perp*(xVal-0.5_PRC*dx(1)) + phases(4))*sin(k_par*(zVal-0.5_PRC*dx(3))))
                ! x,y face; z center
                B(x,y,z,3) = 1._PRC

            enddo
        enddo
    enddo

end subroutine setup_orszag_tang_3d


subroutine setup_cssi(n_e,n_i,u_e,u_i,T_e,T_i,E,B,n_bg,T_bg,delta_e,delta_i,B0_e,B0_i,V0_e,V0_i,noise_level,dimX,BD,xb,xe,dx,xb_loc)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in) :: xb(3), xe(3), dx(3), xb_loc(3), n_bg, T_bg, delta_e, delta_i, B0_e, B0_i, V0_e, V0_i, noise_level
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n_e, n_i, T_e, T_i
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u_e, u_i, E, B

    integer :: x, y, z, seed_size
    integer, allocatable :: my_seed(:)
    real(kind=PRC) :: xVal, yVal, zVal, yVal_global, rnd

    ! current sheet shear instability (CSSI)
    ! Fujimoto & Sydora (2017): Linear theory of the current sheet shear instability
    ! Fujimoto & Sydora (2012): Plasmoid-Induced Turbulence in Collisionless Magnetic Reconnection

    ! coordinate system:
    ! Fujimoto & Sydora 2017 -> setup_cssi: x->z, y->x, z->y

    ! rng different on each process
    call random_seed(size=seed_size)
    allocate(my_seed(seed_size))
    call random_seed(get=my_seed)
    my_seed(1:min(3,seed_size)) = nint(xb_loc(1:min(3,seed_size))/dx(1:min(3,seed_size)))
    call random_seed(put=my_seed)

    do z = -BD(3),dimX(3)+BD(3)
        zVal = xb_loc(3) + (z+0.5_PRC)*dx(3)
        do y = -BD(2),dimX(2)+BD(2)
            yVal = xb_loc(2) + (y+0.5_PRC)*dx(2)
            do x = -BD(1),dimX(1)+BD(1)
                xVal = xb_loc(1) + (x+0.5_PRC)*dx(1)

                call random_number(rnd)

                ! cell center
                n_e(x,y,z) = n_bg
                n_i(x,y,z) = n_bg

                if (abs(yVal) < xe(2)*.8_PRC) then
                    n_e(x,y,z) = n_e(x,y,z) + noise_level*(rnd-0.5)
                    n_i(x,y,z) = n_i(x,y,z) + noise_level*(rnd-0.5)
                end if

                u_e(x,y,z,1) = -V0_e / (cosh(yVal/delta_e))**2
                u_e(x,y,z,2) = 0._PRC
                u_e(x,y,z,3) = 0._PRC
                u_i(x,y,z,1) = -V0_i / (cosh(yVal/delta_i))**2
                u_i(x,y,z,2) = 0._PRC
                u_i(x,y,z,3) = 0._PRC
                
                ! Yee cube
                ! y,z center; x face
                E(x,y,z,1) = 0._PRC
                ! x,z center; y face
                E(x,y,z,2) = 0._PRC
                ! x,y center; z face
                E(x,y,z,3) = 0._PRC
                ! y,z face; x center
                B(x,y,z,1) = 0._PRC
                ! x,z face; y center
                B(x,y,z,2) = 0._PRC
                ! x,y face; z center
                B(x,y,z,3) = -B0_i*tanh((yVal-0.5_PRC*dx(2))/delta_i) - B0_e*tanh((yVal-0.5_PRC*dx(2))/delta_e)

                ! cell center
                T_e(x,y,z) = 0._PRC
                T_i(x,y,z) = 0._PRC
                yVal_global = xb(2)!+0.5_PRC*dx(2)
                do while (yVal_global <= xb_loc(2)+(y+0.5_PRC)*dx(2))
                    T_e(x,y,z) = T_e(x,y,z) - (V0_e/(cosh(yVal_global/delta_e))**2 * &
                                 (-B0_i*tanh(yVal_global/delta_i) - B0_e*tanh(yVal_global/delta_e)))*dx(2)
                    T_i(x,y,z) = T_i(x,y,z) + (V0_i/(cosh(yVal_global/delta_i))**2 * &
                                 (-B0_i*tanh(yVal_global/delta_i) - B0_e*tanh(yVal_global/delta_e)))*dx(2)
                    yVal_global = yVal_global + dx(2)
                end do

                T_e(x,y,z) = T_e(x,y,z)/n_e(x,y,z) + T_bg
                T_i(x,y,z) = T_i(x,y,z)/n_i(x,y,z) + T_bg
            enddo
        enddo
    enddo

end subroutine setup_cssi


subroutine setup_cssi_harris(n_e,n_i,u_e,u_i,n_e_bg,n_i_bg,u_e_bg,u_i_bg,T0_e,T0_i,E,B,&
                             n_harris,n_bg,lambda,noise_level,dimX,BD,xe,dx,xb_loc)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in) :: xe(3), dx(3), xb_loc(3), T0_e, T0_i, n_harris, n_bg, lambda, noise_level
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n_e, n_i, n_e_bg, n_i_bg
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u_e, u_i, u_e_bg, u_i_bg, E, B

    integer :: x, y, z, seed_size
    integer, allocatable :: my_seed(:)
    real(kind=PRC) :: xVal, yVal, zVal, rnd

    ! current sheet shear instability (CSSI)
    ! Fujimoto & Sydora (2017): Linear theory of the current sheet shear instability
    ! Fujimoto & Sydora (2012): Plasmoid-Induced Turbulence in Collisionless Magnetic Reconnection

    ! rng different on each process
    call random_seed(size=seed_size)
    allocate(my_seed(seed_size))
    call random_seed(get=my_seed)
    my_seed(1:min(3,seed_size)) = nint(xb_loc(1:min(3,seed_size))/dx(1:min(3,seed_size)))
    call random_seed(put=my_seed)

    do z = -BD(3),dimX(3)+BD(3)
        zVal = xb_loc(3) + (z+0.5_PRC)*dx(3)
        do y = -BD(2),dimX(2)+BD(2)
            yVal = xb_loc(2) + (y+0.5_PRC)*dx(2)
            do x = -BD(1),dimX(1)+BD(1)
                xVal = xb_loc(1) + (x+0.5_PRC)*dx(1)

                call random_number(rnd)

                ! cell center
                n_e(x,y,z) = n_harris*(1._PRC/cosh(yVal/lambda))**2
                n_i(x,y,z) = n_harris*(1._PRC/cosh(yVal/lambda))**2
                u_e(x,y,z,1) = 2._PRC*T0_e/lambda / n_harris
                u_e(x,y,z,2) = 0._PRC
                u_e(x,y,z,3) = 0._PRC
                u_i(x,y,z,1) = -2._PRC*T0_i/lambda / n_harris
                u_i(x,y,z,2) = 0._PRC
                u_i(x,y,z,3) = 0._PRC
               
                n_e_bg(x,y,z) = n_bg
                n_i_bg(x,y,z) = n_bg
                u_e_bg(x,y,z,1) = 0._PRC
                u_e_bg(x,y,z,2) = 0._PRC
                u_e_bg(x,y,z,3) = 0._PRC
                u_i_bg(x,y,z,1) = 0._PRC
                u_i_bg(x,y,z,2) = 0._PRC
                u_i_bg(x,y,z,3) = 0._PRC
 
                ! Yee cube
                ! y,z center; x face
                E(x,y,z,1) = 0._PRC
                ! x,z center; y face
                E(x,y,z,2) = 0._PRC
                ! x,y center; z face
                E(x,y,z,3) = 0._PRC
                ! y,z face; x center
                B(x,y,z,1) = 0._PRC
                ! x,z face; y center
                B(x,y,z,2) = 0._PRC
                ! x,y face; z center
                B(x,y,z,3) = -tanh((yVal-0.5_PRC*dx(2))/lambda)

                if (abs(yVal) < xe(2)*.9_PRC) then
                    !n_e(x,y,z) = n_e(x,y,z) + noise_level*(rnd-0.5_PRC)
                    !n_i(x,y,z) = n_i(x,y,z) + noise_level*(rnd-0.5_PRC)
                    B(x,y,z,3) = B(x,y,z,3) + noise_level*(rnd-0.5_PRC)
                end if
            enddo
        enddo
    enddo

end subroutine setup_cssi_harris


subroutine setup_lhdi(n_e,n_i,u_e,u_i,T0_e,T0_i,E,B,n_bg,lambda,noise_level,dimX,BD,xe,dx,xb_loc)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in) :: xe(3), dx(3), xb_loc(3), T0_e, T0_i, n_bg, lambda, noise_level
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n_e, n_i
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u_e, u_i, E, B

    integer :: x, y, z, seed_size
    integer, allocatable :: my_seed(:)
    real(kind=PRC) :: xVal, yVal, zVal, rnd

    ! lower hybrid drift instability (LHDI)

    ! in Ng et al. (2019)
    ! https://doi.org/10.1029/2018JA026313
    ! https://arxiv.org/abs/1903.09618

    ! coordinate system:
    ! Ng 2019 -> setup_lhdi: x->z, y->x, z->y

    ! rng different on each process
    call random_seed(size=seed_size)
    allocate(my_seed(seed_size))
    call random_seed(get=my_seed)
    my_seed(1:min(3,seed_size)) = nint(xb_loc(1:min(3,seed_size))/dx(1:min(3,seed_size)))
    call random_seed(put=my_seed)

    do z = -BD(3),dimX(3)+BD(3)
        zVal = xb_loc(3) + (z+0.5_PRC)*dx(3)
        do y = -BD(2),dimX(2)+BD(2)
            yVal = xb_loc(2) + (y+0.5_PRC)*dx(2)
            do x = -BD(1),dimX(1)+BD(1)
                xVal = xb_loc(1) + (x+0.5_PRC)*dx(1)

                call random_number(rnd)

                ! cell center
                n_e(x,y,z) = n_bg + (1._PRC/cosh(yVal/lambda))**2
                n_i(x,y,z) = n_bg + (1._PRC/cosh(yVal/lambda))**2

                u_e(x,y,z,1) = 2._PRC*T0_e/lambda*(1._PRC/cosh(yVal/lambda))**2 / n_e(x,y,z)
                u_e(x,y,z,2) = 0._PRC
                u_e(x,y,z,3) = 0._PRC
                u_i(x,y,z,1) = -2._PRC*T0_i/lambda*(1._PRC/cosh(yVal/lambda))**2 / n_i(x,y,z)
                u_i(x,y,z,2) = 0._PRC
                u_i(x,y,z,3) = 0._PRC
               
                ! Yee cube
                ! y,z center; x face
                E(x,y,z,1) = 0._PRC
                ! x,z center; y face
                E(x,y,z,2) = 0._PRC
                ! x,y center; z face
                E(x,y,z,3) = 0._PRC
                ! y,z face; x center
                B(x,y,z,1) = 0._PRC
                ! x,z face; y center
                B(x,y,z,2) = 0._PRC
                ! x,y face; z center
                B(x,y,z,3) = -tanh((yVal-0.5_PRC*dx(2))/lambda)

                if (abs(yVal) < xe(2)*.9_PRC) then
                    !n_e(x,y,z) = n_e(x,y,z) + noise_level*(rnd-0.5_PRC)
                    !n_i(x,y,z) = n_i(x,y,z) + noise_level*(rnd-0.5_PRC)
                    B(x,y,z,3) = B(x,y,z,3) + noise_level*(rnd-0.5_PRC)
                end if
            enddo
        enddo
    enddo

end subroutine setup_lhdi


subroutine setup_khi(n_e,n_i,u_e,u_i,T_e,T_i,E,B,lambda,T0e_T0i,wpe_wce,m_i,dimX,BD,xe,dx,xb_loc)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in) :: xe(3), dx(3), xb_loc(3), lambda, T0e_T0i, wpe_wce, m_i
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n_e, n_i, T_e, T_i
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u_e, u_i, E, B

    integer :: x, y, z, m, dawndusk
    real(kind=PRC) :: xVal, yVal, zVal, du, n_r, B_r, T_r, pbc, jy, dxEx, phase(10), B0, v0

    ! Kelvin-Helmholtz instability
    
    B0 = 1._PRC  / wpe_wce ! 1/(omega_p,e / Omega_c,e)
    v0 = 1._PRC / (sqrt(m_i)*wpe_wce) ! 1/(sqrt(m_i/m_e)*(omega_p,e / Omega_c,e))
    du = 5._PRC ! velocity shear

    ! sheath/sphere ratios
    n_r = 2._PRC
    B_r = 0.2_PRC
    T_r = 4._PRC

    dawndusk = 1

    pbc = B0**2*(0.5_PRC + 0.5_PRC*(1._PRC-B_r**2)/(T_r*n_r-1._PRC) * n_r*T_r/B_r**2) ! pressure balance constant

    ! phases of perturbation modes
    phase = [0.32901019493129924_PRC, 0.5478100511920081_PRC, 2.5587754659637483_PRC, 0.6767005337449898_PRC, &
             5.662399548792817_PRC, 0.23972655426914796_PRC, 3.3690567796208897_PRC, 2.0872596983548943_PRC, &
             5.35381812450134_PRC, 1.0031884251925596_PRC]

    do z = -BD(3),dimX(3)+BD(3)
        zVal = xb_loc(3) + (z+0.5_PRC)*dx(3)
        do y = -BD(2),dimX(2)+BD(2)
            yVal = xb_loc(2) + (y+0.5_PRC)*dx(2)
            do x = -BD(1),dimX(1)+BD(1)
                xVal = xb_loc(1) + (x+0.5_PRC)*dx(1)

                ! cell center
                u_i(x,y,z,1) = 0._PRC
                u_i(x,y,z,2) = v0*du/2._PRC*(1._PRC + dawndusk*tanh(xVal/lambda))
                u_i(x,y,z,3) = 0._PRC

                jy = ( -dawndusk*B0*(B_r-1._PRC)/(2._PRC*lambda*B_r)*(1._PRC - tanh(xVal/lambda)**2) )
                dxEx = -(v0*du/2._PRC*dawndusk*(1._PRC -tanh(xVal/lambda)**2)/lambda * &
                         B0/B_r*(1._PRC + (B_r-1._PRC)/2._PRC*(1._PRC + dawndusk*tanh(xVal/lambda))) + & ! Bz
                         u_i(x,y,z,2)*dawndusk*B0*(B_r-1._PRC)/(2._PRC*lambda*B_r) * (1._PRC - tanh(xVal/lambda)**2) )
               
                n_i(x,y,z) = 1._PRC/n_r*(1._PRC + (n_r-1._PRC)/2._PRC*(1._PRC + dawndusk*tanh(xVal/lambda)))
                n_e(x,y,z) = n_i(x,y,z) - dxEx

                u_e(x,y,z,1) = 0._PRC
                u_e(x,y,z,2) = (n_i(x,y,z)/n_e(x,y,z)*u_i(x,y,z,2) - jy/n_e(x,y,z) )
                u_e(x,y,z,3) = 0._PRC

                T_i(x,y,z) = (pbc - (B0/B_r*(1._PRC + (B_r-1._PRC)/2._PRC*(1._PRC + dawndusk*tanh(xVal/lambda))))**2/2._PRC) & ! (pbc - Bz**2/2) /
                              /(n_i(x,y,z)+T0e_T0i*n_e(x,y,z))
                T_e(x,y,z) = T_i(x,y,z) * T0e_T0i

                ! Yee cube
                ! y,z face; x center
                B(x,y,z,1) = 0._PRC
                ! x,z face; y center
                B(x,y,z,2) = 0._PRC
                ! x,y face; z center
                B(x,y,z,3) = B0/B_r*(1._PRC + (B_r-1._PRC)/2._PRC*(1._PRC + dawndusk*tanh((xVal-0.5_PRC*dx(1))/lambda)))
                ! y,z center; x face
                E(x,y,z,1) = -v0*du/2._PRC*(1._PRC + dawndusk*tanh((xVal-0.5_PRC*dx(1))/lambda))*(B(x,y,z,3))
                ! x,z center; y face
                E(x,y,z,2) = 0._PRC
                ! x,y center; z face
                E(x,y,z,3) = 0._PRC
                
                do m = 1,10
                    u_i(x,y,z,1) = u_i(x,y,z,1) + v0*exp(-(xVal/lambda)**2)*sin(2._PRC*M_PI*m*yVal/xe(2) + phase(m))*(2._PRC*M_PI/xe(2))
                    u_i(x,y,z,2) = u_i(x,y,z,2) - v0*2._PRC*xVal/lambda**2 * exp(-(xVal/lambda)**2)*cos(2._PRC*M_PI*m*yVal/xe(2) + phase(m))/m
                end do

            enddo
        enddo
    enddo

end subroutine setup_khi



subroutine setup_whistler_wave(n_e,n_i,u_e,u_i,E,B,dimX,BD,dx,xb_loc)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in) :: dx(3), xb_loc(3)
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n_e, n_i
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u_e, u_i, E, B

    integer :: x, y, z
    real(kind=PRC) :: xVal, yVal, zVal
    real(kind=PRC) :: kw, cw, db, du, cosa, sina, beta

    ! in dissertation tht

    kw = 0.62831853_PRC
    cw = 1.169_PRC
    db = 0.05_PRC
    du = 0.05845_PRC
    cosa = 0.44721_PRC
    sina = 0.89443_PRC
    beta = 0.9978_PRC

    n_e = 1._PRC
    n_i = 1._PRC

    do z = -BD(3),dimX(3)+BD(3)
        zVal = xb_loc(3) + (z+0.5_PRC)*dx(3)
        do y = -BD(2),dimX(2)+BD(2)
            yVal = xb_loc(2) + (y+0.5_PRC)*dx(2)
            do x = -BD(1),dimX(1)+BD(1)
                xVal = xb_loc(1) + (x+0.5_PRC)*dx(1)

                ! cell center
                u_e(x,y,z,1) = -sina*du*cos(kw*(cosa*xVal + sina*yVal))
                u_e(x,y,z,2) = cosa*du*cos(kw*(cosa*xVal + sina*yVal))
                u_e(x,y,z,3) = -du*sin(kw*(cosa*xVal + sina*yVal))
                u_i(x,y,z,1) = -(1._PRC-kw*beta/cw)*sina*du*cos(kw*(cosa*xVal + sina*yVal))
                u_i(x,y,z,2) = (1._PRC-kw*beta/cw)*cosa*du*cos(kw*(cosa*xVal + sina*yVal))
                u_i(x,y,z,3) = -du*(1._PRC-kw*beta/cw)*sin(kw*(cosa*xVal + sina*yVal))


                ! Yee cube
                ! y,z center; x face
                E(x,y,z,1) = -sina*cw*db*sin(kw*(cosa*(xVal-0.5_PRC*dx(1)) + sina*yVal))
                ! x,z center; y face
                E(x,y,z,2) = cosa*cw*db*sin(kw*(cosa*xVal + sina*(yVal-0.5_PRC*dx(2))))
                ! x,y center; z face
                E(x,y,z,3) = cw*db*cos(kw*(cosa*xVal + sina*yVal))
                ! y,z face; x center
                B(x,y,z,1) = cosa + sina*db*cos(kw*(cosa*xVal + sina*(yVal-0.5_PRC*dx(2))))
                ! x,z face; y center
                B(x,y,z,2) = -cosa*db*cos(kw*(cosa*(xVal-0.5_PRC*dx(1)) + sina*yVal)) + sina
                ! x,y face; z center
                B(x,y,z,3) = db*sin(kw*(cosa*(xVal-0.5_PRC*dx(1)) + sina*(yVal-0.5_PRC*dx(2))))

            enddo
        enddo
    enddo

end subroutine setup_whistler_wave


subroutine setup_electromagnetic_vacuum_wave(n_e,n_i,u_e,u_i,E,B,dimX,BD,xb,xe,dx,xb_loc)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in) :: xb(3), xe(3), dx(3), xb_loc(3)
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n_e, n_i
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u_e, u_i, E, B

    integer :: x, y, z
    real(kind=PRC) :: xVal, yVal, zVal

    n_e = 1.e-15
    n_i = 1.e-15
    u_e = 0._PRC
    u_i = 0._PRC

    do z = -BD(3),dimX(3)+BD(3)
        zVal = xb_loc(3) + (z+0.5_PRC)*dx(3)
        do y = -BD(2),dimX(2)+BD(2)
            yVal = xb_loc(2) + (y+0.5_PRC)*dx(2)
            do x = -BD(1),dimX(1)+BD(1)
                xVal = xb_loc(1) + (x+0.5_PRC)*dx(1)

                ! Yee cube
                ! y,z center; x face
                E(x,y,z,1) = 0._PRC
                ! x,z center; y face
                E(x,y,z,2) = sin(xVal*2._PRC*M_PI/(xe(1)-xb(1))*2._PRC)
                ! x,y center; z face
                E(x,y,z,3) = 0._PRC
                ! y,z face; x center
                B(x,y,z,1) = 0._PRC
                ! x,z face; y center
                B(x,y,z,2) = 0._PRC
                ! x,y face; z center
                B(x,y,z,3) = sin((xVal-0.5_PRC*dx(1))*2._PRC*M_PI/(xe(1)-xb(1))*2._PRC)

            enddo
        enddo
    enddo

end subroutine setup_electromagnetic_vacuum_wave


subroutine setup_electromagnetic_vacuum_wave_3d_plane_polarized(n_e,n_i,u_e,u_i,E,B,dimX,BD,dx,xb_loc)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in) :: dx(3), xb_loc(3)
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n_e, n_i
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u_e, u_i, E, B

    integer :: x, y, z
    real(kind=PRC) :: xVal, yVal, zVal

    ! in Balsara et al. 2016
    ! https://arxiv.org/abs/1603.06975

    n_e = 1.e-15
    n_i = 1.e-15
    u_e = 0._PRC
    u_i = 0._PRC

    do z = -BD(3),dimX(3)+BD(3)
        zVal = xb_loc(3) + (z+0.5_PRC)*dx(3)
        do y = -BD(2),dimX(2)+BD(2)
            yVal = xb_loc(2) + (y+0.5_PRC)*dx(2)
            do x = -BD(1),dimX(1)+BD(1)
                xVal = xb_loc(1) + (x+0.5_PRC)*dx(1)

                ! Yee cube
                ! y,z center; x face
                E(x,y,z,1) = 0._PRC
                ! x,z center; y face
                E(x,y,z,2) = -cos(2._PRC*M_PI*(xVal + (yVal-0.5_PRC*dx(2)) + zVal)) / sqrt(2._PRC)
                ! x,y center; z face
                E(x,y,z,3) = cos(2._PRC*M_PI*(xVal + yVal + (zVal-0.5_PRC*dx(3)))) / sqrt(2._PRC)
                ! y,z face; x center
                B(x,y,z,1) = cos(2._PRC*M_PI*(xVal + (yVal-0.5_PRC*dx(2)) + (zVal-0.5_PRC*dx(3)))) * sqrt(2._PRC/3._PRC)
                ! x,z face; y center
                B(x,y,z,2) = -cos(2._PRC*M_PI*((xVal-0.5_PRC*dx(1)) + yVal + (zVal-0.5_PRC*dx(3)))) / sqrt(6._PRC)
                ! x,y face; z center
                B(x,y,z,3) = -cos(2._PRC*M_PI*((xVal-0.5_PRC*dx(1)) + (yVal-0.5_PRC*dx(2)) + zVal)) / sqrt(6._PRC)

            enddo
        enddo
    enddo

end subroutine setup_electromagnetic_vacuum_wave_3d_plane_polarized

