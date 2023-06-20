
! This file is part of the muphyII multiphysics plasma simulation project.
! 
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at https://mozilla.org/MPL/2.0/.

! Positive and flux conservative (PFC) method (Filbet, Sonnendrücker, Bertrand 2001)
! https://doi.org/10.1006/jcph.2001.6818

! Backsubstitution method for v x B (Schmitz, Grauer 2006)
! https://arxiv.org/abs/physics/0603208

! Note:
! The order step_x->step_y->step_z is important because boundary cells
! are only filled, if they are needed by the following steps. This is
! also the reason why step_z does not need the dimensionality_x argument.

subroutine step_x_vlasov_pfc(f,dimX,dimV,BD,dx,dv,v_b,dt,dimensionality_x)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3),dimV(3),BD(3), dimensionality_x
    real(kind=PRC), intent(in) :: dt,dx(3),v_b(3),dv(3)
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3)) :: f
    integer :: x,y,z,vx,vy,vz,shift,signum,floor_a,BDY,BDZ
    real(kind=PRC) :: a,dt_dx,alpha,epsilonp,epsilonm,ff,fh,fb,flux_minus_2,flux_minus_1,flux_current

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
                        floor_a = floor(a)
                        alpha = abs(a-aint(a))
                        if (a < 0) then
                            signum = -1
                        else
                            signum = 1
                        end if

                        flux_minus_1 = 0._PRC
                        flux_current = 0._PRC

                        ! calculate flux
                        !$acc loop seq
                        do x = 0,dimX(1)+1
                            flux_minus_2 = flux_minus_1
                            flux_minus_1 = flux_current
                            shift = x + floor_a

                            ff = f(shift+signum,y,z,vx,vy,vz)
                            fh = f(shift,y,z,vx,vy,vz)
                            fb = f(shift-signum,y,z,vx,vy,vz)

                            if (fb-fh > 0._PRC) then
                                epsilonp = min(1._PRC, 2._PRC*fh/(fb-fh))
                            else
                                epsilonp = 1._PRC
                            end if

                            if (ff-fh > 0._PRC) then
                                epsilonm = min(1._PRC, 2._PRC*fh/(ff-fh))
                            else
                                epsilonm = 1._PRC
                            end if

                            flux_current = -signum*alpha*(fh + epsilonp/6._PRC*(1._PRC-alpha)*(2._PRC-alpha)*(fb-fh) &
                                    + epsilonm/6._PRC*(1._PRC-alpha)*(1._PRC+alpha)*(fh-ff))

                            if (x > 1) then
                                f(x-2,y,z,vx,vy,vz) = f(x-2,y,z,vx,vy,vz) + flux_minus_2 - flux_minus_1
                            end if
                        enddo
                        f(dimX(1),y,z,vx,vy,vz) = f(dimX(1),y,z,vx,vy,vz) + flux_minus_1 - flux_current
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel
    !$acc end data

end subroutine step_x_vlasov_pfc


subroutine step_y_vlasov_pfc(f,dimX,dimV,BD,dx,dv,v_b,dt,dimensionality_x)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3),dimV(3),BD(3),dimensionality_x
    real(kind=PRC), intent(in) :: dt,dx(3),v_b(3),dv(3)
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3)) :: f

    integer :: x,y,z,vx,vy,vz,shift,signum,floor_a,BDZ
    real(kind=PRC) :: a,dt_dx,alpha,epsilonp,epsilonm,ff,fh,fb,flux_minus_2,flux_minus_1,flux_current

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
                        floor_a = floor(a)
                        alpha = abs(a-aint(a))
                        if (a < 0) then
                            signum = -1
                        else
                            signum = 1
                        end if

                        flux_minus_1 = 0._PRC
                        flux_current = 0._PRC

                        ! calculate flux
                        !$acc loop seq
                        do y = 0,dimX(2)+1
                            flux_minus_2 = flux_minus_1
                            flux_minus_1 = flux_current
                            shift = y + floor_a

                            ff = f(x,shift+signum,z,vx,vy,vz)
                            fh = f(x,shift,z,vx,vy,vz)
                            fb = f(x,shift-signum,z,vx,vy,vz)

                            if (fb-fh > 0._PRC) then
                                epsilonp = min(1._PRC, 2._PRC*fh/(fb-fh))
                            else
                                epsilonp = 1._PRC
                            end if

                            if (ff-fh > 0._PRC) then
                                epsilonm = min(1._PRC, 2._PRC*fh/(ff-fh))
                            else
                                epsilonm = 1._PRC
                            end if

                            flux_current = -signum*alpha*(fh + epsilonp/6._PRC*(1._PRC-alpha)*(2._PRC-alpha)*(fb-fh) &
                                    + epsilonm/6._PRC*(1._PRC-alpha)*(1._PRC+alpha)*(fh-ff))

                            if (y > 1) then
                                f(x,y-2,z,vx,vy,vz) = f(x,y-2,z,vx,vy,vz) + flux_minus_2 - flux_minus_1
                            end if
                        enddo
                        f(x,dimX(2),z,vx,vy,vz) = f(x,dimX(2),z,vx,vy,vz) + flux_minus_1 - flux_current
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel
    !$acc end data

end subroutine step_y_vlasov_pfc


subroutine step_z_vlasov_pfc(f,dimX,dimV,BD,dx,dv,v_b,dt)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3),dimV(3),BD(3)
    real(kind=PRC), intent(in) :: dt,dx(3),v_b(3),dv(3)
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3)) :: f

    integer :: x,y,z,vx,vy,vz,shift,signum,floor_a
    real(kind=PRC) :: a,dt_dx,alpha,epsilonp,epsilonm,ff,fh,fb,flux_minus_2,flux_minus_1,flux_current

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
                        floor_a = floor(a)
                        alpha = abs(a-aint(a))
                        if (a < 0) then
                            signum = -1
                        else
                            signum = 1
                        end if

                        flux_minus_1 = 0._PRC
                        flux_current = 0._PRC

                        ! calculate flux
                        !$acc loop seq
                        do z = 0,dimX(3)+1
                            flux_minus_2 = flux_minus_1
                            flux_minus_1 = flux_current
                            shift = z + floor_a

                            ff = f(x,y,shift+signum,vx,vy,vz)
                            fh = f(x,y,shift,vx,vy,vz)
                            fb = f(x,y,shift-signum,vx,vy,vz)

                            if (fb-fh > 0._PRC) then
                                epsilonp = min(1._PRC, 2._PRC*fh/(fb-fh))
                            else
                                epsilonp = 1._PRC
                            end if

                            if (ff-fh > 0._PRC) then
                                epsilonm = min(1._PRC, 2._PRC*fh/(ff-fh))
                            else
                                epsilonm = 1._PRC
                            end if

                            flux_current = -signum*alpha*(fh + epsilonp/6._PRC*(1._PRC-alpha)*(2._PRC-alpha)*(fb-fh) &
                                    + epsilonm/6._PRC*(1._PRC-alpha)*(1._PRC+alpha)*(fh-ff))

                            if (z > 1) then
                                f(x,y,z-2,vx,vy,vz) = f(x,y,z-2,vx,vy,vz) + flux_minus_2 - flux_minus_1
                            end if
                        enddo
                        f(x,y,dimX(3),vx,vy,vz) = f(x,y,dimX(3),vx,vy,vz) + flux_minus_1 - flux_current
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel
    !$acc end data

end subroutine step_z_vlasov_pfc


subroutine step_vx_vlasov_pfc(f,E,B,dimX,dimV,BD,dv,v_b,q,m,dt,dimensionality_x,include_boundary_cells)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), dimensionality_x
    logical, intent(in) :: include_boundary_cells
    real(kind=PRC), intent(in) :: dt, dv(3), v_b(3), q, m
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: E, B
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3)) :: f

    integer :: x, y, z, vx, vy, vz, signum, shift, BDX, BDY, BDZ
    real(kind=PRC) :: a, qmd, vxo, vyo, vzo, vxn, tx, ty, tz, sy, sz, temp, &
            alpha, epsilonp, epsilonm, ff, fh, fb, flux_minus_2, flux_minus_1, flux_current

    qmd = 0.5_PRC*q/m*dt

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

                        tx = -qmd*B(x,y,z,1)
                        ty = -qmd*B(x,y,z,2)
                        tz = -qmd*B(x,y,z,3)
                        temp = 2._PRC / (1._PRC + tx*tx + ty*ty + tz*tz)
                        sy = temp*ty
                        sz = temp*tz

                        flux_minus_1 = 0._PRC
                        flux_current = 0._PRC

                        ! calculate flux
                        ! TODO check boundary of velocity space
                        !$acc loop seq
                        do vx = 1+2,dimV(1)-1
                            flux_minus_2 = flux_minus_1
                            flux_minus_1 = flux_current

                            vxo = v_b(1)+(vx-1)*dv(1) - qmd*E(x,y,z,1)
                            vyo = v_b(2)+(vy-0.5_PRC)*dv(2) + qmd*E(x,y,z,2)
                            vzo = v_b(3)+(vz-0.5_PRC)*dv(3) + qmd*E(x,y,z,3)
                            vxn = (vxo - (tx*sy-sz)*vyo - (tx*sz+sy)*vzo)/(1._PRC-ty*sy-tz*sz) - qmd*E(x,y,z,1)
                            a = -(v_b(1)+(vx-1)*dv(1) - vxn) / dv(1)

                            alpha = abs(a-aint(a))
                            if (a < 0) then
                                signum = -1
                            else
                                signum = 1
                            end if
                            shift = vx + floor(a)

                            ff = f(x,y,z,shift+signum,vy,vz)
                            fh = f(x,y,z,shift,vy,vz)
                            fb = f(x,y,z,shift-signum,vy,vz)

                            if (fb-fh > 0._PRC) then
                                epsilonp = min(1._PRC, 2._PRC*fh/(fb-fh))
                            else
                                epsilonp = 1._PRC
                            endif

                            if (ff-fh > 0._PRC) then
                                epsilonm = min(1._PRC, 2._PRC*fh/(ff-fh))
                            else
                                epsilonm = 1._PRC
                            endif

                            flux_current = -signum*alpha*(fh + epsilonp/6._PRC*(1._PRC-alpha)*(2._PRC-alpha)*(fb-fh) &
                                    + epsilonm/6._PRC*(1._PRC-alpha)*(1._PRC+alpha)*(fh-ff))

                            if (vx > 1 + 3) then
                                f(x,y,z,vx-2,vy,vz) = f(x,y,z,vx-2,vy,vz) + flux_minus_2 - flux_minus_1
                            else if (vx == 1 + 3) then
                                f(x,y,z,vx-2,vy,vz) = f(x,y,z,vx-2,vy,vz) - flux_minus_1
                            end if

                        enddo
                        f(x,y,z,dimV(1)-2,vy,vz) = f(x,y,z,dimV(1)-2,vy,vz) + flux_minus_1 - flux_current
                        f(x,y,z,dimV(1)-1,vy,vz) = f(x,y,z,dimV(1)-1,vy,vz) + flux_current
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel
    !$acc end data

end subroutine step_vx_vlasov_pfc


subroutine step_vy_vlasov_pfc(f,E,B,dimX,dimV,BD,dv,v_b,q,m,dt,dimensionality_x,include_boundary_cells)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), dimensionality_x
    logical, intent(in) :: include_boundary_cells
    real(kind=PRC), intent(in) :: dt, dv(3), v_b(3), q, m
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: E, B
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3)) :: f

    integer :: x, y, z, vx, vy, vz, signum, shift, BDX, BDY, BDZ
    real(kind=PRC) :: a, qmd, vxo, vyo, vzo, vyn, tx, ty, tz, sx, sy, sz, temp, &
            alpha, epsilonp, epsilonm, ff, fh, fb, flux_minus_2, flux_minus_1, flux_current

    qmd = 0.5_PRC*q/m*dt

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

                        tx = -qmd*B(x,y,z,1)
                        ty = -qmd*B(x,y,z,2)
                        tz = -qmd*B(x,y,z,3)
                        temp = 2._PRC / (1._PRC + tx*tx + ty*ty + tz*tz)
                        sx = temp*tx
                        sy = temp*ty
                        sz = temp*tz

                        flux_minus_1 = 0._PRC
                        flux_current = 0._PRC

                        ! calculate flux
                        ! TODO check boundary of velocity space
                        !$acc loop seq
                        do vy = 1+2,dimV(2)-1
                            flux_minus_2 = flux_minus_1
                            flux_minus_1 = flux_current

                            vxo = v_b(1)+(vx-0.5_PRC)*dv(1) - qmd*E(x,y,z,1)
                            vyo = v_b(2)+(vy-1)*dv(2) - qmd*E(x,y,z,2)
                            vzo = v_b(3)+(vz-0.5_PRC)*dv(3) + qmd*E(x,y,z,3)

                            vyn = vxo*(ty*sx-sz) + vyo*(1._PRC-tx*sx-tz*sz) + &
                                    ((ty*sz+sx)/(1._PRC-ty*sy-tx*sx))*(vzo - vyo*(tz*sy-sx) - &
                                    vxo*(tz*sx+sy)) - qmd*E(x,y,z,2)

                            a = -(v_b(2)+(vy-1)*dv(2) - vyn) / dv(2)

                            alpha = abs(a-aint(a))
                            if (a < 0) then
                                signum = -1
                            else
                                signum = 1
                            end if
                            shift = vy + floor(a)

                            ff = f(x,y,z,vx,shift+signum,vz)
                            fh = f(x,y,z,vx,shift,vz)
                            fb = f(x,y,z,vx,shift-signum,vz)

                            if (fb-fh > 0._PRC) then
                                epsilonp = min(1._PRC, 2._PRC*fh/(fb-fh))
                            else
                                epsilonp = 1._PRC
                            endif

                            if (ff-fh > 0._PRC) then
                                epsilonm = min(1._PRC, 2._PRC*fh/(ff-fh))
                            else
                                epsilonm = 1._PRC
                            endif

                            flux_current = -signum*alpha*(fh + epsilonp/6._PRC*(1._PRC-alpha)*(2._PRC-alpha)*(fb-fh) &
                                 + epsilonm/6._PRC*(1._PRC-alpha)*(1._PRC+alpha)*(fh-ff))

                            if (vy > 1 + 3) then
                                f(x,y,z,vx,vy-2,vz) = f(x,y,z,vx,vy-2,vz) + flux_minus_2 - flux_minus_1
                            else if (vy == 1 + 3) then
                                f(x,y,z,vx,vy-2,vz) = f(x,y,z,vx,vy-2,vz) - flux_minus_1
                            end if
                        enddo
                        f(x,y,z,vx,dimV(2)-2,vz) = f(x,y,z,vx,dimV(2)-2,vz) + flux_minus_1 - flux_current
                        f(x,y,z,vx,dimV(2)-1,vz) = f(x,y,z,vx,dimV(2)-1,vz) + flux_current
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel
    !$acc end data

end subroutine step_vy_vlasov_pfc


subroutine step_vz_vlasov_pfc(f,E,B,dimX,dimV,BD,dv,v_b,q,m,dt,dimensionality_x,include_boundary_cells)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), dimensionality_x
    logical, intent(in) :: include_boundary_cells
    real(kind=PRC), intent(in) :: dt, dv(3), v_b(3), q, m
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: E, B
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3)) :: f

    integer :: x, y, z, vx, vy, vz, signum, shift, BDX, BDY, BDZ
    real(kind=PRC) :: a, qmd, vxo, vyo, vzo, vzn, tx, ty, tz, sx, sy, temp, &
            alpha, epsilonp, epsilonm, ff, fh, fb, flux_minus_2, flux_minus_1, flux_current

    qmd = 0.5_PRC*q/m*dt

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

                        tx = -qmd*B(x,y,z,1)
                        ty = -qmd*B(x,y,z,2)
                        tz = -qmd*B(x,y,z,3)
                        temp = 2._PRC / (1._PRC + tx*tx + ty*ty + tz*tz)
                        sx = temp*tx
                        sy = temp*ty

                        flux_minus_1 = 0._PRC
                        flux_current = 0._PRC

                        ! calculate flux
                        ! TODO check boundary of velocity space
                        !$acc loop seq
                        do vz = 1+2,dimV(3)-1
                            flux_minus_2 = flux_minus_1
                            flux_minus_1 = flux_current

                            vxo = v_b(1)+(vx-0.5_PRC)*dv(1) - qmd*E(x,y,z,1)
                            vyo = v_b(2)+(vy-0.5_PRC)*dv(2) - qmd*E(x,y,z,2)
                            vzo = v_b(3)+(vz-1)*dv(3) - qmd*E(x,y,z,3)

                            vzn = vzo*(1._PRC-ty*sy-tx*sx) + vyo*(tz*sy-sx) + vxo*(tz*sx+sy) - qmd*E(x,y,z,3)

                            a = -(v_b(3)+(vz-1)*dv(3) - vzn) / dv(3)

                            alpha = abs(a-aint(a))
                            if (a < 0) then
                                signum = -1
                            else
                                signum = 1
                            end if
                            shift = vz + floor(a)

                            ff = f(x,y,z,vx,vy,shift+signum)
                            fh = f(x,y,z,vx,vy,shift)
                            fb = f(x,y,z,vx,vy,shift-signum)

                            if (fb-fh > 0._PRC) then
                                epsilonp = min(1._PRC, 2._PRC*fh/(fb-fh))
                            else
                                epsilonp = 1._PRC
                            endif

                            if (ff-fh > 0._PRC) then
                                epsilonm = min(1._PRC, 2._PRC*fh/(ff-fh))
                            else
                                epsilonm = 1._PRC
                            endif

                            flux_current = -signum*alpha*(fh + epsilonp/6._PRC*(1._PRC-alpha)*(2._PRC-alpha)*(fb-fh) &
                                    + epsilonm/6._PRC*(1._PRC-alpha)*(1._PRC+alpha)*(fh-ff))

                            if (vz>1+3) then
                                f(x,y,z,vx,vy,vz-2) = f(x,y,z,vx,vy,vz-2) + flux_minus_2 - flux_minus_1
                            else if (vz == 1 + 3) then
                                f(x,y,z,vx,vy,vz-2) = f(x,y,z,vx,vy,vz-2) - flux_minus_1
                            end if
                        enddo
                        f(x,y,z,vx,vy,dimV(3)-2) = f(x,y,z,vx,vy,dimV(3)-2) + flux_minus_1 - flux_current
                        f(x,y,z,vx,vy,dimV(3)-1) = f(x,y,z,vx,vy,dimV(3)-1) + flux_current
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$acc end parallel
    !$acc end data

end subroutine step_vz_vlasov_pfc
