
! This file is part of the muphyII multiphysics plasma simulation project.
! 
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at https://mozilla.org/MPL/2.0/.

! _______ Symmetric-like boundary conditions for values on cell faces _______
! For values on the cell faces the value from cell [dimX+1] is
! excactly at the border so that it has no counterpart in case of the
! symmetric/mirrored and antisymmetric/antimirrored boundary conditions.
! Thus it is handled the following way:
! symmetric: constant extrapolation from the neighbour cell
! mirrored: linear extrapolation from the two cells to the left
! antisymmetric/antimirrored: 0

! _______ Zero boundary condition for reflecting/conducting walls _______
! Fields that are zero at a conducting wall are on the cell face.
! Since the face of the first actual cell is excactly at the
! boundary, for the zero boundary condition not only the boundary
! cells are set to zero but also the first actual cell.
! Reflected particles lead to zero average velocity in the first actual
! cell, which is therefore also set to zero when the zero boundary
! condition is used.


! TODO diagonal neighbours for boundary_schwarz and possibly for vlasov

subroutine boundary_dummy_dimension_handling(u,direction,dimX,BD,ncomponents)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), direction, ncomponents
    real(kind=PRC), intent(inout) :: u(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),ncomponents)

    ! For 1D and 2D simulations the unused dimensions can still have
    ! boundary cells which always carry the same value as the only
    ! actual cell in that dimension (periodic boundary conditions).
    ! This function copies the values of the actual cell to the boundary cells.

    integer :: i, x, y, z, c
    
    !$acc data present(u)
    select case (direction)
        case (0) ! x lower
            !$acc parallel
            !$acc loop gang collapse(2)
            do c = 1,ncomponents
                do z = -BD(3),dimX(3)+BD(3)
                    !$acc loop vector collapse(2)
                    do y = -BD(2),dimX(2)+BD(2)
                        do i=1,BD(1)
                            u(-i,y,z,c) = u(0,y,z,c)
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case (1) ! x upper
            !$acc parallel
            !$acc loop gang collapse(2)
            do c = 1,ncomponents
                do z = -BD(3),dimX(3)+BD(3)
                    !$acc loop vector collapse(2)
                    do y = -BD(2),dimX(2)+BD(2)
                        do i=1,BD(1)
                            u(dimX(1)+i,y,z,c) = u(0,y,z,c)
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case (2) ! y lower
            !$acc parallel
            !$acc loop gang collapse(2)
            do c = 1,ncomponents
                do z = -BD(3),dimX(3)+BD(3)
                    !$acc loop vector collapse(2)
                    do i=1,BD(2)
                        do x = -BD(1),dimX(1)+BD(1)
                            u(x,-i,z,c) = u(x,0,z,c)
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case (3) ! y upper
            !$acc parallel
            !$acc loop gang collapse(2)
            do c = 1,ncomponents
                do z = -BD(3),dimX(3)+BD(3)
                    !$acc loop vector collapse(2)
                    do i=1,BD(2)
                        do x = -BD(1),dimX(1)+BD(1)
                            u(x,dimX(2)+i,z,c) = u(x,0,z,c)
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case (4) ! z lower
            !$acc parallel
            !$acc loop gang collapse(3)
            do c = 1,ncomponents
                do i=1,BD(3)
                    do y = -BD(2),dimX(2)+BD(2)
                        !$acc loop vector
                        do x = -BD(1),dimX(1)+BD(1)
                            u(x,y,-i,c) = u(x,y,0,c)
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case (5) ! z upper
            !$acc parallel
            !$acc loop gang collapse(3)
            do c = 1,ncomponents
                do i=1,BD(3)
                    do y = -BD(2),dimX(2)+BD(2)
                        !$acc loop vector
                        do x = -BD(1),dimX(1)+BD(1)
                            u(x,y,dimX(3)+i,c) = u(x,y,0,c)
                        end do
                    end do
                end do
            end do
            !$acc end parallel
    end select
    !$acc end data
            
end subroutine boundary_dummy_dimension_handling


subroutine boundary_extract_buffer_x(u,out_buffer,side,dimX,BD,ncomponents)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), side, ncomponents
    real(kind=PRC), intent(in) :: u(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),ncomponents)
    real(kind=PRC), intent(inout) :: out_buffer(BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),ncomponents)

    integer :: i, y, z, c

    !$acc data present(u) present_or_copyout(out_buffer)
    if (side == 0) then ! lower boundary
        !$acc parallel
        !$acc loop gang collapse(2)
        do c = 1,ncomponents
            do z = -BD(3),dimX(3)+BD(3)
                !$acc loop vector collapse(2)
                do y = -BD(2),dimX(2)+BD(2)
                    do i=1,BD(1)
                        out_buffer(i,y,z,c) = u(i-1,y,z,c)
                    end do
                end do
            end do
        end do
        !$acc end parallel
    else if (side == 1) then ! upper boundary
        !$acc parallel
        !$acc loop gang collapse(2)
        do c = 1,ncomponents
            do z = -BD(3),dimX(3)+BD(3)
                !$acc loop vector collapse(2)
                do y = -BD(2),dimX(2)+BD(2)
                    do i=1,BD(1)
                        out_buffer(i,y,z,c) = u(dimX(1)+1-i,y,z,c)
                    end do
                end do
            end do
        end do
        !$acc end parallel
    end if
    !$acc end data

end subroutine boundary_extract_buffer_x


subroutine boundary_extract_buffer_y(u,out_buffer,side,dimX,BD,ncomponents)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), side, ncomponents
    real(kind=PRC), intent(in) :: u(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),ncomponents)
    real(kind=PRC), intent(inout) :: out_buffer(-BD(1):dimX(1)+BD(1),BD(2),-BD(3):dimX(3)+BD(3),ncomponents)

    integer :: i, x, z, c

    !$acc data present(u) present_or_copyout(out_buffer)
    if (side == 0) then ! lower boundary
        !$acc parallel
        !$acc loop gang collapse(2)
        do c = 1,ncomponents
            do z = -BD(3),dimX(3)+BD(3)
                !$acc loop vector collapse(2)
                do i=1,BD(2)
                    do x = -BD(1),dimX(1)+BD(1)
                        out_buffer(x,i,z,c) = u(x,i-1,z,c)
                    end do
                end do
            end do
        end do
        !$acc end parallel
    else if (side == 1) then ! upper boundary
        !$acc parallel
        !$acc loop gang collapse(2)
        do c = 1,ncomponents
            do z = -BD(3),dimX(3)+BD(3)
                !$acc loop vector collapse(2)
                do i=1,BD(2)
                    do x = -BD(1),dimX(1)+BD(1)
                        out_buffer(x,i,z,c) = u(x,dimX(2)+1-i,z,c)
                    end do
                end do
            end do
        end do
        !$acc end parallel
    end if
    !$acc end data

end subroutine boundary_extract_buffer_y


subroutine boundary_extract_buffer_z(u,out_buffer,side,dimX,BD,ncomponents)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), side, ncomponents
    real(kind=PRC), intent(in) :: u(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),ncomponents)
    real(kind=PRC), intent(inout) :: out_buffer(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),BD(3),ncomponents)

    integer :: i, x, y, c

    !$acc data present(u) present_or_copyout(out_buffer)
    if (side == 0) then ! lower boundary
        !$acc parallel
        !$acc loop gang collapse(3)
        do c = 1,ncomponents
            do i=1,BD(3)
                do y = -BD(2),dimX(2)+BD(2)
                    !$acc loop vector
                    do x = -BD(1),dimX(1)+BD(1)
                        out_buffer(x,y,i,c) = u(x,y,i-1,c)
                    end do
                end do
            end do
        end do
        !$acc end parallel
    else if (side == 1) then ! upper boundary
        !$acc parallel
        !$acc loop gang collapse(3)
        do c = 1,ncomponents
            do i=1,BD(3)
                do y = -BD(2),dimX(2)+BD(2)
                    !$acc loop vector
                    do x = -BD(1),dimX(1)+BD(1)
                        out_buffer(x,y,i,c) = u(x,y,dimX(3)+1-i,c)
                    end do
                end do
            end do
        end do
        !$acc end parallel
    end if
    !$acc end data

end subroutine boundary_extract_buffer_z


subroutine boundary_apply_buffer_x(u,in_buffer,side,dimX,BD,ncomponents)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), side, ncomponents
    real(kind=PRC), intent(inout) :: u(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),ncomponents)
    real(kind=PRC), intent(in) :: in_buffer(BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),ncomponents)

    integer :: i, y, z, c

    !$acc data present(u) present_or_copyin(in_buffer)
    if (side == 0) then ! lower boundary
        !$acc parallel
        !$acc loop gang collapse(2)
        do c = 1,ncomponents
            do z = -BD(3),dimX(3)+BD(3)
                !$acc loop vector collapse(2)
                do y = -BD(2),dimX(2)+BD(2)
                    do i=1,BD(1)
                        u(-i,y,z,c) = in_buffer(i,y,z,c)
                    end do
                end do
            end do
        end do
        !$acc end parallel
    else if (side == 1) then ! upper boundary
        !$acc parallel
        !$acc loop gang collapse(2)
        do c = 1,ncomponents
            do z = -BD(3),dimX(3)+BD(3)
                !$acc loop vector collapse(2)
                do y = -BD(2),dimX(2)+BD(2)
                    do i=1,BD(1)
                        u(dimX(1)+i,y,z,c) = in_buffer(i,y,z,c)
                    end do
                end do
            end do
        end do
        !$acc end parallel
    end if
    !$acc end data

end subroutine boundary_apply_buffer_x


subroutine boundary_apply_buffer_y(u,in_buffer,side,dimX,BD,ncomponents)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), side, ncomponents
    real(kind=PRC), intent(inout) :: u(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),ncomponents)
    real(kind=PRC), intent(in) :: in_buffer(-BD(1):dimX(1)+BD(1),BD(2),-BD(3):dimX(3)+BD(3),ncomponents)

    integer :: i, x, z, c

    !$acc data present(u) present_or_copyin(in_buffer)
    if (side == 0) then ! lower boundary
        !$acc parallel
        !$acc loop gang collapse(2)
        do c = 1,ncomponents
            do z = -BD(3),dimX(3)+BD(3)
                !$acc loop vector collapse(2)
                do i=1,BD(2)
                    do x = -BD(1),dimX(1)+BD(1)
                        u(x,-i,z,c) = in_buffer(x,i,z,c)
                    end do
                end do
            end do
        end do
        !$acc end parallel
    else if (side == 1) then ! upper boundary
        !$acc parallel
        !$acc loop gang collapse(2)
        do c = 1,ncomponents
            do z = -BD(3),dimX(3)+BD(3)
                !$acc loop vector collapse(2)
                do i=1,BD(2)
                    do x = -BD(1),dimX(1)+BD(1)
                        u(x,dimX(2)+i,z,c) = in_buffer(x,i,z,c)
                    end do
                end do
            end do
        end do
        !$acc end parallel
    end if
    !$acc end data

end subroutine boundary_apply_buffer_y


subroutine boundary_apply_buffer_z(u,in_buffer,side,dimX,BD,ncomponents)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), side, ncomponents
    real(kind=PRC), intent(inout) :: u(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),ncomponents)
    real(kind=PRC), intent(in) :: in_buffer(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),BD(3),ncomponents)

    integer :: i, x, y, c

    !$acc data present(u) present_or_copyin(in_buffer)
    if (side == 0) then ! lower boundary
        !$acc parallel
        !$acc loop gang collapse(3)
        do c = 1,ncomponents
            do i=1,BD(3)
                do y = -BD(2),dimX(2)+BD(2)
                    !$acc loop vector
                    do x = -BD(1),dimX(1)+BD(1)
                        u(x,y,-i,c) = in_buffer(x,y,i,c)
                    end do
                end do
            end do
        end do
        !$acc end parallel
    else if (side == 1) then ! upper boundary
        !$acc parallel
        !$acc loop gang collapse(3)
        do c = 1,ncomponents
            do i=1,BD(3)
                do y = -BD(2),dimX(2)+BD(2)
                    !$acc loop vector
                    do x = -BD(1),dimX(1)+BD(1)
                        u(x,y,dimX(3)+i,c) = in_buffer(x,y,i,c)
                    end do
                end do
            end do
        end do
        !$acc end parallel
    end if
    !$acc end data

end subroutine boundary_apply_buffer_z


subroutine boundary_apply_buffer_schwarz_x(u,in_buffer,side,dimX,BD)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), side
    real(kind=PRC), intent(inout) :: u(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3))
    real(kind=PRC), intent(in) :: in_buffer(BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3))

    ! exchange most outer boundary cell
    ! used in additive Schwarz Poisson solver

    !$acc data present(u) copyin(in_buffer)

    if (side == 0) then ! lower boundary
        !$acc kernels
        u(-BD(1),:,:) = in_buffer(BD(1),:,:)
        !$acc end kernels
    else if (side == 1) then ! upper boundary
        !$acc kernels
        u(dimX(1)+BD(1),:,:) = in_buffer(BD(1),:,:)
        !$acc end kernels
    end if

    !$acc end data

end subroutine boundary_apply_buffer_schwarz_x


subroutine boundary_apply_buffer_schwarz_y(u,in_buffer,side,dimX,BD)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), side
    real(kind=PRC), intent(inout) :: u(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3))
    real(kind=PRC), intent(in) :: in_buffer(-BD(1):dimX(1)+BD(1),BD(2),-BD(3):dimX(3)+BD(3))

    ! exchange most outer boundary cell
    ! used in additive Schwarz Poisson solver

    !$acc data present(u) copyin(in_buffer)

    if (side == 0) then ! lower boundary
        !$acc kernels
        u(:,-BD(2),:) = in_buffer(:,BD(2),:)
        !$acc end kernels
    else if (side == 1) then ! upper boundary
        !$acc kernels
        u(:,dimX(2)+BD(2),:) = in_buffer(:,BD(2),:)
        !$acc end kernels
    end if

    !$acc end data

end subroutine boundary_apply_buffer_schwarz_y


subroutine boundary_apply_buffer_schwarz_z(u,in_buffer,side,dimX,BD)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), side
    real(kind=PRC), intent(inout) :: u(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3))
    real(kind=PRC), intent(in) :: in_buffer(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),BD(3))

    ! exchange most outer boundary cell
    ! used in additive Schwarz Poisson solver

    !$acc data present(u) copyin(in_buffer)

    if (side == 0) then ! lower boundary
        !$acc kernels
        u(:,:,-BD(3)) = in_buffer(:,:,BD(3))
        !$acc end kernels
    else if (side == 1) then ! upper boundary
        !$acc kernels
        u(:,:,dimX(3)+BD(3)) = in_buffer(:,:,BD(3))
        !$acc end kernels
    end if

    !$acc end data

end subroutine boundary_apply_buffer_schwarz_z


subroutine boundary_apply_boundary_conditions_x(u,side,dimX,BD,ncomponents,bdCond,orientation)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), side, ncomponents
    real(kind=PRC), intent(inout) :: u(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),ncomponents)
    character(len=6), intent(in) :: bdCond(ncomponents)
    character(len=3), intent(in) :: orientation(ncomponents)

    integer :: i, c

    !$acc data present(u)

    if (side == 0) then ! lower boundary
        do c = 1,ncomponents
            select case (bdCond(c)(1:1))
            case ('z') ! zero
                do i=0,BD(1)
                    !$acc kernels
                    u(-i,:,:,c) = 0._PRC
                    !$acc end kernels
                enddo
            case ('s','m') ! symmetric/mirrored
                if (orientation(c)(1:1) == 'c') then
                    do i=1,BD(1)
                        !$acc kernels
                        u(-i,:,:,c) = u(i-1,:,:,c)
                        !$acc end kernels
                    enddo
                else
                    do i=1,BD(1)
                        !$acc kernels
                        u(-i,:,:,c) = u(i,:,:,c)
                        !$acc end kernels
                    enddo
                endif
            case ('a') ! antisymmetric/antimirrored
                if (orientation(c)(1:1) == 'c') then
                    do i=1,BD(1)
                        !$acc kernels
                        u(-i,:,:,c) = -u(i-1,:,:,c)
                        !$acc end kernels
                    enddo
                else
                    do i=1,BD(1)
                        !$acc kernels
                        u(-i,:,:,c) = -u(i,:,:,c)
                        !$acc end kernels
                    enddo
                endif
            case ('q') ! antiperiodic
                do i=1,BD(1)
                    !$acc kernels
                    u(-i,:,:,c) = -u(-i,:,:,c)
                    !$acc end kernels
                end do
            end select
        end do

    else if (side == 1) then ! upper boundary
        do c = 1,ncomponents
            select case (bdCond(c)(2:2))
            case ('z') ! zero
                do i=1,BD(1)
                    !$acc kernels
                    u(dimX(1)+i,:,:,c) = 0._PRC
                    !$acc end kernels
                enddo
                if (orientation(c)(1:1) == 'c') then
                    !$acc kernels
                    u(dimX(1),:,:,c) = 0._PRC
                    !$acc end kernels
                endif
            case ('s') ! symmetric
                if (orientation(c)(1:1) == 'c') then
                    do i=1,BD(1)
                        !$acc kernels
                        u(dimX(1)+i,:,:,c) = u(dimX(1)+1-i,:,:,c)
                        !$acc end kernels
                    enddo
                else
                    do i=2,BD(1)
                        !$acc kernels
                        u(dimX(1)+i,:,:,c) = u(dimX(1)+2-i,:,:,c)
                        !$acc end kernels
                    enddo
                    !$acc kernels
                    u(dimX(1)+1,:,:,c) = u(dimX(1),:,:,c)
                    !$acc end kernels
                endif
            case ('m') ! mirrored
                if (orientation(c)(1:1) == 'c') then
                    do i=1,BD(1)
                        !$acc kernels
                        u(dimX(1)+i,:,:,c) = u(dimX(1)+1-i,:,:,c)
                        !$acc end kernels
                    enddo
                else
                    do i=2,BD(1)
                        !$acc kernels
                        u(dimX(1)+i,:,:,c) = u(dimX(1)+2-i,:,:,c)
                        !$acc end kernels
                    enddo
                    !$acc kernels
                    u(dimX(1)+1,:,:,c) = 2._PRC*u(dimX(1),:,:,c)-u(dimX(1)-1,:,:,c)
                    !$acc end kernels
                endif
            case ('a') ! antisymmetric/antimirrored
                if (orientation(c)(1:1) == 'c') then
                    do i=1,BD(1)
                        !$acc kernels
                        u(dimX(1)+i,:,:,c) = -u(dimX(1)+1-i,:,:,c)
                        !$acc end kernels
                    enddo
                else
                    do i=2,BD(1)
                        !$acc kernels
                        u(dimX(1)+i,:,:,c) = -u(dimX(1)+2-i,:,:,c)
                        !$acc end kernels
                    enddo
                    !$acc kernels
                    u(dimX(1)+1,:,:,c) = 0._PRC
                    !$acc end kernels
                endif
            case ('q') ! antiperiodic
                do i=1,BD(1)
                    !$acc kernels
                    u(dimX(1)+i,:,:,c) = -u(dimX(1)+i,:,:,c)
                    !$acc end kernels
                end do
            end select
        end do
    end if

    !$acc end data

end subroutine boundary_apply_boundary_conditions_x


subroutine boundary_apply_boundary_conditions_y(u,side,dimX,BD,ncomponents,bdCond,orientation)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), side, ncomponents
    real(kind=PRC), intent(inout) :: u(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),ncomponents)
    character(len=6), intent(in) :: bdCond(ncomponents)
    character(len=3), intent(in) :: orientation(ncomponents)

    integer :: i, c

    !$acc data present(u)

    if (side == 0) then ! lower boundary
        do c = 1,ncomponents
            select case (bdCond(c)(3:3))
            case ('z') ! zero
                do i=0,BD(2)
                    !$acc kernels
                    u(:,-i,:,c) = 0._PRC
                    !$acc end kernels
                enddo
            case ('s','m') ! symmetric/mirrored
                if (orientation(c)(2:2) == 'c') then
                    do i=1,BD(2)
                        !$acc kernels
                        u(:,-i,:,c) = u(:,i-1,:,c)
                        !$acc end kernels
                    enddo
                else
                    do i=1,BD(2)
                        !$acc kernels
                        u(:,-i,:,c) = u(:,i,:,c)
                        !$acc end kernels
                    enddo
                endif
            case ('a') ! antisymmetric/antimirrored
                if (orientation(c)(2:2) == 'c') then
                    do i=1,BD(2)
                        !$acc kernels
                        u(:,-i,:,c) = -u(:,i-1,:,c)
                        !$acc end kernels
                    enddo
                else
                    do i=1,BD(2)
                        !$acc kernels
                        u(:,-i,:,c) = -u(:,i,:,c)
                        !$acc end kernels
                    enddo
                endif
            case ('q') ! antiperiodic
                do i=1,BD(2)
                    !$acc kernels
                    u(:,-i,:,c) = -u(:,-i,:,c)
                    !$acc end kernels
                end do
            end select
        end do

    else if (side == 1) then ! upper boundary
        do c = 1,ncomponents
            select case (bdCond(c)(4:4))
            case ('z') ! zero
                do i=1,BD(2)
                    !$acc kernels
                    u(:,dimX(2)+i,:,c) = 0._PRC
                    !$acc end kernels
                enddo
                if (orientation(c)(2:2) == 'c') then
                    !$acc kernels
                    u(:,dimX(2),:,c) = 0._PRC
                    !$acc end kernels
                endif
            case ('s') ! symmetric
                if (orientation(c)(2:2) == 'c') then
                    do i=1,BD(2)
                        !$acc kernels
                        u(:,dimX(2)+i,:,c) = u(:,dimX(2)+1-i,:,c)
                        !$acc end kernels
                    enddo
                else
                    do i=2,BD(2)
                        !$acc kernels
                        u(:,dimX(2)+i,:,c) = u(:,dimX(2)+2-i,:,c)
                        !$acc end kernels
                    enddo
                    !$acc kernels
                    u(:,dimX(2)+1,:,c) = u(:,dimX(2),:,c)
                    !$acc end kernels
                endif
            case ('m') ! mirrored
                if (orientation(c)(2:2) == 'c') then
                    do i=1,BD(2)
                        !$acc kernels
                        u(:,dimX(2)+i,:,c) = u(:,dimX(2)+1-i,:,c)
                        !$acc end kernels
                    enddo
                else
                    do i=2,BD(2)
                        !$acc kernels
                        u(:,dimX(2)+i,:,c) = u(:,dimX(2)+2-i,:,c)
                        !$acc end kernels
                    enddo
                    !$acc kernels
                    u(:,dimX(2)+1,:,c) = 2._PRC*u(:,dimX(2),:,c)-u(:,dimX(2)-1,:,c)
                    !$acc end kernels
                endif
            case ('a') ! antisymmetric/antimirrored
                if (orientation(c)(2:2) == 'c') then
                    do i=1,BD(2)
                        !$acc kernels
                        u(:,dimX(2)+i,:,c) = -u(:,dimX(2)+1-i,:,c)
                        !$acc end kernels
                    enddo
                else
                    do i=2,BD(2)
                        !$acc kernels
                        u(:,dimX(2)+i,:,c) = -u(:,dimX(2)+2-i,:,c)
                        !$acc end kernels
                    enddo
                    !$acc kernels
                    u(:,dimX(2)+1,:,c) = 0._PRC
                    !$acc end kernels
                endif
            case ('q') ! antiperiodic
                do i=1,BD(2)
                    !$acc kernels
                    u(:,dimX(2)+i,:,c) = -u(:,dimX(2)+i,:,c)
                    !$acc end kernels
                end do
            end select
        end do
    end if

    !$acc end data

end subroutine boundary_apply_boundary_conditions_y


subroutine boundary_apply_boundary_conditions_z(u,side,dimX,BD,ncomponents,bdCond,orientation)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), side, ncomponents
    real(kind=PRC), intent(inout) :: u(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),ncomponents)
    character(len=6), intent(in) :: bdCond(ncomponents)
    character(len=3), intent(in) :: orientation(ncomponents)

    integer :: i, c

    !$acc data present(u)

    if (side == 0) then ! lower boundary
        do c = 1,ncomponents
            select case (bdCond(c)(5:5))
            case ('z') ! zero
                do i=0,BD(3)
                    !$acc kernels
                    u(:,:,-i,c) = 0._PRC
                    !$acc end kernels
                enddo
            case ('s','m') ! symmetric/mirrored
                if (orientation(c)(3:3) == 'c') then
                    do i=1,BD(3)
                        !$acc kernels
                        u(:,:,-i,c) = u(:,:,i-1,c)
                        !$acc end kernels
                    enddo
                else
                    do i=1,BD(3)
                        !$acc kernels
                        u(:,:,-i,c) = u(:,:,i,c)
                        !$acc end kernels
                    enddo
                endif
            case ('a') ! antisymmetric/antimirrored
                if (orientation(c)(3:3) == 'c') then
                    do i=1,BD(3)
                        !$acc kernels
                        u(:,:,-i,c) = -u(:,:,i-1,c)
                        !$acc end kernels
                    enddo
                else
                    do i=1,BD(3)
                        !$acc kernels
                        u(:,:,-i,c) = -u(:,:,i,c)
                        !$acc end kernels
                    enddo
                endif
            case ('q') ! antiperiodic
                do i=1,BD(3)
                    !$acc kernels
                    u(:,:,-i,c) = -u(:,:,-i,c)
                    !$acc end kernels
                end do
            end select
        end do

    else if (side == 1) then ! upper boundary
        do c = 1,ncomponents
            select case (bdCond(c)(6:6))
            case ('z') ! zero
                do i=1,BD(3)
                    !$acc kernels
                    u(:,:,dimX(3)+i,c) = 0._PRC
                    !$acc end kernels
                enddo
                if (orientation(c)(3:3) == 'c') then
                    !$acc kernels
                    u(:,:,dimX(3),c) = 0._PRC
                    !$acc end kernels
                endif
            case ('s') ! symmetric
                if (orientation(c)(3:3) == 'c') then
                    do i=1,BD(3)
                        !$acc kernels
                        u(:,:,dimX(3)+i,c) = u(:,:,dimX(3)+1-i,c)
                        !$acc end kernels
                    enddo
                else
                    do i=2,BD(3)
                        !$acc kernels
                        u(:,:,dimX(3)+i,c) = u(:,:,dimX(3)+2-i,c)
                        !$acc end kernels
                    enddo
                    !$acc kernels
                    u(:,:,dimX(3)+1,c) = u(:,:,dimX(3),c)
                    !$acc end kernels
                endif
            case ('m') ! mirrored
                if (orientation(c)(3:3) == 'c') then
                    do i=1,BD(3)
                        !$acc kernels
                        u(:,:,dimX(3)+i,c) = u(:,:,dimX(3)+1-i,c)
                        !$acc end kernels
                    enddo
                else
                    do i=2,BD(3)
                        !$acc kernels
                        u(:,:,dimX(3)+i,c) = u(:,:,dimX(3)+2-i,c)
                        !$acc end kernels
                    enddo
                    !$acc kernels
                    u(:,:,dimX(3)+1,c) = 2._PRC*u(:,:,dimX(3),c)-u(:,:,dimX(3)-1,c)
                    !$acc end kernels
                endif
            case ('a') ! antisymmetric/antimirrored
                if (orientation(c)(3:3) == 'c') then
                    do i=1,BD(3)
                        !$acc kernels
                        u(:,:,dimX(3)+i,c) = -u(:,:,dimX(3)+1-i,c)
                        !$acc end kernels
                    enddo
                else
                    do i=2,BD(3)
                        !$acc kernels
                        u(:,:,dimX(3)+i,c) = -u(:,:,dimX(3)+2-i,c)
                        !$acc end kernels
                    enddo
                    !$acc kernels
                    u(:,:,dimX(3)+1,c) = 0._PRC
                    !$acc end kernels
                endif
            case ('q') ! antiperiodic
                do i=1,BD(3)
                    !$acc kernels
                    u(:,:,dimX(3)+i,c) = -u(:,:,dimX(3)+i,c)
                    !$acc end kernels
                end do
            end select
        end do
    end if

    !$acc end data

end subroutine boundary_apply_boundary_conditions_z


subroutine boundary_dummy_dimension_handling_f(f,direction,dimX,dimV,BD)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), direction
    real(kind=PRC), intent(inout) :: f(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))

    ! For 1D and 2D simulations the unused dimensions can still have
    ! boundary cells which always carry the same value as the only
    ! actual cell in that dimension (periodic boundary conditions).
    ! This function copies the values of the actual cell to the boundary cells.

    integer :: i,x,y,z,vx,vy,vz

    !$acc data present(f)
    select case (direction)
    case (0) ! x lower
        do i=1,BD(1)
            !$acc parallel
            !$acc loop gang collapse(3)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    do vx = 1,dimV(1)
                        !$acc loop vector collapse(2)
                        do z = -BD(3),dimX(3)+BD(3)
                            do y = -BD(2),dimX(2)+BD(2)
                                f(-i,y,z,vx,vy,vz) = f(0,y,z,vx,vy,vz)
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        end do
    case (1) ! x upper
        do i=1,BD(1)
            !$acc parallel
            !$acc loop gang collapse(3)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    do vx = 1,dimV(1)
                        !$acc loop vector collapse(2)
                        do z = -BD(3),dimX(3)+BD(3)
                            do y = -BD(2),dimX(2)+BD(2)
                                f(dimX(1)+i,y,z,vx,vy,vz) = f(0,y,z,vx,vy,vz)
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        end do
    case (2) ! y lower
        do i=1,BD(2)
            !$acc parallel
            !$acc loop gang collapse(3)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    do vx = 1,dimV(1)
                        !$acc loop vector collapse(2)
                        do z = -BD(3),dimX(3)+BD(3)
                            do x = -BD(1),dimX(1)+BD(1)
                                f(x,-i,z,vx,vy,vz) = f(x,0,z,vx,vy,vz)
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        end do
    case (3) ! y upper
        do i=1,BD(2)
            !$acc parallel
            !$acc loop gang collapse(3)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    do vx = 1,dimV(1)
                        !$acc loop vector collapse(2)
                        do z = -BD(3),dimX(3)+BD(3)
                            do x = -BD(1),dimX(1)+BD(1)
                                f(x,dimX(2)+i,z,vx,vy,vz) = f(x,0,z,vx,vy,vz)
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        end do
    case (4) ! z lower
        do i=1,BD(3)
            !$acc parallel
            !$acc loop gang collapse(3)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    do vx = 1,dimV(1)
                        !$acc loop vector collapse(2)
                        do y = -BD(2),dimX(2)+BD(2)
                            do x = -BD(1),dimX(1)+BD(1)
                                f(x,y,-i,vx,vy,vz) = f(x,y,0,vx,vy,vz)
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        end do
    case (5) ! z upper
        do i=1,BD(3)
            !$acc parallel
            !$acc loop gang collapse(3)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    do vx = 1,dimV(1)
                        !$acc loop vector collapse(2)
                        do y = -BD(2),dimX(2)+BD(2)
                            do x = -BD(1),dimX(1)+BD(1)
                                f(x,y,dimX(3)+i,vx,vy,vz) = f(x,y,0,vx,vy,vz)
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        end do
    end select
    !$acc end data
 
end subroutine boundary_dummy_dimension_handling_f


subroutine boundary_extract_buffer_f_x(f,out_buffer,side,dimX,dimV,BD)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), side
    real(kind=PRC), intent(in) :: f(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))
    real(kind=PRC), intent(inout) :: out_buffer(BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))

    integer :: i,y,z,vx,vy,vz

    !$acc data present(f) present_or_copyout(out_buffer)
    if (side == 0) then ! lower boundary
        !$acc parallel
        !$acc loop gang collapse(3)
        do vz = 1,dimV(3)
            do vy = 1,dimV(2)
                do vx = 1,dimV(1)
                    !$acc loop vector collapse(3)
                    do z = -BD(3),dimX(3)+BD(3)
                        do y = -BD(2),dimX(2)+BD(2)
                            do i=1,BD(1)
                                out_buffer(i,y,z,vx,vy,vz) = f(i-1,y,z,vx,vy,vz)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$acc end parallel
    else if (side == 1) then ! upper boundary
        !$acc parallel
        !$acc loop gang collapse(3)
        do vz = 1,dimV(3)
            do vy = 1,dimV(2)
                do vx = 1,dimV(1)
                    !$acc loop vector collapse(3)
                    do z = -BD(3),dimX(3)+BD(3)
                        do y = -BD(2),dimX(2)+BD(2)
                            do i=1,BD(1)
                                out_buffer(i,y,z,vx,vy,vz) = f(dimX(1)+1-i,y,z,vx,vy,vz)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$acc end parallel
    end if
    !$acc end data
    
end subroutine boundary_extract_buffer_f_x


subroutine boundary_extract_buffer_f_y(f,out_buffer,side,dimX,dimV,BD)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), side
    real(kind=PRC), intent(in) :: f(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))
    real(kind=PRC), intent(inout) :: out_buffer(-BD(1):dimX(1)+BD(1),BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))

    integer :: i,x,z,vx,vy,vz

    !$acc data present(f) present_or_copyout(out_buffer)
    if (side == 0) then ! lower boundary
        !$acc parallel
        !$acc loop gang collapse(3)
        do vz = 1,dimV(3)
            do vy = 1,dimV(2)
                do vx = 1,dimV(1)
                    !$acc loop vector collapse(3)
                    do z = -BD(3),dimX(3)+BD(3)
                        do i=1,BD(2)
                            do x = -BD(1),dimX(1)+BD(1)
                                out_buffer(x,i,z,vx,vy,vz) = f(x,i-1,z,vx,vy,vz)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$acc end parallel
    else if (side == 1) then ! upper boundary
        !$acc parallel
        !$acc loop gang collapse(3)
        do vz = 1,dimV(3)
            do vy = 1,dimV(2)
                do vx = 1,dimV(1)
                    !$acc loop vector collapse(3)
                    do z = -BD(3),dimX(3)+BD(3)
                        do i=1,BD(2)
                            do x = -BD(1),dimX(1)+BD(1)
                                out_buffer(x,i,z,vx,vy,vz) = f(x,dimX(2)+1-i,z,vx,vy,vz)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$acc end parallel
    end if
    !$acc end data

end subroutine boundary_extract_buffer_f_y


subroutine boundary_extract_buffer_f_z(f,out_buffer,side,dimX,dimV,BD)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), side
    real(kind=PRC), intent(in) :: f(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))
    real(kind=PRC), intent(inout) :: out_buffer(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),BD(3),dimV(1),dimV(2),dimV(3))

    integer :: i,x,y,vx,vy,vz

    !$acc data present(f) present_or_copyout(out_buffer)
    if (side == 0) then ! lower boundary
        !$acc parallel
        !$acc loop gang collapse(3)
        do vz = 1,dimV(3)
            do vy = 1,dimV(2)
                do vx = 1,dimV(1)
                    !$acc loop vector collapse(3)
                    do i=1,BD(3)
                        do y = -BD(2),dimX(2)+BD(2)
                            do x = -BD(1),dimX(1)+BD(1)
                                out_buffer(x,y,i,vx,vy,vz) = f(x,y,i-1,vx,vy,vz)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$acc end parallel
    else if (side == 1) then ! upper boundary
        !$acc parallel
        !$acc loop gang collapse(3)
        do vz = 1,dimV(3)
            do vy = 1,dimV(2)
                do vx = 1,dimV(1)
                    !$acc loop vector collapse(3)
                    do i=1,BD(3)
                        do y = -BD(2),dimX(2)+BD(2)
                            do x = -BD(1),dimX(1)+BD(1)
                                out_buffer(x,y,i,vx,vy,vz) = f(x,y,dimX(3)+1-i,vx,vy,vz)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$acc end parallel
    end if
    !$acc end data

end subroutine boundary_extract_buffer_f_z


subroutine boundary_apply_buffer_f_x(f,in_buffer,side,dimX,dimV,BD)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), side
    real(kind=PRC), intent(inout) :: f(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))
    real(kind=PRC), intent(in) :: in_buffer(BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))

    integer :: i,y,z,vx,vy,vz

    !$acc data present(f) present_or_copyin(in_buffer)
    if (side == 0) then ! lower boundary
        !$acc parallel
        !$acc loop gang collapse(3)
        do vz = 1,dimV(3)
            do vy = 1,dimV(2)
                do vx = 1,dimV(1)
                    !$acc loop vector collapse(3)
                    do z = -BD(3),dimX(3)+BD(3)
                        do y = -BD(2),dimX(2)+BD(2)
                            do i=1,BD(1)
                                f(-i,y,z,vx,vy,vz) = in_buffer(i,y,z,vx,vy,vz)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$acc end parallel
    else if (side == 1) then ! upper boundary
        !$acc parallel
        !$acc loop gang collapse(3)
        do vz = 1,dimV(3)
            do vy = 1,dimV(2)
                do vx = 1,dimV(1)
                    !$acc loop vector collapse(3)
                    do z = -BD(3),dimX(3)+BD(3)
                        do y = -BD(2),dimX(2)+BD(2)
                            do i=1,BD(1)
                                f(dimX(1)+i,y,z,vx,vy,vz) = in_buffer(i,y,z,vx,vy,vz)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$acc end parallel
    end if
    !$acc end data

end subroutine boundary_apply_buffer_f_x


subroutine boundary_apply_buffer_f_y(f,in_buffer,side,dimX,dimV,BD)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), side
    real(kind=PRC), intent(inout) :: f(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))
    real(kind=PRC), intent(in) :: in_buffer(-BD(1):dimX(1)+BD(1),BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))

    integer :: i,x,z,vx,vy,vz

    !$acc data present(f) present_or_copyin(in_buffer)
    if (side == 0) then ! lower boundary
        !$acc parallel
        !$acc loop gang collapse(3)
        do vz = 1,dimV(3)
            do vy = 1,dimV(2)
                do vx = 1,dimV(1)
                    !$acc loop vector collapse(3)
                    do z = -BD(3),dimX(3)+BD(3)
                        do i=1,BD(2)
                            do x = -BD(1),dimX(1)+BD(1)
                                f(x,-i,z,vx,vy,vz) = in_buffer(x,i,z,vx,vy,vz)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$acc end parallel
    else if (side == 1) then ! upper boundary
        !$acc parallel
        !$acc loop gang collapse(3)
        do vz = 1,dimV(3)
            do vy = 1,dimV(2)
                do vx = 1,dimV(1)
                    !$acc loop vector collapse(3)
                    do z = -BD(3),dimX(3)+BD(3)
                        do i=1,BD(2)
                            do x = -BD(1),dimX(1)+BD(1)
                                f(x,dimX(2)+i,z,vx,vy,vz) = in_buffer(x,i,z,vx,vy,vz)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$acc end parallel
    end if
    !$acc end data

end subroutine boundary_apply_buffer_f_y


subroutine boundary_apply_buffer_f_z(f,in_buffer,side,dimX,dimV,BD)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), side
    real(kind=PRC), intent(inout) :: f(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))
    real(kind=PRC), intent(in) :: in_buffer(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),BD(3),dimV(1),dimV(2),dimV(3))

    integer :: i,x,y,vx,vy,vz

    !$acc data present(f) present_or_copyin(in_buffer)
    if (side == 0) then ! lower boundary
        !$acc parallel
        !$acc loop gang collapse(3)
        do vz = 1,dimV(3)
            do vy = 1,dimV(2)
                do vx = 1,dimV(1)
                    !$acc loop vector collapse(3)
                    do i=1,BD(3)
                        do y = -BD(2),dimX(2)+BD(2)
                            do x = -BD(1),dimX(1)+BD(1)
                                f(x,y,-i,vx,vy,vz) = in_buffer(x,y,i,vx,vy,vz)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$acc end parallel
    else if (side == 1) then ! upper boundary
        !$acc parallel
        !$acc loop gang collapse(3)
        do vz = 1,dimV(3)
            do vy = 1,dimV(2)
                do vx = 1,dimV(1)
                    !$acc loop vector collapse(3)
                    do i=1,BD(3)
                        do y = -BD(2),dimX(2)+BD(2)
                            do x = -BD(1),dimX(1)+BD(1)
                                f(x,y,dimX(3)+i,vx,vy,vz) = in_buffer(x,y,i,vx,vy,vz)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$acc end parallel
    end if
    !$acc end data

end subroutine boundary_apply_buffer_f_z


subroutine boundary_apply_boundary_conditions_f_x(f,side,dimX,dimV,BD,bdCond)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), side
    real(kind=PRC), intent(inout) :: f(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))
    character(len=6), intent(in) :: bdCond

    integer :: i, y, z, vx, vy, vz, offset
    real(kind=PRC) :: f_tmp
    offset = 0

    ! boundary conditions for particles (velocity)

    !$acc data present(f)
    if (side == 0) then ! lower boundary
        select case (bdCond(1:1))
        case ('z') ! zero
            if (mod(dimV(1),2) == 1) then
                offset = 1
            end if
            !$acc parallel
            !$acc loop gang collapse(2)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    !$acc loop vector collapse(3)
                    do z = -BD(3),dimX(3)+BD(3)
                        do y = -BD(2),dimX(2)+BD(2)
                            do i=0,BD(1)
                                !$acc loop seq
                                do vx = 1,dimV(1)/2-offset
                                    f_tmp = f(0,y,z,dimV(1)/2+1-vx,vy,vz)
                                    f(-i,y,z,dimV(1)/2+1-vx,vy,vz) = 0.5_PRC*(f(0,y,z,dimV(1)/2+1-vx,vy,vz)+f(0,y,z,dimV(1)/2+vx+offset,vy,vz))
                                    f(-i,y,z,dimV(1)/2+vx+offset,vy,vz) = 0.5_PRC*(f_tmp+f(0,y,z,dimV(1)/2+vx+offset,vy,vz))
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case ('s','m') ! symmetric/mirrored
            !$acc parallel
            !$acc loop gang collapse(3)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    do vx = 1,dimV(1)
                        !$acc loop vector collapse(3)
                        do z = -BD(3),dimX(3)+BD(3)
                            do y = -BD(2),dimX(2)+BD(2)
                                do i=1,BD(1)
                                    f(-i,y,z,vx,vy,vz) = f(i-1,y,z,vx,vy,vz)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case ('a') ! antisymmetric
            !$acc parallel
            !$acc loop gang collapse(3)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    do vx = 1,dimV(1)
                        !$acc loop vector collapse(3)
                        do z = -BD(3),dimX(3)+BD(3)
                            do y = -BD(2),dimX(2)+BD(2)
                                do i=1,BD(1)
                                    f(-i,y,z,vx,vy,vz) = f(i-1,y,z,dimV(1)+1-vx,vy,vz)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case ('q') ! antiperiodic
            !$acc parallel
            !$acc loop gang collapse(2)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    !$acc loop vector collapse(3)
                    do z = -BD(3),dimX(3)+BD(3)
                        do y = -BD(2),dimX(2)+BD(2)
                            do i=1,BD(1)
                                !$acc loop seq
                                do vx = 1,dimV(1)/2
                                    f_tmp = f(-i,y,z,dimV(1)+1-vx,vy,vz)
                                    f(-i,y,z,dimV(1)+1-vx,vy,vz) = f(-i,y,z,vx,vy,vz)
                                    f(-i,y,z,vx,vy,vz) = f_tmp
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        end select

    else if (side == 1) then ! upper boundary
        select case (bdCond(2:2))
        case ('z') ! zero
            if (mod(dimV(1),2) == 1) then
                offset = 1
            end if
            !$acc parallel
            !$acc loop gang collapse(2)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    !$acc loop vector collapse(3)
                    do z = -BD(3),dimX(3)+BD(3)
                        do y = -BD(2),dimX(2)+BD(2)
                            do i=0,BD(1)
                                !$acc loop seq
                                do vx = 1,dimV(1)/2-offset
                                    f_tmp = f(dimX(1),y,z,dimV(1)/2+1-vx,vy,vz)
                                    f(dimX(1)+i,y,z,dimV(1)/2+1-vx,vy,vz) = 0.5_PRC*(f(dimX(1),y,z,dimV(1)/2+1-vx,vy,vz)+f(dimX(1),y,z,dimV(1)/2+vx+offset,vy,vz))
                                    f(dimX(1)+i,y,z,dimV(1)/2+vx+offset,vy,vz) = 0.5_PRC*(f_tmp+f(dimX(1),y,z,dimV(1)/2+vx+offset,vy,vz))
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case ('s','m') ! symmetric/mirrored
            !$acc parallel
            !$acc loop gang collapse(3)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    do vx = 1,dimV(1)
                        !$acc loop vector collapse(3)
                        do z = -BD(3),dimX(3)+BD(3)
                            do y = -BD(2),dimX(2)+BD(2)
                                do i=1,BD(1)
                                    f(dimX(1)+i,y,z,vx,vy,vz) = f(dimX(1)+1-i,y,z,vx,vy,vz)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case ('a') ! antisymmetric
            !$acc parallel
            !$acc loop gang collapse(3)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    do vx = 1,dimV(1)
                        !$acc loop vector collapse(3)
                        do z = -BD(3),dimX(3)+BD(3)
                            do y = -BD(2),dimX(2)+BD(2)
                                do i=1,BD(1)
                                    f(dimX(1)+i,y,z,vx,vy,vz) = f(dimX(1)+1-i,y,z,dimV(1)+1-vx,vy,vz)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case ('q') ! antiperiodic
            !$acc parallel
            !$acc loop gang collapse(2)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    !$acc loop vector collapse(3)
                    do z = -BD(3),dimX(3)+BD(3)
                        do y = -BD(2),dimX(2)+BD(2)
                            do i=1,BD(1)
                                !$acc loop seq
                                do vx = 1,dimV(1)/2
                                    f_tmp = f(dimX(1)+i,y,z,dimV(1)+1-vx,vy,vz)
                                    f(dimX(1)+i,y,z,dimV(1)+1-vx,vy,vz) = f(dimX(1)+i,y,z,vx,vy,vz)
                                    f(dimX(1)+i,y,z,vx,vy,vz) = f_tmp
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        end select
    end if
    !$acc end data

end subroutine boundary_apply_boundary_conditions_f_x


subroutine boundary_apply_boundary_conditions_f_y(f,side,dimX,dimV,BD,bdCond)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), side
    real(kind=PRC), intent(inout) :: f(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))
    character(len=6), intent(in) :: bdCond

    integer :: i, x, z, vx, vy, vz, offset
    real(kind=PRC) :: f_tmp
    offset = 0

    ! boundary conditions for particles (velocity)

    !$acc data present(f)
    if (side == 0) then ! lower boundary
        select case (bdCond(3:3))
        case ('z') ! zero
            if (mod(dimV(2),2) == 1) then
                offset = 1
            end if
            !$acc parallel
            !$acc loop gang collapse(2)
            do vz = 1,dimV(3)
                do vx = 1,dimV(1)
                    !$acc loop vector collapse(3)
                    do z = -BD(3),dimX(3)+BD(3)
                        do i=0,BD(2)
                            do x = -BD(1),dimX(1)+BD(1)
                                !$acc loop seq
                                do vy = 1,dimV(2)/2-offset
                                    f_tmp = f(x,0,z,vx,dimV(2)/2+1-vy,vz)
                                    f(x,-i,z,vx,dimV(2)/2+1-vy,vz) = 0.5_PRC*(f(x,0,z,vx,dimV(2)/2+1-vy,vz)+f(x,0,z,vx,dimV(2)/2+vy+offset,vz))
                                    f(x,-i,z,vx,dimV(2)/2+vy+offset,vz) = 0.5_PRC*(f_tmp+f(x,0,z,vx,dimV(2)/2+vy+offset,vz))
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case ('s','m') ! symmetric/mirrored
            !$acc parallel
            !$acc loop gang collapse(3)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    do vx = 1,dimV(1)
                        !$acc loop vector collapse(3)
                        do z = -BD(3),dimX(3)+BD(3)
                            do i=1,BD(2)
                                do x = -BD(1),dimX(1)+BD(1)
                                    f(x,-i,z,vx,vy,vz) = f(x,i-1,z,vx,vy,vz)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case ('a') ! antisymmetric
            !$acc parallel
            !$acc loop gang collapse(3)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    do vx = 1,dimV(1)
                        !$acc loop vector collapse(3)
                        do z = -BD(3),dimX(3)+BD(3)
                            do i=1,BD(2)
                                do x = -BD(1),dimX(1)+BD(1)
                                    f(x,-i,z,vx,vy,vz) = f(x,i-1,z,vx,dimV(2)+1-vy,vz)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case ('q') ! antiperiodic
            !$acc parallel
            !$acc loop gang collapse(2)
            do vz = 1,dimV(3)
                do vx = 1,dimV(1)
                    !$acc loop vector collapse(3)
                    do z = -BD(3),dimX(3)+BD(3)
                        do i=1,BD(2)
                            do x = -BD(1),dimX(1)+BD(1)
                                !$acc loop seq
                                do vy = 1,dimV(2)/2
                                    f_tmp = f(x,-i,z,vx,dimV(2)+1-vy,vz)
                                    f(x,-i,z,vx,dimV(2)+1-vy,vz) = f(x,-i,z,vx,vy,vz)
                                    f(x,-i,z,vx,vy,vz) = f_tmp
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        end select

    else if (side == 1) then ! upper boundary
        select case (bdCond(4:4))
        case ('z') ! zero
            if (mod(dimV(2),2) == 1) then
                offset = 1
            end if
            !$acc parallel
            !$acc loop gang collapse(2)
            do vz = 1,dimV(3)
                do vx = 1,dimV(1)
                    !$acc loop vector collapse(3)
                    do z = -BD(3),dimX(3)+BD(3)
                        do i=0,BD(2)
                            do x = -BD(1),dimX(1)+BD(1)
                                !$acc loop seq
                                do vy = 1,dimV(2)/2-offset
                                    f_tmp = f(x,dimX(2),z,vx,dimV(2)/2+1-vy,vz)
                                    f(x,dimX(2)+i,z,vx,dimV(2)/2+1-vy,vz) = 0.5_PRC*(f(x,dimX(2),z,vx,dimV(2)/2+1-vy,vz)+f(x,dimX(2),z,vx,dimV(2)/2+vy+offset,vz))
                                    f(x,dimX(2)+i,z,vx,dimV(2)/2+vy+offset,vz) = 0.5_PRC*(f_tmp+f(x,dimX(2),z,vx,dimV(2)/2+vy+offset,vz))
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case ('s','m') ! symmetric/mirrored
            !$acc parallel
            !$acc loop gang collapse(3)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    do vx = 1,dimV(1)
                        !$acc loop vector collapse(3)
                        do z = -BD(3),dimX(3)+BD(3)
                            do i=1,BD(2)
                                do x = -BD(1),dimX(1)+BD(1)
                                    f(x,dimX(2)+i,z,vx,vy,vz) = f(x,dimX(2)+1-i,z,vx,vy,vz)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case ('a') ! antisymmetric
            !$acc parallel
            !$acc loop gang collapse(3)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    do vx = 1,dimV(1)
                        !$acc loop vector collapse(3)
                        do z = -BD(3),dimX(3)+BD(3)
                            do i=1,BD(2)
                                do x = -BD(1),dimX(1)+BD(1)
                                    f(x,dimX(2)+i,z,vx,vy,vz) = f(x,dimX(2)+1-i,z,vx,dimV(2)+1-vy,vz)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case ('q') ! antiperiodic
            !$acc parallel
            !$acc loop gang collapse(2)
            do vz = 1,dimV(3)
                do vx = 1,dimV(1)
                    !$acc loop vector collapse(3)
                    do z = -BD(3),dimX(3)+BD(3)
                        do i=1,BD(2)
                            do x = -BD(1),dimX(1)+BD(1)
                                !$acc loop seq
                                do vy = 1,dimV(2)/2
                                    f_tmp = f(x,dimX(2)+i,z,vx,dimV(2)+1-vy,vz)
                                    f(x,dimX(2)+i,z,vx,dimV(2)+1-vy,vz) = f(x,dimX(2)+i,z,vx,vy,vz)
                                    f(x,dimX(2)+i,z,vx,vy,vz) = f_tmp
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        end select
    end if
    !$acc end data

end subroutine boundary_apply_boundary_conditions_f_y


subroutine boundary_apply_boundary_conditions_f_z(f,side,dimX,dimV,BD,bdCond)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), side
    real(kind=PRC), intent(inout) :: f(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3))
    character(len=6), intent(in) :: bdCond

    integer :: i, x, y, vx, vy, vz, offset
    real(kind=PRC) :: f_tmp
    offset = 0

    ! boundary conditions for particles (velocity)

    !$acc data present(f)
    if (side == 0) then ! lower boundary
        select case (bdCond(5:5))
        case ('z') ! zero
            if (mod(dimV(3),2) == 1) then
                offset = 1
            end if
            !$acc parallel
            !$acc loop gang collapse(2)
            do vy = 1,dimV(3)
                do vx = 1,dimV(1)
                    !$acc loop vector collapse(3)
                    do i=0,BD(3)
                        do y = -BD(2),dimX(2)+BD(2)
                            do x = -BD(1),dimX(1)+BD(1)
                                !$acc loop seq
                                do vz = 1,dimV(3)/2-offset
                                    f_tmp = f(x,y,0,vx,vy,dimV(3)/2+1-vz)
                                    f(x,y,-i,vx,vy,dimV(3)/2+1-vz) = 0.5_PRC*(f(x,y,0,vx,vy,dimV(3)/2+1-vz)+f(x,y,0,vx,vy,dimV(3)/2+vz+offset))
                                    f(x,y,-i,vx,vy,dimV(3)/2+vz+offset) = 0.5_PRC*(f_tmp+f(x,y,0,vx,vy,dimV(3)/2+vz+offset))
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case ('s','m') ! symmetric/mirrored
            !$acc parallel
            !$acc loop gang collapse(3)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    do vx = 1,dimV(1)
                        !$acc loop vector collapse(3)
                        do i=1,BD(3)
                            do y = -BD(2),dimX(2)+BD(2)
                                do x = -BD(1),dimX(1)+BD(1)
                                    f(x,y,-i,vx,vy,vz) = f(x,y,i-1,vx,vy,vz)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case ('a') ! antisymmetric
            !$acc parallel
            !$acc loop gang collapse(3)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    do vx = 1,dimV(1)
                        !$acc loop vector collapse(3)
                        do i=1,BD(3)
                            do y = -BD(2),dimX(2)+BD(2)
                                do x = -BD(1),dimX(1)+BD(1)
                                    f(x,y,-i,vx,vy,vz) = f(x,y,i-1,vx,vy,dimV(3)+1-vz)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case ('q') ! antiperiodic
            !$acc parallel
            !$acc loop gang collapse(2)
            do vy = 1,dimV(3)
                do vx = 1,dimV(1)
                    !$acc loop vector collapse(3)
                    do i=1,BD(3)
                        do y = -BD(2),dimX(2)+BD(2)
                            do x = -BD(1),dimX(1)+BD(1)
                                !$acc loop seq
                                do vz = 1,dimV(3)/2
                                    f_tmp = f(x,y,-i,vx,vy,dimV(3)+1-vz)
                                    f(x,y,-i,vx,vy,dimV(3)+1-vz) = f(x,y,-i,vx,vy,vz)
                                    f(x,y,-i,vx,vy,vz) = f_tmp
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        end select

    else if (side == 1) then ! upper boundary
        select case (bdCond(6:6))
        case ('z') ! zero
            if (mod(dimV(3),2) == 1) then
                offset = 1
            end if
            !$acc parallel
            !$acc loop gang collapse(2)
            do vy = 1,dimV(3)
                do vx = 1,dimV(1)
                    !$acc loop vector collapse(3)
                    do i=0,BD(3)
                        do y = -BD(2),dimX(2)+BD(2)
                            do x = -BD(1),dimX(1)+BD(1)
                                !$acc loop seq
                                do vz = 1,dimV(3)/2-offset
                                    f_tmp = f(x,y,dimX(3),vx,vy,dimV(3)/2+1-vz)
                                    f(x,y,dimX(3)+i,vx,vy,dimV(3)/2+1-vz) = 0.5_PRC*(f(x,y,dimX(3),vx,vy,dimV(3)/2+1-vz)+f(x,y,dimX(3),vx,vy,dimV(3)/2+vz+offset))
                                    f(x,y,dimX(3)+i,vx,vy,dimV(3)/2+vz+offset) = 0.5_PRC*(f_tmp+f(x,y,dimX(3),vx,vy,dimV(3)/2+vz+offset))
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case ('s','m') ! symmetric/mirrored
            !$acc parallel
            !$acc loop gang collapse(3)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    do vx = 1,dimV(1)
                        !$acc loop vector collapse(3)
                        do i=1,BD(3)
                            do y = -BD(2),dimX(2)+BD(2)
                                do x = -BD(1),dimX(1)+BD(1)
                                    f(x,y,dimX(3)+i,vx,vy,vz) = f(x,y,dimX(3)+1-i,vx,vy,vz)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case ('a') ! antisymmetric
            !$acc parallel
            !$acc loop gang collapse(3)
            do vz = 1,dimV(3)
                do vy = 1,dimV(2)
                    do vx = 1,dimV(1)
                        !$acc loop vector collapse(3)
                        do i=1,BD(3)
                            do y = -BD(2),dimX(2)+BD(2)
                                do x = -BD(1),dimX(1)+BD(1)
                                    f(x,y,dimX(3)+i,vx,vy,vz) = f(x,y,dimX(3)+1-i,vx,vy,dimV(3)+1-vz)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        case ('q') ! antiperiodic
            !$acc parallel
            !$acc loop gang collapse(2)
            do vy = 1,dimV(3)
                do vx = 1,dimV(1)
                    !$acc loop vector collapse(3)
                    do i=1,BD(3)
                        do y = -BD(2),dimX(2)+BD(2)
                            do x = -BD(1),dimX(1)+BD(1)
                                !$acc loop seq
                                do vz = 1,dimV(3)/2
                                    f_tmp = f(x,y,dimX(3)+i,vx,vy,dimV(3)+1-vz)
                                    f(x,y,dimX(3)+i,vx,vy,dimV(3)+1-vz) = f(x,y,dimX(3)+i,vx,vy,vz)
                                    f(x,y,dimX(3)+i,vx,vy,vz) = f_tmp
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end parallel
        end select
    end if
    !$acc end data

end subroutine boundary_apply_boundary_conditions_f_z

