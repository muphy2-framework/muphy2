
! This file is part of the muphyII multiphysics plasma simulation project.
! 
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at https://mozilla.org/MPL/2.0/.


subroutine utilities_copy_array(target_array, source_array, array_length)
    use definitions
    implicit none

    integer, intent(in) :: array_length
    real(kind=PRC), intent(out) :: target_array(array_length)
    real(kind=PRC), intent(in) :: source_array(array_length)

    !$acc data present_or_copyout(target_array) present_or_copyin(source_array)
    !$acc kernels
    target_array = source_array
    !$acc end kernels
    !$acc end data

end subroutine utilities_copy_array

subroutine utilities_copy_array_no_gpu(target_array, source_array, array_length)
    use definitions
    implicit none

    integer, intent(in) :: array_length
    real(kind=PRC), intent(inout) :: target_array(array_length)
    real(kind=PRC), intent(in) :: source_array(array_length)

    target_array = source_array

end subroutine utilities_copy_array_no_gpu


subroutine utilities_add_arrays(target_array, array1, array2, factor1, factor2, array_length)
    use definitions
    implicit none

    integer, intent(in) :: array_length
    real(kind=PRC), intent(inout) :: target_array(array_length)
    real(kind=PRC), intent(in) :: array1(array_length), array2(array_length), factor1, factor2

    !$acc data present_or_copyout(target_array) present_or_copyin(array1, array2) copyin(factor1, factor2)
    !$acc kernels
    target_array = factor1*array1 + factor2*array2
    !$acc end kernels
    !$acc end data

end subroutine utilities_add_arrays


subroutine utilities_array_add_scalar(array, scalar, array_length)
    use definitions
    implicit none

    integer, intent(in) :: array_length
    real(kind=PRC), intent(inout) :: array(array_length)
    real(kind=PRC), intent(in) :: scalar

    !$acc data present_or_copy(array) copyin(scalar)
    !$acc kernels
    array = array + scalar
    !$acc end kernels
    !$acc end data

end subroutine utilities_array_add_scalar


subroutine utilities_array_fill(array, fill_value, array_length)
    use definitions
    implicit none

    integer, intent(in) :: array_length
    real(kind=PRC), intent(out) :: array(array_length)
    real(kind=PRC), intent(in) :: fill_value

    !$acc data present_or_copy(array) copyin(fill_value)
    !$acc kernels
    array = fill_value
    !$acc end kernels
    !$acc end data

end subroutine utilities_array_fill


subroutine utilities_vector_magnitude(magnitude, vector, dimX, BD)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(inout) :: magnitude(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3))
    real(kind=PRC), intent(in) :: vector(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3)

    integer :: x, y, z

    !$acc data present(magnitude, vector)
    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            !$acc loop vector
            do x = 0,dimX(1)
                magnitude(x,y,z) = sqrt(vector(x,y,z,1)*vector(x,y,z,1) + &
                                        vector(x,y,z,2)*vector(x,y,z,2) + &
                                        vector(x,y,z,3)*vector(x,y,z,3))
            enddo
        enddo
    enddo
    !$acc end parallel
    !$acc end data

end subroutine utilities_vector_magnitude

