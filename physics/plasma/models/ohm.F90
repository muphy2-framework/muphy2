
! This file is part of the muphyII multiphysics plasma simulation project.
! 
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at https://mozilla.org/MPL/2.0/.


! TODO this model is work in progress and not functional

subroutine ohm_calc_E_fc_ideal(E_fc, B_ec, u_cc, dimX, BD)
  use definitions
  implicit none

  integer, intent(in) :: dimX(1:3), BD(3)
  real(kind=PRC), intent(in),  dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u_cc, B_ec
  real(kind=PRC), intent(out), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: E_fc

  integer :: x,y,z
  ! E = - u x B  @fc
  do z = 0,dimX(3)
     do y = 0,dimX(2)
        do x = 0,dimX(1)
           E_fc(x,y,z,1) = -0.25_PRC*( &
                (u_cc(x,y,z,2) + u_cc(x-1,y,z,2)) * (B_ec(x,y,z,3) + B_ec(x,y+1,z,3)) - &
                (u_cc(x,y,z,3) + u_cc(x-1,y,z,3)) * (B_ec(x,y,z,2) + B_ec(x,y,z+1,2)) )
           E_fc(x,y,z,2) = -0.25_PRC*( &              
                (u_cc(x,y,z,3) + u_cc(x,y-1,z,3)) * (B_ec(x,y,z,1) + B_ec(x,y,z+1,1)) - &
                (u_cc(x,y,z,1) + u_cc(x,y-1,z,1)) * (B_ec(x,y,z,3) + B_ec(x+1,y,z,3)) )
           E_fc(x,y,z,3) = -0.25_PRC*( &             
                (u_cc(x,y,z,1) + u_cc(x,y,z-1,1)) * (B_ec(x,y,z,2) + B_ec(x+1,y,z,2)) - &
                (u_cc(x,y,z,2) + u_cc(x,y,z-1,2)) * (B_ec(x,y,z,1) + B_ec(x,y+1,z,1)) )
        enddo
     enddo
  enddo

end subroutine ohm_calc_E_fc_ideal


subroutine ohm_calc_E_fc_hallterm(E_fc, B_ec, J_fc, rho_cc, dimX, BD)
  use definitions
  implicit none

  integer, intent(in) :: dimX(1:3), BD(3)
  real(kind=PRC), intent(in),    dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: B_ec, J_fc
  real(kind=PRC), intent(in),    dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3))   :: rho_cc
  real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: E_fc

  integer :: x,y,z
  ! E += 1/rho J x B  @fc
  do z = 0,dimX(3)
     do y = 0,dimX(2)
        do x = 0,dimX(1)
           E_fc(x,y,z,1) = E_fc(x,y,z,1) + &
                (J_fc(x,y,z,2) * (B_ec(x,y,z,3) + B_ec(x,y+1,z,3)) - &
                J_fc(x,y,z,3) * (B_ec(x,y,z,2) + B_ec(x,y,z+1,2))) / (rho_cc(x,y,z) + rho_cc(x-1,y,z)) 
           E_fc(x,y,z,2) = E_fc(x,y,z,2) + &
                (J_fc(x,y,z,3) * (B_ec(x,y,z,1) + B_ec(x,y,z+1,1)) - &
                J_fc(x,y,z,1) * (B_ec(x,y,z,3) + B_ec(x+1,y,z,3))) / (rho_cc(x,y,z) + rho_cc(x,y-1,z)) 
           E_fc(x,y,z,3) = E_fc(x,y,z,3) + &
                (J_fc(x,y,z,1) * (B_ec(x,y,z,2) + B_ec(x+1,y,z,2)) - &
                J_fc(x,y,z,2) * (B_ec(x,y,z,1) + B_ec(x,y+1,z,1))) / (rho_cc(x,y,z) + rho_cc(x,y,z-1)) 
        enddo
     enddo
  enddo
  !@TODO: rho~0
end subroutine ohm_calc_E_fc_hallterm


subroutine ohm_calc_E_fc_pressureterm(E_fc, P_cc, rho_cc, dx, dimX, BD)
  use definitions
  implicit none

  integer, intent(in) :: dimX(1:3), BD(3)
  real(kind=PRC), intent(in) :: dx(3)
  real(kind=PRC), intent(in),    dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3))   :: P_cc, rho_cc
  real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: E_fc

  integer :: x,y,z
  ! E -= 1/rho \grad P_e  @fc
  ! where fc \grad p is calculated directly from cc-values
  do z = 0,dimX(3)
     do y = 0,dimX(2)
        do x = 0,dimX(1)
           E_fc(x,y,z,1) = E_fc(x,y,z,1) - &
                (P_cc(x,y,z) - P_cc(x-1,y,z))/dx(1) / (0.5_PRC*(rho_cc(x,y,z) + rho_cc(x-1,y,z)))
           E_fc(x,y,z,2) = E_fc(x,y,z,2) - &
                (P_cc(x,y,z) - P_cc(x,y-1,z))/dx(2) / (0.5_PRC*(rho_cc(x,y,z) + rho_cc(x,y-1,z)))
           E_fc(x,y,z,3) = E_fc(x,y,z,3) - &
                (P_cc(x,y,z) - P_cc(x,y,z-1))/dx(3) / (0.5_PRC*(rho_cc(x,y,z) + rho_cc(x,y,z-1)))
        enddo
     enddo
  enddo
  !@TODO: rho~0

end subroutine ohm_calc_E_fc_pressureterm


subroutine ohm_calc_E_fc_resistivityterm(E_fc, J_fc, resistivity, dimX, BD)
  use definitions
  implicit none

  integer, intent(in) :: dimX(1:3), BD(3)
  real(kind=PRC), intent(in) :: resistivity
  real(kind=PRC), intent(in),    dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: J_fc
  real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: E_fc

  E_fc = E_fc + J_fc * resistivity

  !@TODO

end subroutine ohm_calc_E_fc_resistivityterm


subroutine ohm_step_B(E, B, dimX, BD, dx, dt)
  use definitions
  implicit none

  integer, intent(in)           :: dimX(1:3), BD(3)
  real(kind=PRC), intent(in)    :: dx(1:3), dt
  real(kind=PRC), intent(in)    :: E(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3)
  real(kind=PRC), intent(inout) :: B(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3)

  ! dB/dt = - curl E
  B(0:dimX(1), 0:dimX(2), 0:dimX(3), 1) = B(0:dimX(1), 0:dimX(2), 0:dimX(3),1) - dt * (&
       ( E(0:dimX(1), 0:dimX(2), 0:dimX(3), 3) - E(0:dimX(1), -1:dimX(2)-1,  0:dimX(3),   3) )/dx(2) - &
       ( E(0:dimX(1), 0:dimX(2), 0:dimX(3), 2) - E(0:dimX(1),  0:dimX(2),   -1:dimX(3)-1, 2) )/dx(3) )
  B(0:dimX(1), 0:dimX(2), 0:dimX(3), 2) = B(0:dimX(1), 0:dimX(2), 0:dimX(3),2) - dt * ( &
       ( E(0:dimX(1), 0:dimX(2), 0:dimX(3), 1) - E( 0:dimX(1),   0:dimX(2), -1:dimX(3)-1, 1) )/dx(3) - &
       ( E(0:dimX(1), 0:dimX(2), 0:dimX(3), 3) - E(-1:dimX(1)-1, 0:dimX(2),  0:dimX(3),   3) )/dx(1) )
  B(0:dimX(1), 0:dimX(2), 0:dimX(3), 3) = B(0:dimX(1), 0:dimX(2), 0:dimX(3),3) - dt * ( &
       ( E(0:dimX(1), 0:dimX(2), 0:dimX(3), 2) - E(-1:dimX(1)-1,  0:dimX(2),   0:dimX(3), 2) )/dx(1) - &
       ( E(0:dimX(1), 0:dimX(2), 0:dimX(3), 1) - E( 0:dimX(1),   -1:dimX(2)-1, 0:dimX(3), 1) )/dx(2) )

end subroutine ohm_step_B


subroutine ohm_calc_J(J, B, dimX, BD, dx)
  use definitions
  implicit none

  integer, intent(in) :: dimX(1:3), BD(3)
  real(kind=PRC), intent(in)  :: dx(1:3)
  real(kind=PRC), intent(out), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: J
  real(kind=PRC), intent(in),  dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: B

  ! J = curl B
  J(0:dimX(1), 0:dimX(2), 0:dimX(3), 1) = &
       ( B(0:dimX(1), 1:dimX(2)+1, 0:dimX(3),   3) - B(0:dimX(1), 0:dimX(2), 0:dimX(3), 3) )/dx(2) - &
       ( B(0:dimX(1), 0:dimX(2),   1:dimX(3)+1, 2) - B(0:dimX(1), 0:dimX(2), 0:dimX(3), 2) )/dx(3)
  J(0:dimX(1), 0:dimX(2), 0:dimX(3), 2) = &
       ( B(0:dimX(1),   0:dimX(2), 1:dimX(3)+1, 1) - B(0:dimX(1), 0:dimX(2), 0:dimX(3), 1) )/dx(3) - &
       ( B(1:dimX(1)+1, 0:dimX(2), 0:dimX(3),   3) - B(0:dimX(1), 0:dimX(2), 0:dimX(3), 3) )/dx(1)
  J(0:dimX(1), 0:dimX(2), 0:dimX(3), 3) = &
       ( B(1:dimX(1)+1, 0:dimX(2),   0:dimX(3), 2) - B(0:dimX(1), 0:dimX(2), 0:dimX(3), 2) )/dx(1) - &
       ( B(0:dimX(1),   1:dimX(2)+1, 0:dimX(3), 1) - B(0:dimX(1), 0:dimX(2), 0:dimX(3), 1) )/dx(2)

end subroutine ohm_calc_J

