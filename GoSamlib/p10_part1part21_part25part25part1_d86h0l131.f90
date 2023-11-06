module     p10_part1part21_part25part25part1_d86h0l131
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p10_part1part &
   ! &21_part25part25part1/helicity0d86h0l131.f90
   ! generator: buildfortran_tn3.py
   use p10_part1part21_part25part25part1_config, only: ki
   use p10_part1part21_part25part25part1_util, only: cond_t, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, parameter :: ninjaidxt3mu0 = 0
   integer, parameter :: ninjaidxt2mu0 = 1
   integer, parameter :: ninjaidxt1mu0 = 2
   integer, parameter :: ninjaidxt1mu2 = 3
   integer, parameter :: ninjaidxt0mu0 = 4
   integer, parameter :: ninjaidxt0mu2 = 5
   public :: numerator_t3
contains
!---#[ subroutine brack_31:
   pure subroutine brack_31(ninjaA, ninjaE3, ninjaE4, ninjaP, brack)
      use p10_part1part21_part25part25part1_model
      use p10_part1part21_part25part25part1_kinematics
      use p10_part1part21_part25part25part1_color
      use p10_part1part21_part25part25part1_abbrevd86h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(17) :: acd86
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd86(1)=dotproduct(ninjaA,ninjaE3)
      acd86(2)=dotproduct(ninjaE3,spvak5e2)
      acd86(3)=dotproduct(ninjaE3,spvae2k1)
      acd86(4)=abb86(9)
      acd86(5)=dotproduct(ninjaE3,spval3k1)
      acd86(6)=dotproduct(ninjaE3,spvae2l3)
      acd86(7)=abb86(7)
      acd86(8)=dotproduct(ninjaE3,spvak2k1)
      acd86(9)=dotproduct(ninjaE3,spvae2k2)
      acd86(10)=abb86(18)
      acd86(11)=dotproduct(ninjaE3,spvak5l3)
      acd86(12)=dotproduct(ninjaE3,spval3e2)
      acd86(13)=dotproduct(ninjaE3,spvak5k2)
      acd86(14)=dotproduct(ninjaE3,spvak2e2)
      acd86(15)=acd86(4)*acd86(1)*acd86(2)
      acd86(16)=acd86(12)*acd86(11)*acd86(7)
      acd86(17)=acd86(14)*acd86(13)*acd86(10)
      acd86(15)=acd86(17)+acd86(16)+2.0_ki*acd86(15)
      acd86(15)=acd86(3)*acd86(15)
      acd86(16)=acd86(6)*acd86(5)*acd86(7)
      acd86(17)=acd86(9)*acd86(8)*acd86(10)
      acd86(16)=acd86(17)+acd86(16)
      acd86(16)=acd86(2)*acd86(16)
      acd86(15)=acd86(15)+acd86(16)
      brack(ninjaidxt3mu0)=0.0_ki
      brack(ninjaidxt2mu0)=acd86(15)
   end subroutine brack_31
!---#] subroutine brack_31:
!---#[ subroutine brack_32:
   pure subroutine brack_32(ninjaA, ninjaE3, ninjaE4, ninjaP, brack)
      use p10_part1part21_part25part25part1_model
      use p10_part1part21_part25part25part1_kinematics
      use p10_part1part21_part25part25part1_color
      use p10_part1part21_part25part25part1_abbrevd86h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(89) :: acd86
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd86(1)=dotproduct(ninjaA,ninjaE3)
      acd86(2)=dotproduct(ninjaE3,spvak5e2)
      acd86(3)=dotproduct(ninjaE4,spvae2k1)
      acd86(4)=abb86(9)
      acd86(5)=dotproduct(ninjaE3,spvae2k1)
      acd86(6)=dotproduct(ninjaE4,spvak5e2)
      acd86(7)=dotproduct(ninjaA,ninjaE4)
      acd86(8)=dotproduct(ninjaA,spvak5e2)
      acd86(9)=dotproduct(ninjaA,spvae2k1)
      acd86(10)=dotproduct(ninjaE3,spval3k1)
      acd86(11)=dotproduct(ninjaE3,spvae2l3)
      acd86(12)=abb86(7)
      acd86(13)=dotproduct(ninjaE4,spvae2l3)
      acd86(14)=dotproduct(ninjaE4,spval3k1)
      acd86(15)=dotproduct(ninjaE3,spvae2k2)
      acd86(16)=dotproduct(ninjaE4,spvak2k1)
      acd86(17)=abb86(18)
      acd86(18)=dotproduct(ninjaE3,spvak2k1)
      acd86(19)=dotproduct(ninjaE4,spvae2k2)
      acd86(20)=abb86(25)
      acd86(21)=dotproduct(ninjaE3,spvak5l3)
      acd86(22)=dotproduct(ninjaE4,spval3e2)
      acd86(23)=dotproduct(ninjaE3,spval3e2)
      acd86(24)=dotproduct(ninjaE4,spvak5l3)
      acd86(25)=dotproduct(ninjaE3,spvak2e2)
      acd86(26)=dotproduct(ninjaE4,spvak5k2)
      acd86(27)=dotproduct(ninjaE3,spvak5k2)
      acd86(28)=dotproduct(ninjaE4,spvak2e2)
      acd86(29)=abb86(23)
      acd86(30)=dotproduct(ninjaA,ninjaA)
      acd86(31)=dotproduct(ninjaA,spval3k1)
      acd86(32)=dotproduct(ninjaA,spvae2l3)
      acd86(33)=dotproduct(ninjaA,spvak5l3)
      acd86(34)=dotproduct(ninjaA,spval3e2)
      acd86(35)=dotproduct(ninjaA,spvae2k2)
      acd86(36)=dotproduct(ninjaA,spvak2e2)
      acd86(37)=dotproduct(ninjaA,spvak5k2)
      acd86(38)=dotproduct(ninjaA,spvak2k1)
      acd86(39)=abb86(28)
      acd86(40)=abb86(14)
      acd86(41)=abb86(29)
      acd86(42)=abb86(33)
      acd86(43)=abb86(21)
      acd86(44)=abb86(20)
      acd86(45)=abb86(16)
      acd86(46)=abb86(32)
      acd86(47)=abb86(22)
      acd86(48)=abb86(24)
      acd86(49)=abb86(27)
      acd86(50)=abb86(11)
      acd86(51)=abb86(31)
      acd86(52)=abb86(8)
      acd86(53)=abb86(26)
      acd86(54)=abb86(30)
      acd86(55)=abb86(12)
      acd86(56)=abb86(13)
      acd86(57)=abb86(17)
      acd86(58)=abb86(34)
      acd86(59)=dotproduct(ninjaE3,spval4e2)
      acd86(60)=abb86(19)
      acd86(61)=dotproduct(ninjaE3,spvae2l4)
      acd86(62)=abb86(37)
      acd86(63)=2.0_ki*acd86(1)
      acd86(64)=acd86(63)*acd86(3)
      acd86(64)=acd86(64)+acd86(9)
      acd86(64)=acd86(64)*acd86(4)
      acd86(65)=2.0_ki*acd86(4)
      acd86(65)=acd86(7)*acd86(65)*acd86(5)
      acd86(64)=acd86(65)+acd86(64)+acd86(20)
      acd86(65)=acd86(18)*acd86(19)
      acd86(66)=acd86(15)*acd86(16)
      acd86(65)=acd86(65)+acd86(66)
      acd86(66)=acd86(17)*acd86(65)
      acd86(67)=acd86(10)*acd86(13)
      acd86(68)=acd86(11)*acd86(14)
      acd86(67)=acd86(67)+acd86(68)
      acd86(68)=acd86(12)*acd86(67)
      acd86(66)=acd86(68)+acd86(66)+acd86(64)
      acd86(66)=acd86(2)*acd86(66)
      acd86(68)=acd86(27)*acd86(28)
      acd86(69)=acd86(25)*acd86(26)
      acd86(68)=acd86(68)+acd86(69)
      acd86(69)=acd86(17)*acd86(68)
      acd86(70)=acd86(21)*acd86(22)
      acd86(71)=acd86(23)*acd86(24)
      acd86(70)=acd86(70)+acd86(71)
      acd86(71)=acd86(12)*acd86(70)
      acd86(72)=acd86(63)*acd86(6)
      acd86(72)=acd86(72)+acd86(8)
      acd86(72)=acd86(72)*acd86(4)
      acd86(72)=acd86(72)+acd86(29)
      acd86(69)=acd86(71)+acd86(69)+acd86(72)
      acd86(69)=acd86(5)*acd86(69)
      acd86(71)=acd86(15)*acd86(18)
      acd86(73)=acd86(71)*acd86(6)
      acd86(74)=acd86(25)*acd86(27)
      acd86(75)=acd86(74)*acd86(3)
      acd86(73)=acd86(73)+acd86(75)
      acd86(75)=acd86(17)*acd86(73)
      acd86(76)=acd86(11)*acd86(10)
      acd86(77)=acd86(76)*acd86(6)
      acd86(78)=acd86(23)*acd86(21)
      acd86(79)=acd86(78)*acd86(3)
      acd86(77)=acd86(77)+acd86(79)
      acd86(79)=acd86(12)*acd86(77)
      acd86(66)=acd86(66)+acd86(69)+acd86(75)+acd86(79)
      acd86(69)=acd86(15)*acd86(42)
      acd86(75)=acd86(11)*acd86(40)
      acd86(79)=acd86(63)*acd86(20)
      acd86(69)=acd86(79)+acd86(69)+acd86(75)
      acd86(75)=acd86(4)*acd86(9)
      acd86(79)=acd86(63)*acd86(75)
      acd86(80)=acd86(18)*acd86(35)
      acd86(81)=acd86(15)*acd86(38)
      acd86(80)=acd86(80)+acd86(81)
      acd86(81)=acd86(17)*acd86(80)
      acd86(82)=acd86(10)*acd86(32)
      acd86(83)=acd86(11)*acd86(31)
      acd86(82)=acd86(82)+acd86(83)
      acd86(83)=acd86(12)*acd86(82)
      acd86(84)=ninjaP+acd86(30)
      acd86(84)=acd86(4)*acd86(84)
      acd86(84)=acd86(41)+acd86(84)
      acd86(84)=acd86(5)*acd86(84)
      acd86(79)=acd86(84)+acd86(83)+acd86(81)+acd86(79)+acd86(69)
      acd86(79)=acd86(2)*acd86(79)
      acd86(81)=acd86(4)*acd86(8)
      acd86(81)=acd86(81)+acd86(29)
      acd86(83)=acd86(63)*acd86(81)
      acd86(84)=acd86(25)*acd86(45)
      acd86(85)=acd86(23)*acd86(44)
      acd86(83)=acd86(84)+acd86(85)+acd86(83)
      acd86(84)=acd86(27)*acd86(36)
      acd86(85)=acd86(25)*acd86(37)
      acd86(84)=acd86(84)+acd86(85)
      acd86(85)=acd86(17)*acd86(84)
      acd86(86)=acd86(21)*acd86(34)
      acd86(87)=acd86(23)*acd86(33)
      acd86(86)=acd86(86)+acd86(87)
      acd86(87)=acd86(12)*acd86(86)
      acd86(85)=acd86(87)+acd86(85)+acd86(83)
      acd86(85)=acd86(5)*acd86(85)
      acd86(87)=acd86(9)*acd86(74)
      acd86(88)=acd86(8)*acd86(71)
      acd86(87)=acd86(87)+acd86(88)
      acd86(87)=acd86(17)*acd86(87)
      acd86(88)=acd86(9)*acd86(78)
      acd86(89)=acd86(8)*acd86(76)
      acd86(88)=acd86(88)+acd86(89)
      acd86(88)=acd86(12)*acd86(88)
      acd86(74)=acd86(47)*acd86(74)
      acd86(78)=acd86(43)*acd86(78)
      acd86(71)=acd86(46)*acd86(71)
      acd86(76)=acd86(39)*acd86(76)
      acd86(71)=acd86(79)+acd86(85)+acd86(88)+acd86(87)+acd86(76)+acd86(71)+acd&
      &86(74)+acd86(78)
      acd86(74)=acd86(2)*acd86(5)*acd86(4)
      acd86(68)=ninjaP*acd86(68)
      acd86(76)=acd86(36)*acd86(37)
      acd86(68)=acd86(76)+acd86(68)
      acd86(68)=acd86(17)*acd86(68)
      acd86(70)=ninjaP*acd86(70)
      acd86(76)=acd86(34)*acd86(33)
      acd86(70)=acd86(76)+acd86(70)
      acd86(70)=acd86(12)*acd86(70)
      acd86(76)=acd86(30)*acd86(81)
      acd86(78)=acd86(36)*acd86(45)
      acd86(79)=acd86(34)*acd86(44)
      acd86(81)=acd86(8)*acd86(41)
      acd86(72)=ninjaP*acd86(72)
      acd86(68)=acd86(70)+acd86(68)+acd86(72)+acd86(81)+acd86(79)+acd86(53)+acd&
      &86(78)+acd86(76)
      acd86(68)=acd86(5)*acd86(68)
      acd86(64)=ninjaP*acd86(64)
      acd86(65)=ninjaP*acd86(65)
      acd86(70)=acd86(35)*acd86(38)
      acd86(65)=acd86(70)+acd86(65)
      acd86(65)=acd86(17)*acd86(65)
      acd86(67)=ninjaP*acd86(67)
      acd86(70)=acd86(32)*acd86(31)
      acd86(67)=acd86(70)+acd86(67)
      acd86(67)=acd86(12)*acd86(67)
      acd86(70)=acd86(75)+acd86(20)
      acd86(70)=acd86(30)*acd86(70)
      acd86(72)=acd86(35)*acd86(42)
      acd86(75)=acd86(32)*acd86(40)
      acd86(76)=acd86(9)*acd86(41)
      acd86(64)=acd86(67)+acd86(65)+acd86(76)+acd86(75)+acd86(51)+acd86(72)+acd&
      &86(64)+acd86(70)
      acd86(64)=acd86(2)*acd86(64)
      acd86(65)=acd86(9)*acd86(84)
      acd86(67)=acd86(8)*acd86(80)
      acd86(70)=ninjaP*acd86(73)
      acd86(65)=acd86(70)+acd86(65)+acd86(67)
      acd86(65)=acd86(17)*acd86(65)
      acd86(67)=acd86(9)*acd86(86)
      acd86(70)=acd86(8)*acd86(82)
      acd86(72)=ninjaP*acd86(77)
      acd86(67)=acd86(72)+acd86(67)+acd86(70)
      acd86(67)=acd86(12)*acd86(67)
      acd86(70)=acd86(9)*acd86(83)
      acd86(69)=acd86(8)*acd86(69)
      acd86(72)=acd86(61)*acd86(62)
      acd86(73)=acd86(59)*acd86(60)
      acd86(75)=acd86(36)*acd86(47)
      acd86(75)=acd86(57)+acd86(75)
      acd86(75)=acd86(27)*acd86(75)
      acd86(76)=acd86(34)*acd86(43)
      acd86(76)=acd86(52)+acd86(76)
      acd86(76)=acd86(21)*acd86(76)
      acd86(77)=acd86(35)*acd86(46)
      acd86(77)=acd86(58)+acd86(77)
      acd86(77)=acd86(18)*acd86(77)
      acd86(78)=acd86(32)*acd86(39)
      acd86(78)=acd86(49)+acd86(78)
      acd86(78)=acd86(10)*acd86(78)
      acd86(79)=acd86(37)*acd86(47)
      acd86(79)=acd86(56)+acd86(79)
      acd86(79)=acd86(25)*acd86(79)
      acd86(80)=acd86(33)*acd86(43)
      acd86(80)=acd86(54)+acd86(80)
      acd86(80)=acd86(23)*acd86(80)
      acd86(81)=acd86(38)*acd86(46)
      acd86(81)=acd86(55)+acd86(81)
      acd86(81)=acd86(15)*acd86(81)
      acd86(82)=acd86(31)*acd86(39)
      acd86(82)=acd86(50)+acd86(82)
      acd86(82)=acd86(11)*acd86(82)
      acd86(63)=acd86(48)*acd86(63)
      acd86(63)=acd86(64)+acd86(68)+acd86(67)+acd86(65)+acd86(69)+acd86(63)+acd&
      &86(82)+acd86(81)+acd86(80)+acd86(79)+acd86(78)+acd86(77)+acd86(76)+acd86&
      &(75)+acd86(72)+acd86(73)+acd86(70)
      brack(ninjaidxt1mu0)=acd86(71)
      brack(ninjaidxt1mu2)=acd86(74)
      brack(ninjaidxt0mu0)=acd86(63)
      brack(ninjaidxt0mu2)=acd86(66)
   end subroutine brack_32
!---#] subroutine brack_32:
!---#[ subroutine numerator_t3:
   subroutine numerator_t3(ncut, a, b, c, param, deg, coeffs) &
   & bind(c, name="p10_part1part21_part25part25part1_d86h0_ninja_t3")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki => ki_nin
      use p10_part1part21_part25part25part1_globalsl1, only: epspow
      use p10_part1part21_part25part25part1_kinematics
      use p10_part1part21_part25part25part1_abbrevd86h0
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a, b, c
      complex(ki), intent(in) :: param
      complex(ki), dimension(4) :: vecA, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
	     vecA(1:4) = - a(0:3)
	     vecB(1:4) = - b(0:3)
	     vecC(1:4) = - c(0:3)
      if (deg.lt.0) return
      t1 = 0
      call cond_t(epspow.eq.t1,brack_31,vecA,vecB,vecC,param,coeffs)
      if (deg.le.(1+(0))) return
      call cond_t(epspow.eq.t1,brack_32,vecA,vecB,vecC,param,coeffs)
   end subroutine numerator_t3
!---#] subroutine numerator_t3:
end module     p10_part1part21_part25part25part1_d86h0l131
