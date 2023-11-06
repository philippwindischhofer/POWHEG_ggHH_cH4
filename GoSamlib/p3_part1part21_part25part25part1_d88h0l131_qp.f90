module     p3_part1part21_part25part25part1_d88h0l131_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p3_part1part2 &
   ! &1_part25part25part1/helicity0d88h0l131_qp.f90
   ! generator: buildfortran_tn3.py
   use p3_part1part21_part25part25part1_config, only: ki => ki_qp
   use p3_part1part21_part25part25part1_util_qp, only: cond_t, d => metric_tens&
   &or
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
      use p3_part1part21_part25part25part1_model_qp
      use p3_part1part21_part25part25part1_kinematics_qp
      use p3_part1part21_part25part25part1_color_qp
      use p3_part1part21_part25part25part1_abbrevd88h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(17) :: acd88
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd88(1)=dotproduct(ninjaA,ninjaE3)
      acd88(2)=dotproduct(ninjaE3,spvak5e2)
      acd88(3)=dotproduct(ninjaE3,spvae2k1)
      acd88(4)=abb88(9)
      acd88(5)=dotproduct(ninjaE3,spval4k1)
      acd88(6)=dotproduct(ninjaE3,spvae2l4)
      acd88(7)=abb88(7)
      acd88(8)=dotproduct(ninjaE3,spvak2k1)
      acd88(9)=dotproduct(ninjaE3,spvae2k2)
      acd88(10)=abb88(18)
      acd88(11)=dotproduct(ninjaE3,spvak5l4)
      acd88(12)=dotproduct(ninjaE3,spval4e2)
      acd88(13)=dotproduct(ninjaE3,spvak5k2)
      acd88(14)=dotproduct(ninjaE3,spvak2e2)
      acd88(15)=acd88(4)*acd88(1)*acd88(2)
      acd88(16)=acd88(12)*acd88(11)*acd88(7)
      acd88(17)=acd88(14)*acd88(13)*acd88(10)
      acd88(15)=acd88(17)+acd88(16)+2.0_ki*acd88(15)
      acd88(15)=acd88(3)*acd88(15)
      acd88(16)=acd88(6)*acd88(5)*acd88(7)
      acd88(17)=acd88(9)*acd88(8)*acd88(10)
      acd88(16)=acd88(17)+acd88(16)
      acd88(16)=acd88(2)*acd88(16)
      acd88(15)=acd88(15)+acd88(16)
      brack(ninjaidxt3mu0)=0.0_ki
      brack(ninjaidxt2mu0)=acd88(15)
   end subroutine brack_31
!---#] subroutine brack_31:
!---#[ subroutine brack_32:
   pure subroutine brack_32(ninjaA, ninjaE3, ninjaE4, ninjaP, brack)
      use p3_part1part21_part25part25part1_model_qp
      use p3_part1part21_part25part25part1_kinematics_qp
      use p3_part1part21_part25part25part1_color_qp
      use p3_part1part21_part25part25part1_abbrevd88h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(89) :: acd88
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd88(1)=dotproduct(ninjaA,ninjaE3)
      acd88(2)=dotproduct(ninjaE3,spvak5e2)
      acd88(3)=dotproduct(ninjaE4,spvae2k1)
      acd88(4)=abb88(9)
      acd88(5)=dotproduct(ninjaE3,spvae2k1)
      acd88(6)=dotproduct(ninjaE4,spvak5e2)
      acd88(7)=dotproduct(ninjaA,ninjaE4)
      acd88(8)=dotproduct(ninjaA,spvak5e2)
      acd88(9)=dotproduct(ninjaA,spvae2k1)
      acd88(10)=dotproduct(ninjaE3,spval4k1)
      acd88(11)=dotproduct(ninjaE3,spvae2l4)
      acd88(12)=abb88(7)
      acd88(13)=dotproduct(ninjaE4,spvae2l4)
      acd88(14)=dotproduct(ninjaE4,spval4k1)
      acd88(15)=dotproduct(ninjaE3,spvae2k2)
      acd88(16)=dotproduct(ninjaE4,spvak2k1)
      acd88(17)=abb88(18)
      acd88(18)=dotproduct(ninjaE3,spvak2k1)
      acd88(19)=dotproduct(ninjaE4,spvae2k2)
      acd88(20)=abb88(25)
      acd88(21)=dotproduct(ninjaE3,spvak5l4)
      acd88(22)=dotproduct(ninjaE4,spval4e2)
      acd88(23)=dotproduct(ninjaE3,spval4e2)
      acd88(24)=dotproduct(ninjaE4,spvak5l4)
      acd88(25)=dotproduct(ninjaE3,spvak2e2)
      acd88(26)=dotproduct(ninjaE4,spvak5k2)
      acd88(27)=dotproduct(ninjaE3,spvak5k2)
      acd88(28)=dotproduct(ninjaE4,spvak2e2)
      acd88(29)=abb88(23)
      acd88(30)=dotproduct(ninjaA,ninjaA)
      acd88(31)=dotproduct(ninjaA,spval4k1)
      acd88(32)=dotproduct(ninjaA,spvae2l4)
      acd88(33)=dotproduct(ninjaA,spvak5l4)
      acd88(34)=dotproduct(ninjaA,spval4e2)
      acd88(35)=dotproduct(ninjaA,spvae2k2)
      acd88(36)=dotproduct(ninjaA,spvak2e2)
      acd88(37)=dotproduct(ninjaA,spvak5k2)
      acd88(38)=dotproduct(ninjaA,spvak2k1)
      acd88(39)=abb88(28)
      acd88(40)=abb88(14)
      acd88(41)=abb88(30)
      acd88(42)=abb88(33)
      acd88(43)=abb88(21)
      acd88(44)=abb88(20)
      acd88(45)=abb88(16)
      acd88(46)=abb88(32)
      acd88(47)=abb88(22)
      acd88(48)=abb88(24)
      acd88(49)=abb88(27)
      acd88(50)=abb88(11)
      acd88(51)=abb88(31)
      acd88(52)=abb88(8)
      acd88(53)=abb88(26)
      acd88(54)=abb88(29)
      acd88(55)=abb88(12)
      acd88(56)=abb88(13)
      acd88(57)=abb88(17)
      acd88(58)=abb88(34)
      acd88(59)=dotproduct(ninjaE3,spval3e2)
      acd88(60)=abb88(19)
      acd88(61)=dotproduct(ninjaE3,spvae2l3)
      acd88(62)=abb88(36)
      acd88(63)=2.0_ki*acd88(1)
      acd88(64)=acd88(63)*acd88(3)
      acd88(64)=acd88(64)+acd88(9)
      acd88(64)=acd88(64)*acd88(4)
      acd88(65)=2.0_ki*acd88(4)
      acd88(65)=acd88(7)*acd88(65)*acd88(5)
      acd88(64)=acd88(65)+acd88(64)+acd88(20)
      acd88(65)=acd88(18)*acd88(19)
      acd88(66)=acd88(15)*acd88(16)
      acd88(65)=acd88(65)+acd88(66)
      acd88(66)=acd88(17)*acd88(65)
      acd88(67)=acd88(10)*acd88(13)
      acd88(68)=acd88(11)*acd88(14)
      acd88(67)=acd88(67)+acd88(68)
      acd88(68)=acd88(12)*acd88(67)
      acd88(66)=acd88(68)+acd88(66)+acd88(64)
      acd88(66)=acd88(2)*acd88(66)
      acd88(68)=acd88(27)*acd88(28)
      acd88(69)=acd88(25)*acd88(26)
      acd88(68)=acd88(68)+acd88(69)
      acd88(69)=acd88(17)*acd88(68)
      acd88(70)=acd88(21)*acd88(22)
      acd88(71)=acd88(23)*acd88(24)
      acd88(70)=acd88(70)+acd88(71)
      acd88(71)=acd88(12)*acd88(70)
      acd88(72)=acd88(63)*acd88(6)
      acd88(72)=acd88(72)+acd88(8)
      acd88(72)=acd88(72)*acd88(4)
      acd88(72)=acd88(72)+acd88(29)
      acd88(69)=acd88(71)+acd88(69)+acd88(72)
      acd88(69)=acd88(5)*acd88(69)
      acd88(71)=acd88(15)*acd88(18)
      acd88(73)=acd88(71)*acd88(6)
      acd88(74)=acd88(25)*acd88(27)
      acd88(75)=acd88(74)*acd88(3)
      acd88(73)=acd88(73)+acd88(75)
      acd88(75)=acd88(17)*acd88(73)
      acd88(76)=acd88(11)*acd88(10)
      acd88(77)=acd88(76)*acd88(6)
      acd88(78)=acd88(23)*acd88(21)
      acd88(79)=acd88(78)*acd88(3)
      acd88(77)=acd88(77)+acd88(79)
      acd88(79)=acd88(12)*acd88(77)
      acd88(66)=acd88(66)+acd88(69)+acd88(75)+acd88(79)
      acd88(69)=acd88(15)*acd88(42)
      acd88(75)=acd88(11)*acd88(40)
      acd88(79)=acd88(63)*acd88(20)
      acd88(69)=acd88(79)+acd88(69)+acd88(75)
      acd88(75)=acd88(4)*acd88(9)
      acd88(79)=acd88(63)*acd88(75)
      acd88(80)=acd88(18)*acd88(35)
      acd88(81)=acd88(15)*acd88(38)
      acd88(80)=acd88(80)+acd88(81)
      acd88(81)=acd88(17)*acd88(80)
      acd88(82)=acd88(10)*acd88(32)
      acd88(83)=acd88(11)*acd88(31)
      acd88(82)=acd88(82)+acd88(83)
      acd88(83)=acd88(12)*acd88(82)
      acd88(84)=ninjaP+acd88(30)
      acd88(84)=acd88(4)*acd88(84)
      acd88(84)=acd88(41)+acd88(84)
      acd88(84)=acd88(5)*acd88(84)
      acd88(79)=acd88(84)+acd88(83)+acd88(81)+acd88(79)+acd88(69)
      acd88(79)=acd88(2)*acd88(79)
      acd88(81)=acd88(4)*acd88(8)
      acd88(81)=acd88(81)+acd88(29)
      acd88(83)=acd88(63)*acd88(81)
      acd88(84)=acd88(25)*acd88(45)
      acd88(85)=acd88(23)*acd88(44)
      acd88(83)=acd88(84)+acd88(85)+acd88(83)
      acd88(84)=acd88(27)*acd88(36)
      acd88(85)=acd88(25)*acd88(37)
      acd88(84)=acd88(84)+acd88(85)
      acd88(85)=acd88(17)*acd88(84)
      acd88(86)=acd88(21)*acd88(34)
      acd88(87)=acd88(23)*acd88(33)
      acd88(86)=acd88(86)+acd88(87)
      acd88(87)=acd88(12)*acd88(86)
      acd88(85)=acd88(87)+acd88(85)+acd88(83)
      acd88(85)=acd88(5)*acd88(85)
      acd88(87)=acd88(9)*acd88(74)
      acd88(88)=acd88(8)*acd88(71)
      acd88(87)=acd88(87)+acd88(88)
      acd88(87)=acd88(17)*acd88(87)
      acd88(88)=acd88(9)*acd88(78)
      acd88(89)=acd88(8)*acd88(76)
      acd88(88)=acd88(88)+acd88(89)
      acd88(88)=acd88(12)*acd88(88)
      acd88(74)=acd88(47)*acd88(74)
      acd88(78)=acd88(43)*acd88(78)
      acd88(71)=acd88(46)*acd88(71)
      acd88(76)=acd88(39)*acd88(76)
      acd88(71)=acd88(79)+acd88(85)+acd88(88)+acd88(87)+acd88(76)+acd88(71)+acd&
      &88(74)+acd88(78)
      acd88(74)=acd88(2)*acd88(5)*acd88(4)
      acd88(68)=ninjaP*acd88(68)
      acd88(76)=acd88(36)*acd88(37)
      acd88(68)=acd88(76)+acd88(68)
      acd88(68)=acd88(17)*acd88(68)
      acd88(70)=ninjaP*acd88(70)
      acd88(76)=acd88(34)*acd88(33)
      acd88(70)=acd88(76)+acd88(70)
      acd88(70)=acd88(12)*acd88(70)
      acd88(76)=acd88(30)*acd88(81)
      acd88(78)=acd88(36)*acd88(45)
      acd88(79)=acd88(34)*acd88(44)
      acd88(81)=acd88(8)*acd88(41)
      acd88(72)=ninjaP*acd88(72)
      acd88(68)=acd88(70)+acd88(68)+acd88(72)+acd88(81)+acd88(79)+acd88(53)+acd&
      &88(78)+acd88(76)
      acd88(68)=acd88(5)*acd88(68)
      acd88(64)=ninjaP*acd88(64)
      acd88(65)=ninjaP*acd88(65)
      acd88(70)=acd88(35)*acd88(38)
      acd88(65)=acd88(70)+acd88(65)
      acd88(65)=acd88(17)*acd88(65)
      acd88(67)=ninjaP*acd88(67)
      acd88(70)=acd88(32)*acd88(31)
      acd88(67)=acd88(70)+acd88(67)
      acd88(67)=acd88(12)*acd88(67)
      acd88(70)=acd88(75)+acd88(20)
      acd88(70)=acd88(30)*acd88(70)
      acd88(72)=acd88(35)*acd88(42)
      acd88(75)=acd88(32)*acd88(40)
      acd88(76)=acd88(9)*acd88(41)
      acd88(64)=acd88(67)+acd88(65)+acd88(76)+acd88(75)+acd88(51)+acd88(72)+acd&
      &88(64)+acd88(70)
      acd88(64)=acd88(2)*acd88(64)
      acd88(65)=acd88(9)*acd88(84)
      acd88(67)=acd88(8)*acd88(80)
      acd88(70)=ninjaP*acd88(73)
      acd88(65)=acd88(70)+acd88(65)+acd88(67)
      acd88(65)=acd88(17)*acd88(65)
      acd88(67)=acd88(9)*acd88(86)
      acd88(70)=acd88(8)*acd88(82)
      acd88(72)=ninjaP*acd88(77)
      acd88(67)=acd88(72)+acd88(67)+acd88(70)
      acd88(67)=acd88(12)*acd88(67)
      acd88(70)=acd88(9)*acd88(83)
      acd88(69)=acd88(8)*acd88(69)
      acd88(72)=acd88(61)*acd88(62)
      acd88(73)=acd88(59)*acd88(60)
      acd88(75)=acd88(36)*acd88(47)
      acd88(75)=acd88(57)+acd88(75)
      acd88(75)=acd88(27)*acd88(75)
      acd88(76)=acd88(34)*acd88(43)
      acd88(76)=acd88(52)+acd88(76)
      acd88(76)=acd88(21)*acd88(76)
      acd88(77)=acd88(35)*acd88(46)
      acd88(77)=acd88(58)+acd88(77)
      acd88(77)=acd88(18)*acd88(77)
      acd88(78)=acd88(32)*acd88(39)
      acd88(78)=acd88(49)+acd88(78)
      acd88(78)=acd88(10)*acd88(78)
      acd88(79)=acd88(37)*acd88(47)
      acd88(79)=acd88(56)+acd88(79)
      acd88(79)=acd88(25)*acd88(79)
      acd88(80)=acd88(33)*acd88(43)
      acd88(80)=acd88(54)+acd88(80)
      acd88(80)=acd88(23)*acd88(80)
      acd88(81)=acd88(38)*acd88(46)
      acd88(81)=acd88(55)+acd88(81)
      acd88(81)=acd88(15)*acd88(81)
      acd88(82)=acd88(31)*acd88(39)
      acd88(82)=acd88(50)+acd88(82)
      acd88(82)=acd88(11)*acd88(82)
      acd88(63)=acd88(48)*acd88(63)
      acd88(63)=acd88(64)+acd88(68)+acd88(67)+acd88(65)+acd88(69)+acd88(63)+acd&
      &88(82)+acd88(81)+acd88(80)+acd88(79)+acd88(78)+acd88(77)+acd88(76)+acd88&
      &(75)+acd88(72)+acd88(73)+acd88(70)
      brack(ninjaidxt1mu0)=acd88(71)
      brack(ninjaidxt1mu2)=acd88(74)
      brack(ninjaidxt0mu0)=acd88(63)
      brack(ninjaidxt0mu2)=acd88(66)
   end subroutine brack_32
!---#] subroutine brack_32:
!---#[ subroutine numerator_t3:
   subroutine numerator_t3(ncut, a, b, c, param, deg, coeffs) &
   & bind(c, name="p3_part1part21_part25part25part1_d88h0_qp_ninja_t3")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p3_part1part21_part25part25part1_globalsl1_qp, only: epspow
      use p3_part1part21_part25part25part1_kinematics_qp
      use p3_part1part21_part25part25part1_abbrevd88h0_qp
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
end module     p3_part1part21_part25part25part1_d88h0l131_qp
