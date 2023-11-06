module     p2_part21part21_part25part25part21_d56h0l131_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d56h0l131_qp.f90
   ! generator: buildfortran_tn3.py
   use p2_part21part21_part25part25part21_config, only: ki => ki_qp
   use p2_part21part21_part25part25part21_util_qp, only: cond_t, d => metric_te&
   &nsor
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
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd56h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(9) :: acd56
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd56(1)=dotproduct(ninjaE3,spvae1e2)
      acd56(2)=dotproduct(ninjaE3,spvae5e1)
      acd56(3)=dotproduct(ninjaE3,spvae2e5)
      acd56(4)=abb56(7)
      acd56(5)=dotproduct(ninjaE3,spvae2e1)
      acd56(6)=dotproduct(ninjaE3,spvae1e5)
      acd56(7)=dotproduct(ninjaE3,spvae5e2)
      acd56(8)=acd56(3)*acd56(2)*acd56(1)
      acd56(9)=acd56(7)*acd56(6)*acd56(5)
      acd56(8)=acd56(8)+acd56(9)
      acd56(8)=acd56(4)*acd56(8)
      brack(ninjaidxt3mu0)=0.0_ki
      brack(ninjaidxt2mu0)=acd56(8)
   end subroutine brack_31
!---#] subroutine brack_31:
!---#[ subroutine brack_32:
   pure subroutine brack_32(ninjaA, ninjaE3, ninjaE4, ninjaP, brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd56h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(99) :: acd56
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd56(1)=dotproduct(ninjaE3,spvae1e5)
      acd56(2)=dotproduct(ninjaE3,spvae5e2)
      acd56(3)=dotproduct(ninjaE4,spvae2e1)
      acd56(4)=abb56(7)
      acd56(5)=dotproduct(ninjaE3,spvae2e1)
      acd56(6)=dotproduct(ninjaE4,spvae5e2)
      acd56(7)=abb56(25)
      acd56(8)=dotproduct(ninjaE3,spvae5e1)
      acd56(9)=dotproduct(ninjaE3,spvae2e5)
      acd56(10)=dotproduct(ninjaE4,spvae1e2)
      acd56(11)=dotproduct(ninjaE3,spvae1e2)
      acd56(12)=dotproduct(ninjaE4,spvae2e5)
      acd56(13)=abb56(30)
      acd56(14)=dotproduct(ninjaE4,spvae5e1)
      acd56(15)=abb56(27)
      acd56(16)=dotproduct(ninjaE4,spvae1e5)
      acd56(17)=abb56(9)
      acd56(18)=dotproduct(ninjaA,ninjaE3)
      acd56(19)=dotproduct(ninjaA,spvae1e5)
      acd56(20)=dotproduct(ninjaA,spvae5e1)
      acd56(21)=dotproduct(ninjaA,spvae2e5)
      acd56(22)=dotproduct(ninjaA,spvae5e2)
      acd56(23)=dotproduct(ninjaA,spvae1e2)
      acd56(24)=dotproduct(ninjaA,spvae2e1)
      acd56(25)=dotproduct(ninjaE3,spvak1e2)
      acd56(26)=abb56(5)
      acd56(27)=abb56(6)
      acd56(28)=abb56(18)
      acd56(29)=dotproduct(ninjaE3,spvak2e1)
      acd56(30)=abb56(35)
      acd56(31)=dotproduct(ninjaE3,spvak5e1)
      acd56(32)=abb56(32)
      acd56(33)=abb56(4)
      acd56(34)=abb56(15)
      acd56(35)=dotproduct(ninjaE3,spvae2k1)
      acd56(36)=abb56(29)
      acd56(37)=dotproduct(ninjaE3,spvae1k2)
      acd56(38)=abb56(34)
      acd56(39)=dotproduct(ninjaE3,spvae1k5)
      acd56(40)=abb56(31)
      acd56(41)=abb56(16)
      acd56(42)=abb56(19)
      acd56(43)=abb56(26)
      acd56(44)=abb56(14)
      acd56(45)=abb56(23)
      acd56(46)=abb56(28)
      acd56(47)=dotproduct(ninjaA,ninjaA)
      acd56(48)=abb56(13)
      acd56(49)=dotproduct(ninjaA,spvak1e2)
      acd56(50)=dotproduct(ninjaA,spvae2k1)
      acd56(51)=dotproduct(ninjaA,spvak2e1)
      acd56(52)=dotproduct(ninjaA,spvae1k2)
      acd56(53)=dotproduct(ninjaA,spvak5e1)
      acd56(54)=dotproduct(ninjaA,spvae1k5)
      acd56(55)=abb56(2)
      acd56(56)=abb56(21)
      acd56(57)=abb56(20)
      acd56(58)=abb56(12)
      acd56(59)=abb56(17)
      acd56(60)=abb56(11)
      acd56(61)=abb56(8)
      acd56(62)=abb56(10)
      acd56(63)=abb56(24)
      acd56(64)=abb56(33)
      acd56(65)=acd56(11)*acd56(4)
      acd56(66)=acd56(65)*acd56(12)
      acd56(67)=acd56(9)*acd56(4)
      acd56(68)=acd56(67)*acd56(10)
      acd56(66)=acd56(66)+acd56(68)
      acd56(66)=acd56(66)*acd56(8)
      acd56(68)=acd56(5)*acd56(4)
      acd56(69)=acd56(68)*acd56(6)
      acd56(70)=acd56(2)*acd56(4)
      acd56(71)=acd56(70)*acd56(3)
      acd56(69)=acd56(69)+acd56(71)
      acd56(69)=acd56(69)*acd56(1)
      acd56(71)=acd56(2)*acd56(68)*acd56(16)
      acd56(72)=acd56(9)*acd56(65)*acd56(14)
      acd56(66)=acd56(69)+acd56(66)+acd56(71)+acd56(72)
      acd56(69)=acd56(7)*acd56(1)
      acd56(71)=acd56(13)*acd56(8)
      acd56(72)=acd56(15)*acd56(9)
      acd56(73)=acd56(17)*acd56(2)
      acd56(69)=acd56(69)+acd56(71)+acd56(72)+acd56(73)
      acd56(71)=acd56(69)+acd56(66)
      acd56(72)=acd56(26)*acd56(25)
      acd56(73)=acd56(30)*acd56(29)
      acd56(74)=acd56(32)*acd56(31)
      acd56(72)=acd56(74)+acd56(72)+acd56(73)
      acd56(73)=acd56(22)*acd56(68)
      acd56(74)=2.0_ki*acd56(18)
      acd56(75)=acd56(7)*acd56(74)
      acd56(76)=acd56(24)*acd56(70)
      acd56(77)=acd56(27)*acd56(2)
      acd56(78)=acd56(28)*acd56(5)
      acd56(73)=acd56(78)+acd56(77)+acd56(76)+acd56(75)+acd56(73)+acd56(72)
      acd56(73)=acd56(1)*acd56(73)
      acd56(75)=acd56(36)*acd56(35)
      acd56(76)=acd56(38)*acd56(37)
      acd56(77)=acd56(40)*acd56(39)
      acd56(75)=acd56(77)+acd56(75)+acd56(76)
      acd56(76)=acd56(21)*acd56(65)
      acd56(77)=acd56(13)*acd56(74)
      acd56(78)=acd56(23)*acd56(67)
      acd56(79)=acd56(33)*acd56(9)
      acd56(80)=acd56(34)*acd56(11)
      acd56(76)=acd56(80)+acd56(79)+acd56(78)+acd56(77)+acd56(76)+acd56(75)
      acd56(76)=acd56(8)*acd56(76)
      acd56(77)=acd56(44)*acd56(35)
      acd56(78)=acd56(45)*acd56(37)
      acd56(79)=acd56(46)*acd56(39)
      acd56(68)=acd56(68)*acd56(19)
      acd56(68)=acd56(77)+acd56(78)+acd56(79)+acd56(68)
      acd56(77)=acd56(17)*acd56(74)
      acd56(77)=acd56(77)+acd56(68)
      acd56(77)=acd56(2)*acd56(77)
      acd56(78)=acd56(41)*acd56(25)
      acd56(79)=acd56(42)*acd56(29)
      acd56(80)=acd56(43)*acd56(31)
      acd56(65)=acd56(65)*acd56(20)
      acd56(65)=acd56(78)+acd56(79)+acd56(80)+acd56(65)
      acd56(78)=acd56(15)*acd56(74)
      acd56(78)=acd56(78)+acd56(65)
      acd56(78)=acd56(9)*acd56(78)
      acd56(73)=acd56(78)+acd56(77)+acd56(76)+acd56(73)
      acd56(65)=acd56(21)*acd56(65)
      acd56(68)=acd56(22)*acd56(68)
      acd56(66)=ninjaP*acd56(66)
      acd56(69)=acd56(47)*acd56(69)
      acd56(72)=acd56(19)*acd56(72)
      acd56(75)=acd56(20)*acd56(75)
      acd56(76)=ninjaP*acd56(1)
      acd56(77)=acd56(19)*acd56(74)
      acd56(76)=acd56(76)+acd56(77)
      acd56(76)=acd56(7)*acd56(76)
      acd56(77)=ninjaP*acd56(8)
      acd56(78)=acd56(20)*acd56(74)
      acd56(77)=acd56(77)+acd56(78)
      acd56(77)=acd56(13)*acd56(77)
      acd56(78)=ninjaP*acd56(9)
      acd56(79)=acd56(21)*acd56(74)
      acd56(78)=acd56(78)+acd56(79)
      acd56(78)=acd56(15)*acd56(78)
      acd56(79)=ninjaP*acd56(2)
      acd56(80)=acd56(22)*acd56(74)
      acd56(79)=acd56(79)+acd56(80)
      acd56(79)=acd56(17)*acd56(79)
      acd56(67)=acd56(20)*acd56(67)
      acd56(80)=acd56(21)*acd56(8)
      acd56(81)=acd56(4)*acd56(80)
      acd56(67)=acd56(67)+acd56(81)
      acd56(67)=acd56(23)*acd56(67)
      acd56(70)=acd56(19)*acd56(70)
      acd56(81)=acd56(22)*acd56(1)
      acd56(82)=acd56(4)*acd56(81)
      acd56(70)=acd56(70)+acd56(82)
      acd56(70)=acd56(24)*acd56(70)
      acd56(82)=acd56(19)*acd56(2)
      acd56(81)=acd56(82)+acd56(81)
      acd56(81)=acd56(27)*acd56(81)
      acd56(82)=acd56(19)*acd56(5)
      acd56(83)=acd56(24)*acd56(1)
      acd56(82)=acd56(82)+acd56(83)
      acd56(82)=acd56(28)*acd56(82)
      acd56(83)=acd56(20)*acd56(9)
      acd56(80)=acd56(83)+acd56(80)
      acd56(80)=acd56(33)*acd56(80)
      acd56(83)=acd56(20)*acd56(11)
      acd56(84)=acd56(23)*acd56(8)
      acd56(83)=acd56(83)+acd56(84)
      acd56(83)=acd56(34)*acd56(83)
      acd56(84)=acd56(26)*acd56(1)
      acd56(85)=acd56(41)*acd56(9)
      acd56(84)=acd56(84)+acd56(85)
      acd56(84)=acd56(49)*acd56(84)
      acd56(85)=acd56(36)*acd56(8)
      acd56(86)=acd56(44)*acd56(2)
      acd56(85)=acd56(85)+acd56(86)
      acd56(85)=acd56(50)*acd56(85)
      acd56(86)=acd56(30)*acd56(1)
      acd56(87)=acd56(42)*acd56(9)
      acd56(86)=acd56(86)+acd56(87)
      acd56(86)=acd56(51)*acd56(86)
      acd56(87)=acd56(38)*acd56(8)
      acd56(88)=acd56(45)*acd56(2)
      acd56(87)=acd56(87)+acd56(88)
      acd56(87)=acd56(52)*acd56(87)
      acd56(88)=acd56(32)*acd56(1)
      acd56(89)=acd56(43)*acd56(9)
      acd56(88)=acd56(88)+acd56(89)
      acd56(88)=acd56(53)*acd56(88)
      acd56(89)=acd56(40)*acd56(8)
      acd56(90)=acd56(46)*acd56(2)
      acd56(89)=acd56(89)+acd56(90)
      acd56(89)=acd56(54)*acd56(89)
      acd56(74)=acd56(48)*acd56(74)
      acd56(90)=acd56(55)*acd56(1)
      acd56(91)=acd56(56)*acd56(8)
      acd56(92)=acd56(57)*acd56(9)
      acd56(93)=acd56(58)*acd56(25)
      acd56(94)=acd56(59)*acd56(2)
      acd56(95)=acd56(60)*acd56(11)
      acd56(96)=acd56(61)*acd56(5)
      acd56(97)=acd56(62)*acd56(35)
      acd56(98)=acd56(63)*acd56(29)
      acd56(99)=acd56(64)*acd56(37)
      acd56(65)=acd56(99)+acd56(98)+acd56(97)+acd56(96)+acd56(95)+acd56(94)+acd&
      &56(93)+acd56(92)+acd56(91)+acd56(90)+acd56(74)+acd56(89)+acd56(88)+acd56&
      &(87)+acd56(86)+acd56(85)+acd56(84)+acd56(83)+acd56(80)+acd56(82)+acd56(8&
      &1)+acd56(69)+acd56(70)+acd56(67)+acd56(79)+acd56(78)+acd56(77)+acd56(76)&
      &+acd56(66)+acd56(68)+acd56(65)+acd56(75)+acd56(72)
      brack(ninjaidxt1mu0)=acd56(73)
      brack(ninjaidxt1mu2)=0.0_ki
      brack(ninjaidxt0mu0)=acd56(65)
      brack(ninjaidxt0mu2)=acd56(71)
   end subroutine brack_32
!---#] subroutine brack_32:
!---#[ subroutine numerator_t3:
   subroutine numerator_t3(ncut, a, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d56h0_qp_ninja_t3")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd56h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a, b, c
      complex(ki), intent(in) :: param
      complex(ki), dimension(4) :: vecA, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      qshift = -k2+k5
	     vecA(1:4) = - a(0:3) - qshift(1:4)
	     vecB(1:4) = - b(0:3)
	     vecC(1:4) = - c(0:3)
      if (deg.lt.0) return
      t1 = 0
      call cond_t(epspow.eq.t1,brack_31,vecA,vecB,vecC,param,coeffs)
      if (deg.le.(1+(0))) return
      call cond_t(epspow.eq.t1,brack_32,vecA,vecB,vecC,param,coeffs)
   end subroutine numerator_t3
!---#] subroutine numerator_t3:
end module     p2_part21part21_part25part25part21_d56h0l131_qp
