module     p2_part21part21_part25part25part21_d206h0l132_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d206h0l132_qp.f90
   ! generator: buildfortran_tn3.py
   use p2_part21part21_part25part25part21_config, only: ki => ki_qp
   use p2_part21part21_part25part25part21_util_qp, only: cond_t, d => metric_te&
   &nsor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, parameter :: ninjaidxt2x0mu0 = 0
   integer, parameter :: ninjaidxt1x0mu0 = 1
   integer, parameter :: ninjaidxt1x1mu0 = 2
   integer, parameter :: ninjaidxt0x0mu0 = 3
   integer, parameter :: ninjaidxt0x0mu2 = 4
   integer, parameter :: ninjaidxt0x1mu0 = 5
   integer, parameter :: ninjaidxt0x2mu0 = 6
   public :: numerator_t2
contains
!---#[ subroutine brack_21:
   pure subroutine brack_21(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd206h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(21) :: acd206
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd206(1)=dotproduct(ninjaA0,ninjaE3)
      acd206(2)=dotproduct(ninjaE3,spvae1e5)
      acd206(3)=dotproduct(ninjaE3,spvae5e2)
      acd206(4)=dotproduct(ninjaE3,spvae2e1)
      acd206(5)=abb206(81)
      acd206(6)=dotproduct(ninjaE3,spvae5e1)
      acd206(7)=dotproduct(ninjaE3,spvae2e5)
      acd206(8)=dotproduct(ninjaE3,spvae1e2)
      acd206(9)=dotproduct(ninjaE3,spvak2e1)
      acd206(10)=dotproduct(ninjaE3,spvae2k2)
      acd206(11)=abb206(66)
      acd206(12)=dotproduct(ninjaE3,spval4e1)
      acd206(13)=dotproduct(ninjaE3,spvae2l4)
      acd206(14)=dotproduct(ninjaE3,spvae1k2)
      acd206(15)=dotproduct(ninjaE3,spvak2e2)
      acd206(16)=dotproduct(ninjaE3,spvae1l4)
      acd206(17)=dotproduct(ninjaE3,spval4e2)
      acd206(18)=acd206(16)*acd206(17)
      acd206(19)=2.0_ki*acd206(1)
      acd206(20)=-acd206(8)*acd206(19)
      acd206(18)=acd206(18)+acd206(20)
      acd206(20)=acd206(6)*acd206(7)
      acd206(18)=acd206(18)*acd206(20)
      acd206(21)=acd206(12)*acd206(13)
      acd206(19)=-acd206(4)*acd206(19)
      acd206(19)=acd206(21)+acd206(19)
      acd206(21)=acd206(2)*acd206(3)
      acd206(19)=acd206(19)*acd206(21)
      acd206(18)=acd206(18)+acd206(19)
      acd206(18)=acd206(5)*acd206(18)
      acd206(19)=acd206(14)*acd206(15)*acd206(20)
      acd206(20)=acd206(9)*acd206(10)*acd206(21)
      acd206(19)=acd206(19)+acd206(20)
      acd206(19)=acd206(11)*acd206(19)
      acd206(18)=acd206(18)+acd206(19)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd206(18)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd206h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(89) :: acd206
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd206(1)=dotproduct(ninjaE3,spvae1e5)
      acd206(2)=dotproduct(ninjaE3,spvae5e2)
      acd206(3)=dotproduct(ninjaE3,spvae2e1)
      acd206(4)=abb206(81)
      acd206(5)=dotproduct(ninjaE3,spvae1e2)
      acd206(6)=dotproduct(ninjaE3,spvae5e1)
      acd206(7)=dotproduct(ninjaE3,spvae2e5)
      acd206(8)=dotproduct(ninjaA0,ninjaA1)
      acd206(9)=dotproduct(ninjaA0,ninjaE3)
      acd206(10)=dotproduct(ninjaA1,spvae1e5)
      acd206(11)=dotproduct(ninjaA1,spvae5e2)
      acd206(12)=dotproduct(ninjaA1,spvae1e2)
      acd206(13)=dotproduct(ninjaA1,spvae5e1)
      acd206(14)=dotproduct(ninjaA1,spvae2e5)
      acd206(15)=dotproduct(ninjaA1,spvae2e1)
      acd206(16)=dotproduct(ninjaA1,spval4e1)
      acd206(17)=dotproduct(ninjaE3,spvae2l4)
      acd206(18)=dotproduct(ninjaE3,spval4e1)
      acd206(19)=dotproduct(ninjaE3,spvak2e1)
      acd206(20)=dotproduct(ninjaE3,spvae2k2)
      acd206(21)=abb206(66)
      acd206(22)=dotproduct(ninjaE3,spvae1k2)
      acd206(23)=dotproduct(ninjaE3,spvak2e2)
      acd206(24)=dotproduct(ninjaE3,spvae1l4)
      acd206(25)=dotproduct(ninjaE3,spval4e2)
      acd206(26)=dotproduct(ninjaA1,spvak2e1)
      acd206(27)=dotproduct(ninjaA1,spvae1k2)
      acd206(28)=dotproduct(ninjaA1,spvak2e2)
      acd206(29)=dotproduct(ninjaA1,spvae2k2)
      acd206(30)=dotproduct(ninjaA1,spvae1l4)
      acd206(31)=dotproduct(ninjaA1,spvae2l4)
      acd206(32)=dotproduct(ninjaA1,spval4e2)
      acd206(33)=dotproduct(ninjaA1,ninjaA1)
      acd206(34)=dotproduct(ninjaA0,ninjaA0)
      acd206(35)=dotproduct(ninjaA0,spvae1e5)
      acd206(36)=dotproduct(ninjaA0,spvae5e2)
      acd206(37)=dotproduct(ninjaA0,spvae1e2)
      acd206(38)=dotproduct(ninjaA0,spvae5e1)
      acd206(39)=dotproduct(ninjaA0,spvae2e5)
      acd206(40)=dotproduct(ninjaA0,spvae2e1)
      acd206(41)=abb206(42)
      acd206(42)=abb206(76)
      acd206(43)=abb206(71)
      acd206(44)=abb206(84)
      acd206(45)=abb206(83)
      acd206(46)=abb206(52)
      acd206(47)=dotproduct(ninjaA0,spval4e1)
      acd206(48)=dotproduct(ninjaA0,spvak2e1)
      acd206(49)=dotproduct(ninjaA0,spvae1k2)
      acd206(50)=dotproduct(ninjaA0,spvak2e2)
      acd206(51)=dotproduct(ninjaA0,spvae2k2)
      acd206(52)=dotproduct(ninjaA0,spvae1l4)
      acd206(53)=dotproduct(ninjaA0,spvae2l4)
      acd206(54)=dotproduct(ninjaA0,spval4e2)
      acd206(55)=abb206(6)
      acd206(56)=abb206(7)
      acd206(57)=abb206(12)
      acd206(58)=abb206(29)
      acd206(59)=abb206(38)
      acd206(60)=abb206(36)
      acd206(61)=abb206(37)
      acd206(62)=abb206(17)
      acd206(63)=abb206(31)
      acd206(64)=acd206(2)*acd206(1)
      acd206(65)=acd206(64)*acd206(4)
      acd206(66)=acd206(65)*acd206(3)
      acd206(67)=acd206(7)*acd206(6)
      acd206(68)=acd206(67)*acd206(4)
      acd206(69)=acd206(68)*acd206(5)
      acd206(66)=acd206(66)+acd206(69)
      acd206(69)=2.0_ki*acd206(9)
      acd206(70)=-acd206(15)*acd206(69)
      acd206(71)=acd206(16)*acd206(17)
      acd206(72)=acd206(31)*acd206(18)
      acd206(70)=acd206(72)+acd206(71)+acd206(70)
      acd206(70)=acd206(65)*acd206(70)
      acd206(71)=-acd206(12)*acd206(69)
      acd206(72)=acd206(30)*acd206(25)
      acd206(73)=acd206(32)*acd206(24)
      acd206(71)=acd206(73)+acd206(72)+acd206(71)
      acd206(71)=acd206(68)*acd206(71)
      acd206(72)=acd206(23)*acd206(21)
      acd206(73)=acd206(27)*acd206(72)
      acd206(74)=acd206(28)*acd206(22)*acd206(21)
      acd206(73)=acd206(74)+acd206(73)
      acd206(73)=acd206(67)*acd206(73)
      acd206(74)=acd206(20)*acd206(21)
      acd206(75)=acd206(26)*acd206(74)
      acd206(76)=acd206(19)*acd206(21)
      acd206(77)=acd206(29)*acd206(76)
      acd206(75)=acd206(77)+acd206(75)
      acd206(75)=acd206(64)*acd206(75)
      acd206(77)=-ninjaP1-2.0_ki*acd206(8)
      acd206(77)=acd206(66)*acd206(77)
      acd206(78)=acd206(17)*acd206(18)
      acd206(79)=acd206(69)*acd206(3)
      acd206(78)=acd206(78)-acd206(79)
      acd206(79)=-acd206(78)*acd206(2)*acd206(4)
      acd206(80)=acd206(20)*acd206(19)
      acd206(81)=acd206(80)*acd206(21)
      acd206(82)=acd206(81)*acd206(2)
      acd206(79)=-acd206(82)+acd206(79)
      acd206(82)=-acd206(10)*acd206(79)
      acd206(83)=-acd206(78)*acd206(1)*acd206(4)
      acd206(81)=acd206(81)*acd206(1)
      acd206(81)=-acd206(81)+acd206(83)
      acd206(83)=-acd206(11)*acd206(81)
      acd206(84)=acd206(69)*acd206(5)
      acd206(85)=acd206(25)*acd206(24)
      acd206(84)=acd206(84)-acd206(85)
      acd206(85)=acd206(84)*acd206(7)*acd206(4)
      acd206(86)=acd206(23)*acd206(22)
      acd206(87)=acd206(86)*acd206(21)
      acd206(88)=acd206(87)*acd206(7)
      acd206(85)=-acd206(88)+acd206(85)
      acd206(88)=-acd206(13)*acd206(85)
      acd206(89)=acd206(84)*acd206(6)*acd206(4)
      acd206(87)=acd206(87)*acd206(6)
      acd206(87)=-acd206(87)+acd206(89)
      acd206(89)=-acd206(14)*acd206(87)
      acd206(70)=acd206(89)+acd206(88)+acd206(82)+acd206(83)+acd206(71)+acd206(&
      &70)+acd206(77)+acd206(75)+acd206(73)
      acd206(71)=-acd206(33)-ninjaP2
      acd206(71)=acd206(66)*acd206(71)
      acd206(73)=acd206(41)*acd206(69)
      acd206(74)=acd206(48)*acd206(74)
      acd206(75)=acd206(51)*acd206(76)
      acd206(76)=acd206(55)*acd206(18)
      acd206(77)=acd206(57)*acd206(19)
      acd206(73)=acd206(77)+acd206(76)+acd206(75)+acd206(74)+acd206(73)
      acd206(73)=acd206(64)*acd206(73)
      acd206(74)=acd206(50)*acd206(21)
      acd206(74)=acd206(60)+acd206(74)
      acd206(74)=acd206(22)*acd206(74)
      acd206(75)=acd206(46)*acd206(69)
      acd206(72)=acd206(49)*acd206(72)
      acd206(76)=acd206(61)*acd206(24)
      acd206(72)=acd206(76)+acd206(72)+acd206(75)+acd206(74)
      acd206(72)=acd206(67)*acd206(72)
      acd206(74)=acd206(42)*acd206(1)
      acd206(75)=acd206(43)*acd206(2)
      acd206(74)=acd206(75)+acd206(74)
      acd206(74)=acd206(78)*acd206(74)
      acd206(75)=-acd206(40)*acd206(69)
      acd206(76)=acd206(47)*acd206(17)
      acd206(77)=acd206(53)*acd206(18)
      acd206(75)=acd206(77)+acd206(76)+acd206(75)
      acd206(65)=acd206(65)*acd206(75)
      acd206(69)=-acd206(37)*acd206(69)
      acd206(75)=acd206(52)*acd206(25)
      acd206(76)=acd206(54)*acd206(24)
      acd206(69)=acd206(76)+acd206(75)+acd206(69)
      acd206(68)=acd206(68)*acd206(69)
      acd206(69)=acd206(62)*acd206(6)
      acd206(75)=acd206(63)*acd206(7)
      acd206(69)=acd206(75)+acd206(69)
      acd206(69)=acd206(86)*acd206(69)
      acd206(75)=acd206(58)*acd206(1)
      acd206(76)=acd206(59)*acd206(2)
      acd206(75)=acd206(76)+acd206(75)
      acd206(75)=acd206(80)*acd206(75)
      acd206(76)=-ninjaP0-acd206(34)
      acd206(76)=acd206(66)*acd206(76)
      acd206(77)=-acd206(35)*acd206(79)
      acd206(78)=-acd206(36)*acd206(81)
      acd206(79)=-acd206(38)*acd206(85)
      acd206(80)=-acd206(39)*acd206(87)
      acd206(81)=-acd206(44)*acd206(6)*acd206(84)
      acd206(82)=-acd206(45)*acd206(7)*acd206(84)
      acd206(67)=acd206(67)*acd206(5)
      acd206(64)=acd206(64)*acd206(3)
      acd206(64)=-acd206(64)-acd206(67)
      acd206(64)=acd206(56)*acd206(64)
      acd206(64)=acd206(64)+acd206(82)+acd206(81)+acd206(80)+acd206(79)+acd206(&
      &77)+acd206(78)+acd206(73)+acd206(72)+acd206(68)+acd206(65)+acd206(76)+ac&
      &d206(74)+acd206(75)+acd206(69)
      brack(ninjaidxt0x0mu0)=acd206(64)
      brack(ninjaidxt0x0mu2)=-acd206(66)
      brack(ninjaidxt0x1mu0)=acd206(70)
      brack(ninjaidxt0x2mu0)=acd206(71)
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d206h0_qp_ninja_t2")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd206h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a0, a1, b, c
      complex(ki), dimension(0:2), intent(in) :: param
      complex(ki), dimension(4) :: vecA0, vecA1, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      qshift = -k3-k4
	     vecA0(1:4) = - a0(0:3) - qshift(1:4)
	     vecA1(1:4) = - a1(0:3)
	     vecB(1:4) = - b(0:3)
	     vecC(1:4) = - c(0:3)
      if (deg.lt.0) return
      t1 = 0
      call cond_t(epspow.eq.t1,brack_21,vecA0,vecA1,vecB,vecC,param,coeffs)
      if (deg.le.(1+(0))) return
      call cond_t(epspow.eq.t1,brack_22,vecA0,vecA1,vecB,vecC,param,coeffs)
   end subroutine numerator_t2
!---#] subroutine numerator_t2:
end module     p2_part21part21_part25part25part21_d206h0l132_qp
