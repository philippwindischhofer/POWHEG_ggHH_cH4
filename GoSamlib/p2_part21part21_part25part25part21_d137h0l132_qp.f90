module     p2_part21part21_part25part25part21_d137h0l132_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d137h0l132_qp.f90
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
      use p2_part21part21_part25part25part21_abbrevd137h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(22) :: acd137
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd137(1)=dotproduct(k2,ninjaE3)
      acd137(2)=dotproduct(ninjaE3,spvae1k2)
      acd137(3)=dotproduct(ninjaE3,spvak2e1)
      acd137(4)=abb137(9)
      acd137(5)=dotproduct(ninjaA0,ninjaE3)
      acd137(6)=abb137(11)
      acd137(7)=dotproduct(ninjaE3,spvak5e1)
      acd137(8)=dotproduct(ninjaE3,spvae1k5)
      acd137(9)=dotproduct(ninjaE3,spvak2k5)
      acd137(10)=abb137(7)
      acd137(11)=dotproduct(ninjaE3,spvak2l3)
      acd137(12)=dotproduct(ninjaE3,spval3e1)
      acd137(13)=abb137(8)
      acd137(14)=dotproduct(ninjaE3,spvae1l3)
      acd137(15)=dotproduct(ninjaE3,spval3k5)
      acd137(16)=dotproduct(ninjaE3,spvak5k2)
      acd137(17)=dotproduct(ninjaE3,spval3k2)
      acd137(18)=dotproduct(ninjaE3,spvak5l3)
      acd137(19)=acd137(15)*acd137(7)
      acd137(20)=acd137(17)*acd137(3)
      acd137(19)=acd137(20)+acd137(19)
      acd137(19)=acd137(19)*acd137(14)
      acd137(20)=acd137(11)*acd137(2)
      acd137(21)=acd137(18)*acd137(8)
      acd137(20)=acd137(21)+acd137(20)
      acd137(20)=acd137(20)*acd137(12)
      acd137(19)=acd137(20)+acd137(19)
      acd137(19)=acd137(13)*acd137(19)
      acd137(20)=-acd137(9)*acd137(7)*acd137(2)
      acd137(21)=-acd137(16)*acd137(8)*acd137(3)
      acd137(20)=acd137(21)+acd137(20)
      acd137(20)=acd137(10)*acd137(20)
      acd137(21)=acd137(3)*acd137(2)
      acd137(22)=acd137(8)*acd137(7)
      acd137(22)=acd137(21)+acd137(22)
      acd137(22)=acd137(6)*acd137(5)*acd137(22)
      acd137(21)=acd137(4)*acd137(1)*acd137(21)
      acd137(19)=2.0_ki*acd137(22)+acd137(21)+acd137(20)+acd137(19)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd137(19)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd137h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(93) :: acd137
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd137(1)=dotproduct(ninjaE3,spvae1k2)
      acd137(2)=dotproduct(ninjaE3,spvak2e1)
      acd137(3)=abb137(11)
      acd137(4)=dotproduct(ninjaE3,spvak5e1)
      acd137(5)=dotproduct(ninjaE3,spvae1k5)
      acd137(6)=dotproduct(k2,ninjaA1)
      acd137(7)=abb137(9)
      acd137(8)=dotproduct(k2,ninjaE3)
      acd137(9)=dotproduct(ninjaA1,spvae1k2)
      acd137(10)=dotproduct(ninjaA1,spvak2e1)
      acd137(11)=dotproduct(ninjaA0,ninjaA1)
      acd137(12)=dotproduct(ninjaA0,ninjaE3)
      acd137(13)=dotproduct(ninjaA1,spvak5e1)
      acd137(14)=dotproduct(ninjaA1,spvae1k5)
      acd137(15)=dotproduct(ninjaA1,spvak2k5)
      acd137(16)=abb137(7)
      acd137(17)=dotproduct(ninjaE3,spvak2k5)
      acd137(18)=dotproduct(ninjaE3,spvak2l3)
      acd137(19)=dotproduct(ninjaE3,spval3e1)
      acd137(20)=abb137(8)
      acd137(21)=dotproduct(ninjaE3,spvae1l3)
      acd137(22)=dotproduct(ninjaE3,spval3k5)
      acd137(23)=dotproduct(ninjaA1,spvak5k2)
      acd137(24)=dotproduct(ninjaE3,spvak5k2)
      acd137(25)=dotproduct(ninjaE3,spval3k2)
      acd137(26)=dotproduct(ninjaE3,spvak5l3)
      acd137(27)=dotproduct(ninjaA1,spvak2l3)
      acd137(28)=dotproduct(ninjaA1,spval3e1)
      acd137(29)=dotproduct(ninjaA1,spval3k2)
      acd137(30)=dotproduct(ninjaA1,spvae1l3)
      acd137(31)=dotproduct(ninjaA1,spval3k5)
      acd137(32)=dotproduct(ninjaA1,spvak5l3)
      acd137(33)=dotproduct(ninjaA1,ninjaA1)
      acd137(34)=dotproduct(k2,ninjaA0)
      acd137(35)=dotproduct(ninjaA0,spvae1k2)
      acd137(36)=dotproduct(ninjaA0,spvak2e1)
      acd137(37)=abb137(53)
      acd137(38)=abb137(30)
      acd137(39)=dotproduct(ninjaA0,ninjaA0)
      acd137(40)=dotproduct(ninjaA0,spvak5e1)
      acd137(41)=dotproduct(ninjaA0,spvae1k5)
      acd137(42)=abb137(26)
      acd137(43)=abb137(19)
      acd137(44)=abb137(29)
      acd137(45)=abb137(21)
      acd137(46)=dotproduct(ninjaA0,spvak2k5)
      acd137(47)=dotproduct(ninjaA0,spvak5k2)
      acd137(48)=dotproduct(ninjaA0,spvak2l3)
      acd137(49)=dotproduct(ninjaA0,spval3e1)
      acd137(50)=dotproduct(ninjaA0,spval3k2)
      acd137(51)=dotproduct(ninjaA0,spvae1l3)
      acd137(52)=dotproduct(ninjaA0,spval3k5)
      acd137(53)=dotproduct(ninjaA0,spvak5l3)
      acd137(54)=abb137(20)
      acd137(55)=abb137(17)
      acd137(56)=abb137(27)
      acd137(57)=abb137(37)
      acd137(58)=abb137(22)
      acd137(59)=abb137(35)
      acd137(60)=abb137(10)
      acd137(61)=abb137(14)
      acd137(62)=abb137(39)
      acd137(63)=abb137(32)
      acd137(64)=abb137(38)
      acd137(65)=abb137(45)
      acd137(66)=abb137(47)
      acd137(67)=abb137(34)
      acd137(68)=acd137(2)*acd137(1)
      acd137(69)=acd137(68)*acd137(3)
      acd137(70)=acd137(5)*acd137(3)
      acd137(71)=acd137(70)*acd137(4)
      acd137(69)=acd137(69)+acd137(71)
      acd137(71)=2.0_ki*acd137(12)
      acd137(72)=acd137(71)*acd137(3)
      acd137(73)=acd137(8)*acd137(7)
      acd137(73)=acd137(72)+acd137(73)
      acd137(74)=acd137(2)*acd137(73)
      acd137(75)=acd137(19)*acd137(20)
      acd137(76)=acd137(75)*acd137(18)
      acd137(77)=acd137(17)*acd137(16)
      acd137(78)=acd137(77)*acd137(4)
      acd137(74)=acd137(76)-acd137(78)+acd137(74)
      acd137(76)=acd137(9)*acd137(74)
      acd137(73)=-acd137(1)*acd137(73)
      acd137(78)=acd137(21)*acd137(20)
      acd137(79)=acd137(78)*acd137(25)
      acd137(80)=acd137(24)*acd137(16)
      acd137(81)=acd137(80)*acd137(5)
      acd137(73)=-acd137(79)+acd137(81)+acd137(73)
      acd137(79)=-acd137(10)*acd137(73)
      acd137(70)=acd137(70)*acd137(71)
      acd137(81)=acd137(78)*acd137(22)
      acd137(77)=acd137(77)*acd137(1)
      acd137(70)=-acd137(77)+acd137(70)+acd137(81)
      acd137(77)=acd137(13)*acd137(70)
      acd137(72)=acd137(4)*acd137(72)
      acd137(81)=acd137(75)*acd137(26)
      acd137(80)=acd137(80)*acd137(2)
      acd137(72)=acd137(72)+acd137(81)-acd137(80)
      acd137(80)=acd137(14)*acd137(72)
      acd137(81)=acd137(16)*acd137(4)
      acd137(82)=-acd137(15)*acd137(81)
      acd137(83)=acd137(27)*acd137(75)
      acd137(82)=acd137(83)+acd137(82)
      acd137(82)=acd137(1)*acd137(82)
      acd137(83)=acd137(29)*acd137(2)
      acd137(84)=acd137(31)*acd137(4)
      acd137(83)=acd137(84)+acd137(83)
      acd137(83)=acd137(78)*acd137(83)
      acd137(84)=ninjaP1+2.0_ki*acd137(11)
      acd137(84)=acd137(69)*acd137(84)
      acd137(85)=acd137(18)*acd137(1)
      acd137(86)=acd137(26)*acd137(5)
      acd137(85)=acd137(85)+acd137(86)
      acd137(85)=acd137(85)*acd137(20)
      acd137(86)=acd137(28)*acd137(85)
      acd137(87)=acd137(22)*acd137(4)
      acd137(88)=acd137(25)*acd137(2)
      acd137(87)=acd137(87)+acd137(88)
      acd137(87)=acd137(87)*acd137(20)
      acd137(88)=acd137(30)*acd137(87)
      acd137(89)=acd137(6)*acd137(7)*acd137(68)
      acd137(90)=acd137(16)*acd137(5)
      acd137(91)=-acd137(23)*acd137(2)*acd137(90)
      acd137(92)=acd137(75)*acd137(5)
      acd137(93)=acd137(32)*acd137(92)
      acd137(76)=acd137(93)+acd137(91)+acd137(89)+acd137(88)+acd137(86)+acd137(&
      &80)+acd137(77)+acd137(76)+acd137(79)+acd137(84)+acd137(83)+acd137(82)
      acd137(77)=acd137(33)+ninjaP2
      acd137(77)=acd137(69)*acd137(77)
      acd137(79)=acd137(38)*acd137(8)
      acd137(80)=acd137(44)*acd137(71)
      acd137(82)=-acd137(47)*acd137(90)
      acd137(83)=acd137(50)*acd137(78)
      acd137(84)=acd137(60)*acd137(24)
      acd137(86)=acd137(61)*acd137(5)
      acd137(88)=acd137(62)*acd137(21)
      acd137(79)=acd137(84)+acd137(88)+acd137(86)+acd137(83)+acd137(82)+acd137(&
      &80)+acd137(79)
      acd137(79)=acd137(2)*acd137(79)
      acd137(80)=acd137(37)*acd137(8)
      acd137(82)=acd137(42)*acd137(71)
      acd137(81)=-acd137(46)*acd137(81)
      acd137(75)=acd137(48)*acd137(75)
      acd137(83)=acd137(54)*acd137(17)
      acd137(84)=acd137(55)*acd137(4)
      acd137(86)=acd137(57)*acd137(19)
      acd137(75)=acd137(86)+acd137(84)+acd137(83)+acd137(75)+acd137(81)+acd137(&
      &82)+acd137(80)
      acd137(75)=acd137(1)*acd137(75)
      acd137(80)=acd137(63)*acd137(5)
      acd137(81)=acd137(64)*acd137(18)
      acd137(82)=acd137(65)*acd137(26)
      acd137(80)=acd137(82)+acd137(81)+acd137(80)
      acd137(80)=acd137(19)*acd137(80)
      acd137(81)=acd137(59)*acd137(4)
      acd137(82)=acd137(66)*acd137(25)
      acd137(83)=acd137(67)*acd137(22)
      acd137(81)=acd137(83)+acd137(82)+acd137(81)
      acd137(81)=acd137(21)*acd137(81)
      acd137(74)=acd137(35)*acd137(74)
      acd137(73)=-acd137(36)*acd137(73)
      acd137(70)=acd137(40)*acd137(70)
      acd137(72)=acd137(41)*acd137(72)
      acd137(78)=acd137(52)*acd137(78)
      acd137(82)=acd137(58)*acd137(5)
      acd137(78)=acd137(82)+acd137(78)
      acd137(78)=acd137(4)*acd137(78)
      acd137(82)=acd137(43)*acd137(4)
      acd137(83)=acd137(45)*acd137(5)
      acd137(82)=acd137(83)+acd137(82)
      acd137(71)=acd137(71)*acd137(82)
      acd137(82)=acd137(34)*acd137(7)
      acd137(82)=acd137(56)+acd137(82)
      acd137(68)=acd137(68)*acd137(82)
      acd137(82)=ninjaP0+acd137(39)
      acd137(82)=acd137(69)*acd137(82)
      acd137(83)=acd137(49)*acd137(85)
      acd137(84)=acd137(51)*acd137(87)
      acd137(85)=acd137(53)*acd137(92)
      acd137(68)=acd137(85)+acd137(84)+acd137(83)+acd137(72)+acd137(70)+acd137(&
      &74)+acd137(73)+acd137(75)+acd137(81)+acd137(80)+acd137(82)+acd137(68)+ac&
      &d137(71)+acd137(78)+acd137(79)
      brack(ninjaidxt0x0mu0)=acd137(68)
      brack(ninjaidxt0x0mu2)=acd137(69)
      brack(ninjaidxt0x1mu0)=acd137(76)
      brack(ninjaidxt0x2mu0)=acd137(77)
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d137h0_qp_ninja_t2")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd137h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a0, a1, b, c
      complex(ki), dimension(0:2), intent(in) :: param
      complex(ki), dimension(4) :: vecA0, vecA1, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      qshift = k3+k4
	     vecA0(1:4) = + a0(0:3) - qshift(1:4)
	     vecA1(1:4) = + a1(0:3)
	     vecB(1:4) = + b(0:3)
	     vecC(1:4) = + c(0:3)
      if (deg.lt.0) return
      t1 = 0
      call cond_t(epspow.eq.t1,brack_21,vecA0,vecA1,vecB,vecC,param,coeffs)
      if (deg.le.(1+(0))) return
      call cond_t(epspow.eq.t1,brack_22,vecA0,vecA1,vecB,vecC,param,coeffs)
   end subroutine numerator_t2
!---#] subroutine numerator_t2:
end module     p2_part21part21_part25part25part21_d137h0l132_qp
