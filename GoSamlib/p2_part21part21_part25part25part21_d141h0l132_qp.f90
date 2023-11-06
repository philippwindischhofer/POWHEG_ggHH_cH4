module     p2_part21part21_part25part25part21_d141h0l132_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d141h0l132_qp.f90
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
      use p2_part21part21_part25part25part21_abbrevd141h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(22) :: acd141
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd141(1)=dotproduct(k2,ninjaE3)
      acd141(2)=dotproduct(ninjaE3,spvae1k2)
      acd141(3)=dotproduct(ninjaE3,spvak2e1)
      acd141(4)=abb141(9)
      acd141(5)=dotproduct(ninjaA0,ninjaE3)
      acd141(6)=abb141(11)
      acd141(7)=dotproduct(ninjaE3,spvak5e1)
      acd141(8)=dotproduct(ninjaE3,spvae1k5)
      acd141(9)=dotproduct(ninjaE3,spvak2k5)
      acd141(10)=abb141(7)
      acd141(11)=dotproduct(ninjaE3,spvak2l4)
      acd141(12)=dotproduct(ninjaE3,spval4e1)
      acd141(13)=abb141(8)
      acd141(14)=dotproduct(ninjaE3,spvae1l4)
      acd141(15)=dotproduct(ninjaE3,spval4k5)
      acd141(16)=dotproduct(ninjaE3,spvak5k2)
      acd141(17)=dotproduct(ninjaE3,spval4k2)
      acd141(18)=dotproduct(ninjaE3,spvak5l4)
      acd141(19)=acd141(15)*acd141(7)
      acd141(20)=acd141(17)*acd141(3)
      acd141(19)=acd141(20)+acd141(19)
      acd141(19)=acd141(19)*acd141(14)
      acd141(20)=acd141(11)*acd141(2)
      acd141(21)=acd141(18)*acd141(8)
      acd141(20)=acd141(21)+acd141(20)
      acd141(20)=acd141(20)*acd141(12)
      acd141(19)=acd141(20)+acd141(19)
      acd141(19)=acd141(13)*acd141(19)
      acd141(20)=-acd141(9)*acd141(7)*acd141(2)
      acd141(21)=-acd141(16)*acd141(8)*acd141(3)
      acd141(20)=acd141(21)+acd141(20)
      acd141(20)=acd141(10)*acd141(20)
      acd141(21)=acd141(3)*acd141(2)
      acd141(22)=acd141(8)*acd141(7)
      acd141(22)=acd141(21)+acd141(22)
      acd141(22)=acd141(6)*acd141(5)*acd141(22)
      acd141(21)=acd141(4)*acd141(1)*acd141(21)
      acd141(19)=2.0_ki*acd141(22)+acd141(21)+acd141(20)+acd141(19)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd141(19)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd141h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(93) :: acd141
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd141(1)=dotproduct(ninjaE3,spvae1k2)
      acd141(2)=dotproduct(ninjaE3,spvak2e1)
      acd141(3)=abb141(11)
      acd141(4)=dotproduct(ninjaE3,spvak5e1)
      acd141(5)=dotproduct(ninjaE3,spvae1k5)
      acd141(6)=dotproduct(k2,ninjaA1)
      acd141(7)=abb141(9)
      acd141(8)=dotproduct(k2,ninjaE3)
      acd141(9)=dotproduct(ninjaA1,spvae1k2)
      acd141(10)=dotproduct(ninjaA1,spvak2e1)
      acd141(11)=dotproduct(ninjaA0,ninjaA1)
      acd141(12)=dotproduct(ninjaA0,ninjaE3)
      acd141(13)=dotproduct(ninjaA1,spvak5e1)
      acd141(14)=dotproduct(ninjaA1,spvae1k5)
      acd141(15)=dotproduct(ninjaA1,spvak2k5)
      acd141(16)=abb141(7)
      acd141(17)=dotproduct(ninjaE3,spvak2k5)
      acd141(18)=dotproduct(ninjaE3,spvak2l4)
      acd141(19)=dotproduct(ninjaE3,spval4e1)
      acd141(20)=abb141(8)
      acd141(21)=dotproduct(ninjaE3,spvae1l4)
      acd141(22)=dotproduct(ninjaE3,spval4k5)
      acd141(23)=dotproduct(ninjaA1,spvak5k2)
      acd141(24)=dotproduct(ninjaE3,spvak5k2)
      acd141(25)=dotproduct(ninjaE3,spval4k2)
      acd141(26)=dotproduct(ninjaE3,spvak5l4)
      acd141(27)=dotproduct(ninjaA1,spvak2l4)
      acd141(28)=dotproduct(ninjaA1,spval4e1)
      acd141(29)=dotproduct(ninjaA1,spval4k2)
      acd141(30)=dotproduct(ninjaA1,spvae1l4)
      acd141(31)=dotproduct(ninjaA1,spval4k5)
      acd141(32)=dotproduct(ninjaA1,spvak5l4)
      acd141(33)=dotproduct(ninjaA1,ninjaA1)
      acd141(34)=dotproduct(k2,ninjaA0)
      acd141(35)=dotproduct(ninjaA0,spvae1k2)
      acd141(36)=dotproduct(ninjaA0,spvak2e1)
      acd141(37)=abb141(53)
      acd141(38)=abb141(30)
      acd141(39)=dotproduct(ninjaA0,ninjaA0)
      acd141(40)=dotproduct(ninjaA0,spvak5e1)
      acd141(41)=dotproduct(ninjaA0,spvae1k5)
      acd141(42)=abb141(26)
      acd141(43)=abb141(19)
      acd141(44)=abb141(29)
      acd141(45)=abb141(21)
      acd141(46)=dotproduct(ninjaA0,spvak2k5)
      acd141(47)=dotproduct(ninjaA0,spvak5k2)
      acd141(48)=dotproduct(ninjaA0,spvak2l4)
      acd141(49)=dotproduct(ninjaA0,spval4e1)
      acd141(50)=dotproduct(ninjaA0,spval4k2)
      acd141(51)=dotproduct(ninjaA0,spvae1l4)
      acd141(52)=dotproduct(ninjaA0,spval4k5)
      acd141(53)=dotproduct(ninjaA0,spvak5l4)
      acd141(54)=abb141(20)
      acd141(55)=abb141(17)
      acd141(56)=abb141(27)
      acd141(57)=abb141(37)
      acd141(58)=abb141(22)
      acd141(59)=abb141(35)
      acd141(60)=abb141(10)
      acd141(61)=abb141(14)
      acd141(62)=abb141(39)
      acd141(63)=abb141(32)
      acd141(64)=abb141(38)
      acd141(65)=abb141(44)
      acd141(66)=abb141(47)
      acd141(67)=abb141(34)
      acd141(68)=acd141(2)*acd141(1)
      acd141(69)=acd141(68)*acd141(3)
      acd141(70)=acd141(5)*acd141(3)
      acd141(71)=acd141(70)*acd141(4)
      acd141(69)=acd141(69)+acd141(71)
      acd141(71)=2.0_ki*acd141(12)
      acd141(72)=acd141(71)*acd141(3)
      acd141(73)=acd141(8)*acd141(7)
      acd141(73)=acd141(72)+acd141(73)
      acd141(74)=acd141(2)*acd141(73)
      acd141(75)=acd141(19)*acd141(20)
      acd141(76)=acd141(75)*acd141(18)
      acd141(77)=acd141(17)*acd141(16)
      acd141(78)=acd141(77)*acd141(4)
      acd141(74)=acd141(76)-acd141(78)+acd141(74)
      acd141(76)=acd141(9)*acd141(74)
      acd141(73)=-acd141(1)*acd141(73)
      acd141(78)=acd141(21)*acd141(20)
      acd141(79)=acd141(78)*acd141(25)
      acd141(80)=acd141(24)*acd141(16)
      acd141(81)=acd141(80)*acd141(5)
      acd141(73)=-acd141(79)+acd141(81)+acd141(73)
      acd141(79)=-acd141(10)*acd141(73)
      acd141(70)=acd141(70)*acd141(71)
      acd141(81)=acd141(78)*acd141(22)
      acd141(77)=acd141(77)*acd141(1)
      acd141(70)=-acd141(77)+acd141(70)+acd141(81)
      acd141(77)=acd141(13)*acd141(70)
      acd141(72)=acd141(4)*acd141(72)
      acd141(81)=acd141(75)*acd141(26)
      acd141(80)=acd141(80)*acd141(2)
      acd141(72)=acd141(72)+acd141(81)-acd141(80)
      acd141(80)=acd141(14)*acd141(72)
      acd141(81)=acd141(16)*acd141(4)
      acd141(82)=-acd141(15)*acd141(81)
      acd141(83)=acd141(27)*acd141(75)
      acd141(82)=acd141(83)+acd141(82)
      acd141(82)=acd141(1)*acd141(82)
      acd141(83)=acd141(29)*acd141(2)
      acd141(84)=acd141(31)*acd141(4)
      acd141(83)=acd141(84)+acd141(83)
      acd141(83)=acd141(78)*acd141(83)
      acd141(84)=ninjaP1+2.0_ki*acd141(11)
      acd141(84)=acd141(69)*acd141(84)
      acd141(85)=acd141(18)*acd141(1)
      acd141(86)=acd141(26)*acd141(5)
      acd141(85)=acd141(85)+acd141(86)
      acd141(85)=acd141(85)*acd141(20)
      acd141(86)=acd141(28)*acd141(85)
      acd141(87)=acd141(22)*acd141(4)
      acd141(88)=acd141(25)*acd141(2)
      acd141(87)=acd141(87)+acd141(88)
      acd141(87)=acd141(87)*acd141(20)
      acd141(88)=acd141(30)*acd141(87)
      acd141(89)=acd141(6)*acd141(7)*acd141(68)
      acd141(90)=acd141(16)*acd141(5)
      acd141(91)=-acd141(23)*acd141(2)*acd141(90)
      acd141(92)=acd141(75)*acd141(5)
      acd141(93)=acd141(32)*acd141(92)
      acd141(76)=acd141(93)+acd141(91)+acd141(89)+acd141(88)+acd141(86)+acd141(&
      &80)+acd141(77)+acd141(76)+acd141(79)+acd141(84)+acd141(83)+acd141(82)
      acd141(77)=acd141(33)+ninjaP2
      acd141(77)=acd141(69)*acd141(77)
      acd141(79)=acd141(38)*acd141(8)
      acd141(80)=acd141(44)*acd141(71)
      acd141(82)=-acd141(47)*acd141(90)
      acd141(83)=acd141(50)*acd141(78)
      acd141(84)=acd141(60)*acd141(24)
      acd141(86)=acd141(61)*acd141(5)
      acd141(88)=acd141(62)*acd141(21)
      acd141(79)=acd141(84)+acd141(88)+acd141(86)+acd141(83)+acd141(82)+acd141(&
      &80)+acd141(79)
      acd141(79)=acd141(2)*acd141(79)
      acd141(80)=acd141(37)*acd141(8)
      acd141(82)=acd141(42)*acd141(71)
      acd141(81)=-acd141(46)*acd141(81)
      acd141(75)=acd141(48)*acd141(75)
      acd141(83)=acd141(54)*acd141(17)
      acd141(84)=acd141(55)*acd141(4)
      acd141(86)=acd141(57)*acd141(19)
      acd141(75)=acd141(86)+acd141(84)+acd141(83)+acd141(75)+acd141(81)+acd141(&
      &82)+acd141(80)
      acd141(75)=acd141(1)*acd141(75)
      acd141(80)=acd141(63)*acd141(5)
      acd141(81)=acd141(64)*acd141(18)
      acd141(82)=acd141(65)*acd141(26)
      acd141(80)=acd141(82)+acd141(81)+acd141(80)
      acd141(80)=acd141(19)*acd141(80)
      acd141(81)=acd141(59)*acd141(4)
      acd141(82)=acd141(66)*acd141(25)
      acd141(83)=acd141(67)*acd141(22)
      acd141(81)=acd141(83)+acd141(82)+acd141(81)
      acd141(81)=acd141(21)*acd141(81)
      acd141(74)=acd141(35)*acd141(74)
      acd141(73)=-acd141(36)*acd141(73)
      acd141(70)=acd141(40)*acd141(70)
      acd141(72)=acd141(41)*acd141(72)
      acd141(78)=acd141(52)*acd141(78)
      acd141(82)=acd141(58)*acd141(5)
      acd141(78)=acd141(82)+acd141(78)
      acd141(78)=acd141(4)*acd141(78)
      acd141(82)=acd141(43)*acd141(4)
      acd141(83)=acd141(45)*acd141(5)
      acd141(82)=acd141(83)+acd141(82)
      acd141(71)=acd141(71)*acd141(82)
      acd141(82)=acd141(34)*acd141(7)
      acd141(82)=acd141(56)+acd141(82)
      acd141(68)=acd141(68)*acd141(82)
      acd141(82)=ninjaP0+acd141(39)
      acd141(82)=acd141(69)*acd141(82)
      acd141(83)=acd141(49)*acd141(85)
      acd141(84)=acd141(51)*acd141(87)
      acd141(85)=acd141(53)*acd141(92)
      acd141(68)=acd141(85)+acd141(84)+acd141(83)+acd141(72)+acd141(70)+acd141(&
      &74)+acd141(73)+acd141(75)+acd141(81)+acd141(80)+acd141(82)+acd141(68)+ac&
      &d141(71)+acd141(78)+acd141(79)
      brack(ninjaidxt0x0mu0)=acd141(68)
      brack(ninjaidxt0x0mu2)=acd141(69)
      brack(ninjaidxt0x1mu0)=acd141(76)
      brack(ninjaidxt0x2mu0)=acd141(77)
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d141h0_qp_ninja_t2")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd141h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a0, a1, b, c
      complex(ki), dimension(0:2), intent(in) :: param
      complex(ki), dimension(4) :: vecA0, vecA1, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      qshift = k3+k5+k4
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
end module     p2_part21part21_part25part25part21_d141h0l132_qp
