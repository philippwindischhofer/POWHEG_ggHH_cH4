module     p2_part21part21_part25part25part21_d192h0l132_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d192h0l132_qp.f90
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
      use p2_part21part21_part25part25part21_abbrevd192h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(21) :: acd192
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd192(1)=dotproduct(ninjaA0,ninjaE3)
      acd192(2)=dotproduct(ninjaE3,spvae1e2)
      acd192(3)=dotproduct(ninjaE3,spvae5e1)
      acd192(4)=dotproduct(ninjaE3,spvae2e5)
      acd192(5)=abb192(31)
      acd192(6)=dotproduct(ninjaE3,spvae2e1)
      acd192(7)=dotproduct(ninjaE3,spvae1e5)
      acd192(8)=dotproduct(ninjaE3,spvae5e2)
      acd192(9)=dotproduct(ninjaE3,spvae2l3)
      acd192(10)=dotproduct(ninjaE3,spval3e5)
      acd192(11)=dotproduct(ninjaE3,spvae2k2)
      acd192(12)=dotproduct(ninjaE3,spvak2e5)
      acd192(13)=abb192(86)
      acd192(14)=dotproduct(ninjaE3,spval3e2)
      acd192(15)=dotproduct(ninjaE3,spvae5l3)
      acd192(16)=dotproduct(ninjaE3,spvak2e2)
      acd192(17)=dotproduct(ninjaE3,spvae5k2)
      acd192(18)=-acd192(14)*acd192(15)
      acd192(19)=2.0_ki*acd192(1)
      acd192(20)=acd192(8)*acd192(19)
      acd192(18)=acd192(18)+acd192(20)
      acd192(20)=acd192(6)*acd192(7)
      acd192(18)=acd192(18)*acd192(20)
      acd192(21)=-acd192(9)*acd192(10)
      acd192(19)=acd192(4)*acd192(19)
      acd192(19)=acd192(21)+acd192(19)
      acd192(21)=acd192(2)*acd192(3)
      acd192(19)=acd192(19)*acd192(21)
      acd192(18)=acd192(18)+acd192(19)
      acd192(18)=acd192(5)*acd192(18)
      acd192(19)=-acd192(16)*acd192(17)*acd192(20)
      acd192(20)=-acd192(11)*acd192(12)*acd192(21)
      acd192(19)=acd192(19)+acd192(20)
      acd192(19)=acd192(13)*acd192(19)
      acd192(18)=acd192(18)+acd192(19)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd192(18)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd192h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(92) :: acd192
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd192(1)=dotproduct(ninjaE3,spvae2e1)
      acd192(2)=dotproduct(ninjaE3,spvae1e5)
      acd192(3)=dotproduct(ninjaE3,spvae5e2)
      acd192(4)=abb192(31)
      acd192(5)=dotproduct(ninjaE3,spvae1e2)
      acd192(6)=dotproduct(ninjaE3,spvae5e1)
      acd192(7)=dotproduct(ninjaE3,spvae2e5)
      acd192(8)=dotproduct(ninjaA0,ninjaA1)
      acd192(9)=dotproduct(ninjaA0,ninjaE3)
      acd192(10)=dotproduct(ninjaA1,spvae2e1)
      acd192(11)=dotproduct(ninjaA1,spvae1e5)
      acd192(12)=dotproduct(ninjaA1,spvae1e2)
      acd192(13)=dotproduct(ninjaA1,spvae5e1)
      acd192(14)=dotproduct(ninjaA1,spvae2e5)
      acd192(15)=dotproduct(ninjaA1,spvae5e2)
      acd192(16)=dotproduct(ninjaA1,spvak2e2)
      acd192(17)=dotproduct(ninjaE3,spvae5k2)
      acd192(18)=abb192(86)
      acd192(19)=dotproduct(ninjaE3,spvak2e2)
      acd192(20)=dotproduct(ninjaE3,spval3e2)
      acd192(21)=dotproduct(ninjaE3,spvae5l3)
      acd192(22)=dotproduct(ninjaA1,spvae5k2)
      acd192(23)=dotproduct(ninjaA1,spval3e2)
      acd192(24)=dotproduct(ninjaA1,spvae2k2)
      acd192(25)=dotproduct(ninjaE3,spvak2e5)
      acd192(26)=dotproduct(ninjaE3,spvae2k2)
      acd192(27)=dotproduct(ninjaE3,spvae2l3)
      acd192(28)=dotproduct(ninjaE3,spval3e5)
      acd192(29)=dotproduct(ninjaA1,spvae2l3)
      acd192(30)=dotproduct(ninjaA1,spval3e5)
      acd192(31)=dotproduct(ninjaA1,spvae5l3)
      acd192(32)=dotproduct(ninjaA1,spvak2e5)
      acd192(33)=dotproduct(ninjaA1,ninjaA1)
      acd192(34)=dotproduct(ninjaA0,ninjaA0)
      acd192(35)=dotproduct(ninjaA0,spvae2e1)
      acd192(36)=dotproduct(ninjaA0,spvae1e5)
      acd192(37)=dotproduct(ninjaA0,spvae1e2)
      acd192(38)=dotproduct(ninjaA0,spvae5e1)
      acd192(39)=dotproduct(ninjaA0,spvae2e5)
      acd192(40)=dotproduct(ninjaA0,spvae5e2)
      acd192(41)=abb192(38)
      acd192(42)=abb192(63)
      acd192(43)=abb192(27)
      acd192(44)=abb192(62)
      acd192(45)=dotproduct(ninjaA0,spvak2e2)
      acd192(46)=dotproduct(ninjaA0,spvae5k2)
      acd192(47)=dotproduct(ninjaA0,spval3e2)
      acd192(48)=dotproduct(ninjaA0,spvae2k2)
      acd192(49)=dotproduct(ninjaA0,spvae2l3)
      acd192(50)=dotproduct(ninjaA0,spval3e5)
      acd192(51)=dotproduct(ninjaA0,spvae5l3)
      acd192(52)=dotproduct(ninjaA0,spvak2e5)
      acd192(53)=abb192(7)
      acd192(54)=abb192(8)
      acd192(55)=abb192(13)
      acd192(56)=abb192(30)
      acd192(57)=abb192(59)
      acd192(58)=dotproduct(ninjaE3,spvak2k1)
      acd192(59)=abb192(89)
      acd192(60)=dotproduct(ninjaE3,spvak2k5)
      acd192(61)=abb192(91)
      acd192(62)=abb192(23)
      acd192(63)=abb192(87)
      acd192(64)=abb192(25)
      acd192(65)=abb192(26)
      acd192(66)=dotproduct(ninjaE3,spvak1k2)
      acd192(67)=abb192(34)
      acd192(68)=dotproduct(ninjaE3,spvak5k2)
      acd192(69)=abb192(78)
      acd192(70)=acd192(2)*acd192(1)
      acd192(71)=acd192(70)*acd192(4)
      acd192(72)=acd192(71)*acd192(3)
      acd192(73)=acd192(6)*acd192(5)
      acd192(74)=acd192(73)*acd192(4)
      acd192(75)=acd192(74)*acd192(7)
      acd192(72)=acd192(72)+acd192(75)
      acd192(75)=2.0_ki*acd192(9)
      acd192(76)=acd192(15)*acd192(75)
      acd192(77)=-acd192(23)*acd192(21)
      acd192(78)=-acd192(31)*acd192(20)
      acd192(76)=acd192(78)+acd192(77)+acd192(76)
      acd192(71)=acd192(71)*acd192(76)
      acd192(76)=acd192(14)*acd192(75)
      acd192(77)=-acd192(29)*acd192(28)
      acd192(78)=-acd192(30)*acd192(27)
      acd192(76)=acd192(78)+acd192(77)+acd192(76)
      acd192(74)=acd192(74)*acd192(76)
      acd192(76)=acd192(25)*acd192(18)
      acd192(77)=-acd192(24)*acd192(76)
      acd192(78)=acd192(26)*acd192(18)
      acd192(79)=-acd192(32)*acd192(78)
      acd192(77)=acd192(79)+acd192(77)
      acd192(77)=acd192(73)*acd192(77)
      acd192(79)=acd192(17)*acd192(18)
      acd192(80)=-acd192(16)*acd192(79)
      acd192(81)=acd192(19)*acd192(18)
      acd192(82)=-acd192(22)*acd192(81)
      acd192(80)=acd192(82)+acd192(80)
      acd192(80)=acd192(70)*acd192(80)
      acd192(82)=ninjaP1+2.0_ki*acd192(8)
      acd192(82)=acd192(72)*acd192(82)
      acd192(83)=acd192(21)*acd192(20)
      acd192(84)=acd192(75)*acd192(3)
      acd192(83)=acd192(83)-acd192(84)
      acd192(84)=acd192(83)*acd192(2)*acd192(4)
      acd192(85)=acd192(79)*acd192(19)
      acd192(86)=acd192(85)*acd192(2)
      acd192(84)=acd192(86)+acd192(84)
      acd192(86)=-acd192(10)*acd192(84)
      acd192(87)=acd192(83)*acd192(1)*acd192(4)
      acd192(85)=acd192(85)*acd192(1)
      acd192(85)=acd192(85)+acd192(87)
      acd192(87)=-acd192(11)*acd192(85)
      acd192(88)=acd192(75)*acd192(7)
      acd192(89)=acd192(28)*acd192(27)
      acd192(88)=acd192(88)-acd192(89)
      acd192(89)=-acd192(88)*acd192(6)*acd192(4)
      acd192(90)=acd192(76)*acd192(26)
      acd192(91)=acd192(90)*acd192(6)
      acd192(89)=acd192(91)+acd192(89)
      acd192(91)=-acd192(12)*acd192(89)
      acd192(92)=-acd192(88)*acd192(5)*acd192(4)
      acd192(90)=acd192(90)*acd192(5)
      acd192(90)=acd192(90)+acd192(92)
      acd192(92)=-acd192(13)*acd192(90)
      acd192(71)=acd192(92)+acd192(91)+acd192(86)+acd192(87)+acd192(74)+acd192(&
      &71)+acd192(82)+acd192(80)+acd192(77)
      acd192(74)=acd192(33)+ninjaP2
      acd192(74)=acd192(72)*acd192(74)
      acd192(77)=acd192(39)*acd192(75)
      acd192(80)=-acd192(49)*acd192(28)
      acd192(82)=-acd192(50)*acd192(27)
      acd192(77)=acd192(82)+acd192(80)+acd192(77)
      acd192(77)=acd192(4)*acd192(77)
      acd192(80)=acd192(43)*acd192(75)
      acd192(76)=-acd192(48)*acd192(76)
      acd192(78)=-acd192(52)*acd192(78)
      acd192(82)=acd192(62)*acd192(26)
      acd192(86)=acd192(64)*acd192(27)
      acd192(87)=acd192(65)*acd192(28)
      acd192(91)=acd192(67)*acd192(66)
      acd192(92)=acd192(69)*acd192(68)
      acd192(76)=acd192(92)+acd192(91)+acd192(87)+acd192(86)+acd192(82)+acd192(&
      &78)+acd192(76)+acd192(80)+acd192(77)
      acd192(76)=acd192(73)*acd192(76)
      acd192(77)=acd192(40)*acd192(75)
      acd192(78)=-acd192(47)*acd192(21)
      acd192(80)=-acd192(51)*acd192(20)
      acd192(77)=acd192(80)+acd192(78)+acd192(77)
      acd192(77)=acd192(4)*acd192(77)
      acd192(75)=acd192(41)*acd192(75)
      acd192(78)=-acd192(45)*acd192(79)
      acd192(79)=-acd192(46)*acd192(81)
      acd192(80)=acd192(53)*acd192(19)
      acd192(81)=acd192(55)*acd192(20)
      acd192(82)=acd192(56)*acd192(21)
      acd192(86)=acd192(59)*acd192(58)
      acd192(87)=acd192(61)*acd192(60)
      acd192(75)=acd192(87)+acd192(86)+acd192(82)+acd192(81)+acd192(80)+acd192(&
      &79)+acd192(78)+acd192(75)+acd192(77)
      acd192(75)=acd192(70)*acd192(75)
      acd192(77)=acd192(42)*acd192(83)
      acd192(78)=acd192(54)*acd192(17)*acd192(19)
      acd192(77)=acd192(78)+acd192(77)
      acd192(77)=acd192(1)*acd192(77)
      acd192(78)=ninjaP0+acd192(34)
      acd192(78)=acd192(72)*acd192(78)
      acd192(79)=-acd192(35)*acd192(84)
      acd192(80)=-acd192(36)*acd192(85)
      acd192(81)=-acd192(37)*acd192(89)
      acd192(82)=-acd192(38)*acd192(90)
      acd192(83)=-acd192(44)*acd192(5)*acd192(88)
      acd192(73)=acd192(73)*acd192(7)
      acd192(70)=acd192(70)*acd192(3)
      acd192(70)=acd192(70)+acd192(73)
      acd192(70)=acd192(57)*acd192(70)
      acd192(73)=acd192(63)*acd192(25)*acd192(26)*acd192(5)
      acd192(70)=acd192(73)+acd192(70)+acd192(83)+acd192(82)+acd192(81)+acd192(&
      &79)+acd192(80)+acd192(75)+acd192(76)+acd192(78)+acd192(77)
      brack(ninjaidxt0x0mu0)=acd192(70)
      brack(ninjaidxt0x0mu2)=acd192(72)
      brack(ninjaidxt0x1mu0)=acd192(71)
      brack(ninjaidxt0x2mu0)=acd192(74)
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d192h0_qp_ninja_t2")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd192h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a0, a1, b, c
      complex(ki), dimension(0:2), intent(in) :: param
      complex(ki), dimension(4) :: vecA0, vecA1, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      qshift = k2
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
end module     p2_part21part21_part25part25part21_d192h0l132_qp
