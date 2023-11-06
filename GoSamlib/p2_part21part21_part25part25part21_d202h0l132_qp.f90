module     p2_part21part21_part25part25part21_d202h0l132_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d202h0l132_qp.f90
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
      use p2_part21part21_part25part25part21_abbrevd202h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(21) :: acd202
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd202(1)=dotproduct(ninjaA0,ninjaE3)
      acd202(2)=dotproduct(ninjaE3,spvae1e5)
      acd202(3)=dotproduct(ninjaE3,spvae5e2)
      acd202(4)=dotproduct(ninjaE3,spvae2e1)
      acd202(5)=abb202(81)
      acd202(6)=dotproduct(ninjaE3,spvae5e1)
      acd202(7)=dotproduct(ninjaE3,spvae2e5)
      acd202(8)=dotproduct(ninjaE3,spvae1e2)
      acd202(9)=dotproduct(ninjaE3,spvak2e1)
      acd202(10)=dotproduct(ninjaE3,spvae2k2)
      acd202(11)=abb202(63)
      acd202(12)=dotproduct(ninjaE3,spval3e1)
      acd202(13)=dotproduct(ninjaE3,spvae2l3)
      acd202(14)=dotproduct(ninjaE3,spvae1k2)
      acd202(15)=dotproduct(ninjaE3,spvak2e2)
      acd202(16)=dotproduct(ninjaE3,spvae1l3)
      acd202(17)=dotproduct(ninjaE3,spval3e2)
      acd202(18)=acd202(16)*acd202(17)
      acd202(19)=2.0_ki*acd202(1)
      acd202(20)=-acd202(8)*acd202(19)
      acd202(18)=acd202(18)+acd202(20)
      acd202(20)=acd202(6)*acd202(7)
      acd202(18)=acd202(18)*acd202(20)
      acd202(21)=acd202(12)*acd202(13)
      acd202(19)=-acd202(4)*acd202(19)
      acd202(19)=acd202(21)+acd202(19)
      acd202(21)=acd202(2)*acd202(3)
      acd202(19)=acd202(19)*acd202(21)
      acd202(18)=acd202(18)+acd202(19)
      acd202(18)=acd202(5)*acd202(18)
      acd202(19)=acd202(14)*acd202(15)*acd202(20)
      acd202(20)=acd202(9)*acd202(10)*acd202(21)
      acd202(19)=acd202(19)+acd202(20)
      acd202(19)=acd202(11)*acd202(19)
      acd202(18)=acd202(18)+acd202(19)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd202(18)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd202h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(89) :: acd202
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd202(1)=dotproduct(ninjaE3,spvae1e5)
      acd202(2)=dotproduct(ninjaE3,spvae5e2)
      acd202(3)=dotproduct(ninjaE3,spvae2e1)
      acd202(4)=abb202(81)
      acd202(5)=dotproduct(ninjaE3,spvae1e2)
      acd202(6)=dotproduct(ninjaE3,spvae5e1)
      acd202(7)=dotproduct(ninjaE3,spvae2e5)
      acd202(8)=dotproduct(ninjaA0,ninjaA1)
      acd202(9)=dotproduct(ninjaA0,ninjaE3)
      acd202(10)=dotproduct(ninjaA1,spvae1e5)
      acd202(11)=dotproduct(ninjaA1,spvae5e2)
      acd202(12)=dotproduct(ninjaA1,spvae1e2)
      acd202(13)=dotproduct(ninjaA1,spvae5e1)
      acd202(14)=dotproduct(ninjaA1,spvae2e5)
      acd202(15)=dotproduct(ninjaA1,spvae2e1)
      acd202(16)=dotproduct(ninjaA1,spval3e1)
      acd202(17)=dotproduct(ninjaE3,spvae2l3)
      acd202(18)=dotproduct(ninjaE3,spval3e1)
      acd202(19)=dotproduct(ninjaE3,spvak2e1)
      acd202(20)=dotproduct(ninjaE3,spvae2k2)
      acd202(21)=abb202(63)
      acd202(22)=dotproduct(ninjaE3,spvae1k2)
      acd202(23)=dotproduct(ninjaE3,spvak2e2)
      acd202(24)=dotproduct(ninjaE3,spvae1l3)
      acd202(25)=dotproduct(ninjaE3,spval3e2)
      acd202(26)=dotproduct(ninjaA1,spvak2e1)
      acd202(27)=dotproduct(ninjaA1,spvae1k2)
      acd202(28)=dotproduct(ninjaA1,spvak2e2)
      acd202(29)=dotproduct(ninjaA1,spvae2k2)
      acd202(30)=dotproduct(ninjaA1,spvae1l3)
      acd202(31)=dotproduct(ninjaA1,spvae2l3)
      acd202(32)=dotproduct(ninjaA1,spval3e2)
      acd202(33)=dotproduct(ninjaA1,ninjaA1)
      acd202(34)=dotproduct(ninjaA0,ninjaA0)
      acd202(35)=dotproduct(ninjaA0,spvae1e5)
      acd202(36)=dotproduct(ninjaA0,spvae5e2)
      acd202(37)=dotproduct(ninjaA0,spvae1e2)
      acd202(38)=dotproduct(ninjaA0,spvae5e1)
      acd202(39)=dotproduct(ninjaA0,spvae2e5)
      acd202(40)=dotproduct(ninjaA0,spvae2e1)
      acd202(41)=abb202(42)
      acd202(42)=abb202(76)
      acd202(43)=abb202(71)
      acd202(44)=abb202(84)
      acd202(45)=abb202(83)
      acd202(46)=abb202(52)
      acd202(47)=dotproduct(ninjaA0,spval3e1)
      acd202(48)=dotproduct(ninjaA0,spvak2e1)
      acd202(49)=dotproduct(ninjaA0,spvae1k2)
      acd202(50)=dotproduct(ninjaA0,spvak2e2)
      acd202(51)=dotproduct(ninjaA0,spvae2k2)
      acd202(52)=dotproduct(ninjaA0,spvae1l3)
      acd202(53)=dotproduct(ninjaA0,spvae2l3)
      acd202(54)=dotproduct(ninjaA0,spval3e2)
      acd202(55)=abb202(6)
      acd202(56)=abb202(7)
      acd202(57)=abb202(12)
      acd202(58)=abb202(29)
      acd202(59)=abb202(38)
      acd202(60)=abb202(36)
      acd202(61)=abb202(37)
      acd202(62)=abb202(17)
      acd202(63)=abb202(31)
      acd202(64)=acd202(2)*acd202(1)
      acd202(65)=acd202(64)*acd202(4)
      acd202(66)=acd202(65)*acd202(3)
      acd202(67)=acd202(7)*acd202(6)
      acd202(68)=acd202(67)*acd202(4)
      acd202(69)=acd202(68)*acd202(5)
      acd202(66)=acd202(66)+acd202(69)
      acd202(69)=2.0_ki*acd202(9)
      acd202(70)=-acd202(15)*acd202(69)
      acd202(71)=acd202(16)*acd202(17)
      acd202(72)=acd202(31)*acd202(18)
      acd202(70)=acd202(72)+acd202(71)+acd202(70)
      acd202(70)=acd202(65)*acd202(70)
      acd202(71)=-acd202(12)*acd202(69)
      acd202(72)=acd202(30)*acd202(25)
      acd202(73)=acd202(32)*acd202(24)
      acd202(71)=acd202(73)+acd202(72)+acd202(71)
      acd202(71)=acd202(68)*acd202(71)
      acd202(72)=acd202(23)*acd202(21)
      acd202(73)=acd202(27)*acd202(72)
      acd202(74)=acd202(28)*acd202(22)*acd202(21)
      acd202(73)=acd202(74)+acd202(73)
      acd202(73)=acd202(67)*acd202(73)
      acd202(74)=acd202(20)*acd202(21)
      acd202(75)=acd202(26)*acd202(74)
      acd202(76)=acd202(19)*acd202(21)
      acd202(77)=acd202(29)*acd202(76)
      acd202(75)=acd202(77)+acd202(75)
      acd202(75)=acd202(64)*acd202(75)
      acd202(77)=-ninjaP1-2.0_ki*acd202(8)
      acd202(77)=acd202(66)*acd202(77)
      acd202(78)=acd202(17)*acd202(18)
      acd202(79)=acd202(69)*acd202(3)
      acd202(78)=acd202(78)-acd202(79)
      acd202(79)=-acd202(78)*acd202(2)*acd202(4)
      acd202(80)=acd202(20)*acd202(19)
      acd202(81)=acd202(80)*acd202(21)
      acd202(82)=acd202(81)*acd202(2)
      acd202(79)=-acd202(82)+acd202(79)
      acd202(82)=-acd202(10)*acd202(79)
      acd202(83)=-acd202(78)*acd202(1)*acd202(4)
      acd202(81)=acd202(81)*acd202(1)
      acd202(81)=-acd202(81)+acd202(83)
      acd202(83)=-acd202(11)*acd202(81)
      acd202(84)=acd202(69)*acd202(5)
      acd202(85)=acd202(25)*acd202(24)
      acd202(84)=acd202(84)-acd202(85)
      acd202(85)=acd202(84)*acd202(7)*acd202(4)
      acd202(86)=acd202(23)*acd202(22)
      acd202(87)=acd202(86)*acd202(21)
      acd202(88)=acd202(87)*acd202(7)
      acd202(85)=-acd202(88)+acd202(85)
      acd202(88)=-acd202(13)*acd202(85)
      acd202(89)=acd202(84)*acd202(6)*acd202(4)
      acd202(87)=acd202(87)*acd202(6)
      acd202(87)=-acd202(87)+acd202(89)
      acd202(89)=-acd202(14)*acd202(87)
      acd202(70)=acd202(89)+acd202(88)+acd202(82)+acd202(83)+acd202(71)+acd202(&
      &70)+acd202(77)+acd202(75)+acd202(73)
      acd202(71)=-acd202(33)-ninjaP2
      acd202(71)=acd202(66)*acd202(71)
      acd202(73)=acd202(41)*acd202(69)
      acd202(74)=acd202(48)*acd202(74)
      acd202(75)=acd202(51)*acd202(76)
      acd202(76)=acd202(55)*acd202(18)
      acd202(77)=acd202(57)*acd202(19)
      acd202(73)=acd202(77)+acd202(76)+acd202(75)+acd202(74)+acd202(73)
      acd202(73)=acd202(64)*acd202(73)
      acd202(74)=acd202(50)*acd202(21)
      acd202(74)=acd202(60)+acd202(74)
      acd202(74)=acd202(22)*acd202(74)
      acd202(75)=acd202(46)*acd202(69)
      acd202(72)=acd202(49)*acd202(72)
      acd202(76)=acd202(61)*acd202(24)
      acd202(72)=acd202(76)+acd202(72)+acd202(75)+acd202(74)
      acd202(72)=acd202(67)*acd202(72)
      acd202(74)=acd202(42)*acd202(1)
      acd202(75)=acd202(43)*acd202(2)
      acd202(74)=acd202(75)+acd202(74)
      acd202(74)=acd202(78)*acd202(74)
      acd202(75)=-acd202(40)*acd202(69)
      acd202(76)=acd202(47)*acd202(17)
      acd202(77)=acd202(53)*acd202(18)
      acd202(75)=acd202(77)+acd202(76)+acd202(75)
      acd202(65)=acd202(65)*acd202(75)
      acd202(69)=-acd202(37)*acd202(69)
      acd202(75)=acd202(52)*acd202(25)
      acd202(76)=acd202(54)*acd202(24)
      acd202(69)=acd202(76)+acd202(75)+acd202(69)
      acd202(68)=acd202(68)*acd202(69)
      acd202(69)=acd202(62)*acd202(6)
      acd202(75)=acd202(63)*acd202(7)
      acd202(69)=acd202(75)+acd202(69)
      acd202(69)=acd202(86)*acd202(69)
      acd202(75)=acd202(58)*acd202(1)
      acd202(76)=acd202(59)*acd202(2)
      acd202(75)=acd202(76)+acd202(75)
      acd202(75)=acd202(80)*acd202(75)
      acd202(76)=-ninjaP0-acd202(34)
      acd202(76)=acd202(66)*acd202(76)
      acd202(77)=-acd202(35)*acd202(79)
      acd202(78)=-acd202(36)*acd202(81)
      acd202(79)=-acd202(38)*acd202(85)
      acd202(80)=-acd202(39)*acd202(87)
      acd202(81)=-acd202(44)*acd202(6)*acd202(84)
      acd202(82)=-acd202(45)*acd202(7)*acd202(84)
      acd202(67)=acd202(67)*acd202(5)
      acd202(64)=acd202(64)*acd202(3)
      acd202(64)=-acd202(64)-acd202(67)
      acd202(64)=acd202(56)*acd202(64)
      acd202(64)=acd202(64)+acd202(82)+acd202(81)+acd202(80)+acd202(79)+acd202(&
      &77)+acd202(78)+acd202(73)+acd202(72)+acd202(68)+acd202(65)+acd202(76)+ac&
      &d202(74)+acd202(75)+acd202(69)
      brack(ninjaidxt0x0mu0)=acd202(64)
      brack(ninjaidxt0x0mu2)=-acd202(66)
      brack(ninjaidxt0x1mu0)=acd202(70)
      brack(ninjaidxt0x2mu0)=acd202(71)
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d202h0_qp_ninja_t2")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd202h0_qp
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
end module     p2_part21part21_part25part25part21_d202h0l132_qp
