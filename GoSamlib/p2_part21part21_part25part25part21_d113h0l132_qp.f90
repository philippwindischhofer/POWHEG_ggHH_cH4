module     p2_part21part21_part25part25part21_d113h0l132_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d113h0l132_qp.f90
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
      use p2_part21part21_part25part25part21_abbrevd113h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(11) :: acd113
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd113(1)=dotproduct(ninjaA0,ninjaE3)
      acd113(2)=dotproduct(ninjaE3,spvae1e5)
      acd113(3)=dotproduct(ninjaE3,spvae5e1)
      acd113(4)=abb113(47)
      acd113(5)=dotproduct(ninjaE3,spvak1e5)
      acd113(6)=dotproduct(ninjaE3,spvae5k1)
      acd113(7)=abb113(79)
      acd113(8)=dotproduct(ninjaE3,spvak2e5)
      acd113(9)=dotproduct(ninjaE3,spvae5k2)
      acd113(10)=acd113(6)*acd113(5)
      acd113(11)=-acd113(9)*acd113(8)
      acd113(10)=acd113(11)+acd113(10)
      acd113(10)=acd113(10)*acd113(7)
      acd113(11)=acd113(4)*acd113(3)*acd113(2)
      acd113(10)=acd113(11)+acd113(10)
      acd113(10)=2.0_ki*acd113(1)*acd113(10)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd113(10)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd113h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(75) :: acd113
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd113(1)=dotproduct(ninjaE3,spvae5k1)
      acd113(2)=dotproduct(ninjaE3,spvak1e5)
      acd113(3)=abb113(79)
      acd113(4)=dotproduct(ninjaE3,spvae1e5)
      acd113(5)=dotproduct(ninjaE3,spvae5e1)
      acd113(6)=abb113(47)
      acd113(7)=dotproduct(ninjaE3,spvae5k2)
      acd113(8)=dotproduct(ninjaE3,spvak2e5)
      acd113(9)=dotproduct(ninjaA0,ninjaA1)
      acd113(10)=dotproduct(ninjaA0,ninjaE3)
      acd113(11)=dotproduct(ninjaA1,spvae5k1)
      acd113(12)=dotproduct(ninjaA1,spvae1e5)
      acd113(13)=dotproduct(ninjaA1,spvae5e1)
      acd113(14)=dotproduct(ninjaA1,spvak1e5)
      acd113(15)=dotproduct(ninjaA1,spvae5k2)
      acd113(16)=dotproduct(ninjaA1,spvak2e5)
      acd113(17)=dotproduct(ninjaA1,ninjaA1)
      acd113(18)=dotproduct(k2,ninjaE3)
      acd113(19)=abb113(49)
      acd113(20)=abb113(54)
      acd113(21)=abb113(51)
      acd113(22)=abb113(45)
      acd113(23)=dotproduct(ninjaA0,ninjaA0)
      acd113(24)=dotproduct(ninjaA0,spvae5k1)
      acd113(25)=dotproduct(ninjaA0,spvae1e5)
      acd113(26)=dotproduct(ninjaA0,spvae5e1)
      acd113(27)=dotproduct(ninjaA0,spvak1e5)
      acd113(28)=dotproduct(ninjaA0,spvae5k2)
      acd113(29)=dotproduct(ninjaA0,spvak2e5)
      acd113(30)=abb113(44)
      acd113(31)=abb113(16)
      acd113(32)=abb113(20)
      acd113(33)=abb113(37)
      acd113(34)=abb113(57)
      acd113(35)=abb113(61)
      acd113(36)=dotproduct(ninjaE3,spvak2l3)
      acd113(37)=abb113(7)
      acd113(38)=abb113(15)
      acd113(39)=abb113(22)
      acd113(40)=abb113(69)
      acd113(41)=dotproduct(ninjaE3,spval4k2)
      acd113(42)=abb113(75)
      acd113(43)=dotproduct(ninjaE3,spval4l3)
      acd113(44)=abb113(56)
      acd113(45)=dotproduct(ninjaE3,spval3l4)
      acd113(46)=abb113(9)
      acd113(47)=abb113(43)
      acd113(48)=abb113(59)
      acd113(49)=dotproduct(ninjaE3,spvak2l4)
      acd113(50)=abb113(14)
      acd113(51)=abb113(29)
      acd113(52)=dotproduct(ninjaE3,spval3k2)
      acd113(53)=abb113(18)
      acd113(54)=abb113(19)
      acd113(55)=abb113(28)
      acd113(56)=abb113(66)
      acd113(57)=abb113(31)
      acd113(58)=abb113(34)
      acd113(59)=abb113(83)
      acd113(60)=acd113(6)*acd113(5)
      acd113(61)=acd113(60)*acd113(4)
      acd113(62)=acd113(8)*acd113(3)
      acd113(63)=acd113(62)*acd113(7)
      acd113(64)=acd113(2)*acd113(3)
      acd113(65)=acd113(64)*acd113(1)
      acd113(61)=acd113(65)+acd113(61)-acd113(63)
      acd113(63)=acd113(11)*acd113(64)
      acd113(65)=acd113(12)*acd113(60)
      acd113(66)=acd113(6)*acd113(4)
      acd113(67)=acd113(13)*acd113(66)
      acd113(68)=acd113(1)*acd113(3)
      acd113(69)=acd113(14)*acd113(68)
      acd113(70)=-acd113(15)*acd113(62)
      acd113(71)=acd113(7)*acd113(3)
      acd113(72)=-acd113(16)*acd113(71)
      acd113(63)=acd113(72)+acd113(70)+acd113(69)+acd113(67)+acd113(65)+acd113(&
      &63)
      acd113(65)=2.0_ki*acd113(10)
      acd113(63)=acd113(65)*acd113(63)
      acd113(67)=2.0_ki*acd113(9)+ninjaP1
      acd113(67)=acd113(61)*acd113(67)
      acd113(63)=acd113(63)+acd113(67)
      acd113(67)=acd113(17)+ninjaP2
      acd113(67)=acd113(61)*acd113(67)
      acd113(64)=acd113(24)*acd113(64)
      acd113(60)=acd113(25)*acd113(60)
      acd113(66)=acd113(26)*acd113(66)
      acd113(68)=acd113(27)*acd113(68)
      acd113(62)=-acd113(28)*acd113(62)
      acd113(69)=-acd113(29)*acd113(71)
      acd113(70)=acd113(30)*acd113(1)
      acd113(71)=acd113(31)*acd113(4)
      acd113(72)=acd113(32)*acd113(5)
      acd113(73)=acd113(33)*acd113(2)
      acd113(74)=acd113(34)*acd113(7)
      acd113(75)=acd113(35)*acd113(8)
      acd113(60)=acd113(75)+acd113(74)+acd113(73)+acd113(72)+acd113(71)+acd113(&
      &70)+acd113(69)+acd113(62)+acd113(68)+acd113(66)+acd113(60)+acd113(64)
      acd113(60)=acd113(65)*acd113(60)
      acd113(62)=acd113(19)*acd113(1)
      acd113(64)=acd113(20)*acd113(4)
      acd113(65)=acd113(21)*acd113(5)
      acd113(66)=acd113(22)*acd113(2)
      acd113(62)=acd113(66)+acd113(65)+acd113(64)+acd113(62)
      acd113(62)=acd113(18)*acd113(62)
      acd113(64)=acd113(23)+ninjaP0
      acd113(64)=acd113(61)*acd113(64)
      acd113(65)=acd113(37)*acd113(1)
      acd113(66)=acd113(38)*acd113(5)
      acd113(68)=acd113(39)*acd113(7)
      acd113(65)=acd113(68)+acd113(66)+acd113(65)
      acd113(65)=acd113(36)*acd113(65)
      acd113(66)=acd113(44)*acd113(1)
      acd113(68)=acd113(56)*acd113(5)
      acd113(69)=acd113(59)*acd113(7)
      acd113(66)=acd113(69)+acd113(68)+acd113(66)
      acd113(66)=acd113(43)*acd113(66)
      acd113(68)=acd113(46)*acd113(4)
      acd113(69)=acd113(47)*acd113(2)
      acd113(70)=acd113(48)*acd113(8)
      acd113(68)=acd113(70)+acd113(69)+acd113(68)
      acd113(68)=acd113(45)*acd113(68)
      acd113(69)=acd113(53)*acd113(4)
      acd113(70)=acd113(57)*acd113(2)
      acd113(71)=acd113(58)*acd113(8)
      acd113(69)=acd113(71)+acd113(70)+acd113(69)
      acd113(69)=acd113(52)*acd113(69)
      acd113(70)=acd113(42)*acd113(1)
      acd113(71)=acd113(55)*acd113(5)
      acd113(70)=acd113(71)+acd113(70)
      acd113(70)=acd113(41)*acd113(70)
      acd113(71)=acd113(50)*acd113(4)
      acd113(72)=acd113(54)*acd113(2)
      acd113(71)=acd113(72)+acd113(71)
      acd113(71)=acd113(49)*acd113(71)
      acd113(72)=acd113(2)*acd113(1)
      acd113(73)=acd113(8)*acd113(7)
      acd113(72)=-acd113(72)+acd113(73)
      acd113(72)=acd113(40)*acd113(72)
      acd113(73)=acd113(51)*acd113(5)*acd113(4)
      acd113(60)=acd113(73)+acd113(72)+acd113(60)+acd113(62)+acd113(69)+acd113(&
      &68)+acd113(66)+acd113(65)+acd113(64)+acd113(71)+acd113(70)
      brack(ninjaidxt0x0mu0)=acd113(60)
      brack(ninjaidxt0x0mu2)=acd113(61)
      brack(ninjaidxt0x1mu0)=acd113(63)
      brack(ninjaidxt0x2mu0)=acd113(67)
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d113h0_qp_ninja_t2")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd113h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a0, a1, b, c
      complex(ki), dimension(0:2), intent(in) :: param
      complex(ki), dimension(4) :: vecA0, vecA1, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      qshift = k5+k4
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
end module     p2_part21part21_part25part25part21_d113h0l132_qp
