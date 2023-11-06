module     p9_part21part21_part25part25part21_d109h0l132_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d109h0l132_qp.f90
   ! generator: buildfortran_tn3.py
   use p9_part21part21_part25part25part21_config, only: ki => ki_qp
   use p9_part21part21_part25part25part21_util_qp, only: cond_t, d => metric_te&
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
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_color_qp
      use p9_part21part21_part25part25part21_abbrevd109h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(11) :: acd109
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd109(1)=dotproduct(ninjaA0,ninjaE3)
      acd109(2)=dotproduct(ninjaE3,spvae1e5)
      acd109(3)=dotproduct(ninjaE3,spvae5e1)
      acd109(4)=abb109(47)
      acd109(5)=dotproduct(ninjaE3,spvak1e5)
      acd109(6)=dotproduct(ninjaE3,spvae5k1)
      acd109(7)=abb109(79)
      acd109(8)=dotproduct(ninjaE3,spvak2e5)
      acd109(9)=dotproduct(ninjaE3,spvae5k2)
      acd109(10)=acd109(6)*acd109(5)
      acd109(11)=-acd109(9)*acd109(8)
      acd109(10)=acd109(11)+acd109(10)
      acd109(10)=acd109(10)*acd109(7)
      acd109(11)=acd109(4)*acd109(3)*acd109(2)
      acd109(10)=acd109(11)+acd109(10)
      acd109(10)=2.0_ki*acd109(1)*acd109(10)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd109(10)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_color_qp
      use p9_part21part21_part25part25part21_abbrevd109h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(75) :: acd109
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd109(1)=dotproduct(ninjaE3,spvae5e1)
      acd109(2)=dotproduct(ninjaE3,spvae1e5)
      acd109(3)=abb109(47)
      acd109(4)=dotproduct(ninjaE3,spvae5k1)
      acd109(5)=dotproduct(ninjaE3,spvak1e5)
      acd109(6)=abb109(79)
      acd109(7)=dotproduct(ninjaE3,spvae5k2)
      acd109(8)=dotproduct(ninjaE3,spvak2e5)
      acd109(9)=dotproduct(ninjaA0,ninjaA1)
      acd109(10)=dotproduct(ninjaA0,ninjaE3)
      acd109(11)=dotproduct(ninjaA1,spvae5e1)
      acd109(12)=dotproduct(ninjaA1,spvae5k1)
      acd109(13)=dotproduct(ninjaA1,spvae1e5)
      acd109(14)=dotproduct(ninjaA1,spvak1e5)
      acd109(15)=dotproduct(ninjaA1,spvae5k2)
      acd109(16)=dotproduct(ninjaA1,spvak2e5)
      acd109(17)=dotproduct(ninjaA1,ninjaA1)
      acd109(18)=dotproduct(k2,ninjaE3)
      acd109(19)=abb109(51)
      acd109(20)=abb109(49)
      acd109(21)=abb109(54)
      acd109(22)=abb109(45)
      acd109(23)=dotproduct(ninjaA0,ninjaA0)
      acd109(24)=dotproduct(ninjaA0,spvae5e1)
      acd109(25)=dotproduct(ninjaA0,spvae5k1)
      acd109(26)=dotproduct(ninjaA0,spvae1e5)
      acd109(27)=dotproduct(ninjaA0,spvak1e5)
      acd109(28)=dotproduct(ninjaA0,spvae5k2)
      acd109(29)=dotproduct(ninjaA0,spvak2e5)
      acd109(30)=abb109(20)
      acd109(31)=abb109(44)
      acd109(32)=abb109(16)
      acd109(33)=abb109(37)
      acd109(34)=abb109(61)
      acd109(35)=abb109(63)
      acd109(36)=dotproduct(ninjaE3,spval3l4)
      acd109(37)=abb109(8)
      acd109(38)=abb109(9)
      acd109(39)=abb109(46)
      acd109(40)=abb109(29)
      acd109(41)=dotproduct(ninjaE3,spvak2l4)
      acd109(42)=abb109(22)
      acd109(43)=dotproduct(ninjaE3,spval3k2)
      acd109(44)=abb109(25)
      acd109(45)=abb109(69)
      acd109(46)=abb109(55)
      acd109(47)=abb109(31)
      acd109(48)=dotproduct(ninjaE3,spvak2l3)
      acd109(49)=abb109(12)
      acd109(50)=abb109(14)
      acd109(51)=dotproduct(ninjaE3,spval4k2)
      acd109(52)=abb109(30)
      acd109(53)=dotproduct(ninjaE3,spval4l3)
      acd109(54)=abb109(59)
      acd109(55)=abb109(87)
      acd109(56)=abb109(57)
      acd109(57)=abb109(60)
      acd109(58)=abb109(88)
      acd109(59)=abb109(83)
      acd109(60)=acd109(3)*acd109(2)
      acd109(61)=acd109(60)*acd109(1)
      acd109(62)=acd109(8)*acd109(6)
      acd109(63)=acd109(62)*acd109(7)
      acd109(64)=acd109(5)*acd109(6)
      acd109(65)=acd109(64)*acd109(4)
      acd109(61)=acd109(65)+acd109(61)-acd109(63)
      acd109(63)=acd109(11)*acd109(60)
      acd109(65)=acd109(12)*acd109(64)
      acd109(66)=acd109(3)*acd109(1)
      acd109(67)=acd109(13)*acd109(66)
      acd109(68)=acd109(4)*acd109(6)
      acd109(69)=acd109(14)*acd109(68)
      acd109(70)=-acd109(15)*acd109(62)
      acd109(71)=acd109(7)*acd109(6)
      acd109(72)=-acd109(16)*acd109(71)
      acd109(63)=acd109(72)+acd109(70)+acd109(69)+acd109(67)+acd109(65)+acd109(&
      &63)
      acd109(65)=2.0_ki*acd109(10)
      acd109(63)=acd109(65)*acd109(63)
      acd109(67)=2.0_ki*acd109(9)+ninjaP1
      acd109(67)=acd109(61)*acd109(67)
      acd109(63)=acd109(63)+acd109(67)
      acd109(67)=acd109(17)+ninjaP2
      acd109(67)=acd109(61)*acd109(67)
      acd109(60)=acd109(24)*acd109(60)
      acd109(64)=acd109(25)*acd109(64)
      acd109(66)=acd109(26)*acd109(66)
      acd109(68)=acd109(27)*acd109(68)
      acd109(62)=-acd109(28)*acd109(62)
      acd109(69)=-acd109(29)*acd109(71)
      acd109(70)=acd109(30)*acd109(1)
      acd109(71)=acd109(31)*acd109(4)
      acd109(72)=acd109(32)*acd109(2)
      acd109(73)=acd109(33)*acd109(5)
      acd109(74)=acd109(34)*acd109(7)
      acd109(75)=acd109(35)*acd109(8)
      acd109(60)=acd109(75)+acd109(74)+acd109(73)+acd109(72)+acd109(71)+acd109(&
      &70)+acd109(69)+acd109(62)+acd109(68)+acd109(66)+acd109(64)+acd109(60)
      acd109(60)=acd109(65)*acd109(60)
      acd109(62)=acd109(19)*acd109(1)
      acd109(64)=acd109(20)*acd109(4)
      acd109(65)=acd109(21)*acd109(2)
      acd109(66)=acd109(22)*acd109(5)
      acd109(62)=acd109(66)+acd109(65)+acd109(64)+acd109(62)
      acd109(62)=acd109(18)*acd109(62)
      acd109(64)=acd109(23)+ninjaP0
      acd109(64)=acd109(61)*acd109(64)
      acd109(65)=acd109(37)*acd109(1)
      acd109(66)=acd109(38)*acd109(4)
      acd109(68)=acd109(39)*acd109(7)
      acd109(65)=acd109(68)+acd109(66)+acd109(65)
      acd109(65)=acd109(36)*acd109(65)
      acd109(66)=acd109(42)*acd109(1)
      acd109(68)=acd109(46)*acd109(4)
      acd109(69)=acd109(57)*acd109(7)
      acd109(66)=acd109(69)+acd109(68)+acd109(66)
      acd109(66)=acd109(41)*acd109(66)
      acd109(68)=acd109(52)*acd109(2)
      acd109(69)=acd109(55)*acd109(5)
      acd109(70)=acd109(58)*acd109(8)
      acd109(68)=acd109(70)+acd109(69)+acd109(68)
      acd109(68)=acd109(51)*acd109(68)
      acd109(69)=acd109(54)*acd109(2)
      acd109(70)=acd109(56)*acd109(5)
      acd109(71)=acd109(59)*acd109(8)
      acd109(69)=acd109(71)+acd109(70)+acd109(69)
      acd109(69)=acd109(53)*acd109(69)
      acd109(70)=acd109(44)*acd109(1)
      acd109(71)=acd109(47)*acd109(4)
      acd109(70)=acd109(71)+acd109(70)
      acd109(70)=acd109(43)*acd109(70)
      acd109(71)=acd109(49)*acd109(2)
      acd109(72)=acd109(50)*acd109(5)
      acd109(71)=acd109(72)+acd109(71)
      acd109(71)=acd109(48)*acd109(71)
      acd109(72)=acd109(5)*acd109(4)
      acd109(73)=acd109(8)*acd109(7)
      acd109(72)=-acd109(72)+acd109(73)
      acd109(72)=acd109(45)*acd109(72)
      acd109(73)=acd109(40)*acd109(2)*acd109(1)
      acd109(60)=acd109(73)+acd109(72)+acd109(60)+acd109(62)+acd109(69)+acd109(&
      &68)+acd109(66)+acd109(65)+acd109(64)+acd109(71)+acd109(70)
      brack(ninjaidxt0x0mu0)=acd109(60)
      brack(ninjaidxt0x0mu2)=acd109(61)
      brack(ninjaidxt0x1mu0)=acd109(63)
      brack(ninjaidxt0x2mu0)=acd109(67)
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p9_part21part21_part25part25part21_d109h0_qp_ninja_t2")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p9_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_abbrevd109h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a0, a1, b, c
      complex(ki), dimension(0:2), intent(in) :: param
      complex(ki), dimension(4) :: vecA0, vecA1, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      qshift = k3+k5
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
end module     p9_part21part21_part25part25part21_d109h0l132_qp
