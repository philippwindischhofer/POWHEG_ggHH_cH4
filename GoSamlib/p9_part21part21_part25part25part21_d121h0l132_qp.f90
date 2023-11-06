module     p9_part21part21_part25part25part21_d121h0l132_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d121h0l132_qp.f90
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
      use p9_part21part21_part25part25part21_abbrevd121h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(1) :: acd121
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=0.0_ki
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_color_qp
      use p9_part21part21_part25part25part21_abbrevd121h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(80) :: acd121
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd121(1)=dotproduct(ninjaA0,ninjaE3)
      acd121(2)=abb121(9)
      acd121(3)=dotproduct(ninjaE3,spvae1k2)
      acd121(4)=abb121(8)
      acd121(5)=dotproduct(ninjaE3,spvak2e2)
      acd121(6)=abb121(27)
      acd121(7)=dotproduct(ninjaE3,spvae5k2)
      acd121(8)=abb121(62)
      acd121(9)=dotproduct(ninjaE3,spvae2k2)
      acd121(10)=abb121(11)
      acd121(11)=dotproduct(ninjaE3,spvae2l3)
      acd121(12)=abb121(58)
      acd121(13)=dotproduct(ninjaE3,spvak1k2)
      acd121(14)=abb121(15)
      acd121(15)=dotproduct(ninjaE3,spvak5k2)
      acd121(16)=abb121(18)
      acd121(17)=dotproduct(ninjaE3,spval4e2)
      acd121(18)=abb121(33)
      acd121(19)=dotproduct(ninjaE3,spvae2l4)
      acd121(20)=abb121(31)
      acd121(21)=dotproduct(ninjaE3,spval3e2)
      acd121(22)=abb121(61)
      acd121(23)=dotproduct(ninjaE3,spvak2e1)
      acd121(24)=abb121(23)
      acd121(25)=dotproduct(ninjaE3,spvak2k1)
      acd121(26)=abb121(49)
      acd121(27)=dotproduct(ninjaE3,spvak2e5)
      acd121(28)=abb121(65)
      acd121(29)=dotproduct(ninjaE3,spvak2k5)
      acd121(30)=abb121(64)
      acd121(31)=abb121(7)
      acd121(32)=abb121(22)
      acd121(33)=abb121(10)
      acd121(34)=abb121(14)
      acd121(35)=abb121(12)
      acd121(36)=abb121(43)
      acd121(37)=abb121(60)
      acd121(38)=abb121(21)
      acd121(39)=abb121(59)
      acd121(40)=abb121(19)
      acd121(41)=abb121(28)
      acd121(42)=abb121(45)
      acd121(43)=abb121(70)
      acd121(44)=abb121(71)
      acd121(45)=abb121(72)
      acd121(46)=abb121(35)
      acd121(47)=abb121(34)
      acd121(48)=abb121(68)
      acd121(49)=abb121(56)
      acd121(50)=abb121(69)
      acd121(51)=abb121(41)
      acd121(52)=abb121(55)
      acd121(53)=abb121(37)
      acd121(54)=dotproduct(ninjaE3,spvae1e2)
      acd121(55)=dotproduct(ninjaE3,spvae2e1)
      acd121(56)=abb121(38)
      acd121(57)=dotproduct(ninjaE3,spvak1e2)
      acd121(58)=dotproduct(ninjaE3,spvae2k1)
      acd121(59)=abb121(50)
      acd121(60)=dotproduct(ninjaE3,spvak5e2)
      acd121(61)=dotproduct(ninjaE3,spvae2k5)
      acd121(62)=dotproduct(ninjaE3,spvae2e5)
      acd121(63)=dotproduct(ninjaE3,spvae5e2)
      acd121(64)=abb121(52)
      acd121(65)=acd121(29)*acd121(30)
      acd121(66)=acd121(27)*acd121(28)
      acd121(67)=acd121(25)*acd121(26)
      acd121(68)=acd121(23)*acd121(24)
      acd121(69)=acd121(19)*acd121(20)
      acd121(70)=acd121(17)*acd121(18)
      acd121(71)=acd121(15)*acd121(16)
      acd121(72)=acd121(13)*acd121(14)
      acd121(73)=acd121(7)*acd121(8)
      acd121(74)=acd121(3)*acd121(4)
      acd121(75)=acd121(21)*acd121(22)
      acd121(76)=acd121(11)*acd121(12)
      acd121(77)=acd121(9)*acd121(10)
      acd121(78)=acd121(5)*acd121(6)
      acd121(79)=2.0_ki*acd121(1)
      acd121(80)=acd121(2)*acd121(79)
      acd121(65)=acd121(80)+acd121(78)+acd121(77)+acd121(76)+acd121(75)+acd121(&
      &74)+acd121(73)+acd121(72)+acd121(71)+acd121(70)+acd121(69)+acd121(68)+ac&
      &d121(67)+acd121(65)+acd121(66)
      acd121(65)=acd121(65)*acd121(79)
      acd121(66)=acd121(19)*acd121(38)
      acd121(67)=acd121(15)*acd121(37)
      acd121(68)=acd121(13)*acd121(36)
      acd121(69)=acd121(7)*acd121(33)
      acd121(70)=acd121(3)*acd121(31)
      acd121(71)=acd121(11)*acd121(35)
      acd121(72)=acd121(9)*acd121(34)
      acd121(66)=acd121(72)+acd121(71)+acd121(70)+acd121(69)+acd121(68)+acd121(&
      &66)+acd121(67)
      acd121(66)=acd121(5)*acd121(66)
      acd121(67)=acd121(29)*acd121(45)
      acd121(68)=acd121(27)*acd121(44)
      acd121(69)=acd121(25)*acd121(43)
      acd121(70)=acd121(23)*acd121(42)
      acd121(71)=acd121(17)*acd121(40)
      acd121(72)=acd121(21)*acd121(41)
      acd121(67)=acd121(72)+acd121(71)+acd121(70)+acd121(69)+acd121(67)+acd121(&
      &68)
      acd121(67)=acd121(9)*acd121(67)
      acd121(68)=acd121(19)*acd121(53)
      acd121(69)=acd121(15)*acd121(52)
      acd121(70)=acd121(13)*acd121(51)
      acd121(71)=acd121(7)*acd121(39)
      acd121(72)=acd121(3)*acd121(32)
      acd121(68)=acd121(72)+acd121(71)+acd121(70)+acd121(68)+acd121(69)
      acd121(68)=acd121(21)*acd121(68)
      acd121(69)=acd121(29)*acd121(50)
      acd121(70)=acd121(27)*acd121(49)
      acd121(71)=acd121(25)*acd121(48)
      acd121(72)=acd121(23)*acd121(47)
      acd121(73)=acd121(17)*acd121(46)
      acd121(69)=acd121(73)+acd121(72)+acd121(71)+acd121(69)+acd121(70)
      acd121(69)=acd121(11)*acd121(69)
      acd121(70)=-acd121(60)*acd121(61)
      acd121(71)=-acd121(57)*acd121(58)
      acd121(70)=acd121(70)+acd121(71)
      acd121(70)=acd121(59)*acd121(70)
      acd121(71)=acd121(62)*acd121(63)*acd121(64)
      acd121(72)=acd121(54)*acd121(55)*acd121(56)
      acd121(65)=acd121(65)+acd121(66)+acd121(67)+acd121(69)+acd121(68)+acd121(&
      &70)+acd121(71)+acd121(72)
      brack(ninjaidxt0x0mu0)=acd121(65)
      brack(ninjaidxt0x0mu2)=0.0_ki
      brack(ninjaidxt0x1mu0)=0.0_ki
      brack(ninjaidxt0x2mu0)=0.0_ki
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p9_part21part21_part25part25part21_d121h0_qp_ninja_t2")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p9_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_abbrevd121h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a0, a1, b, c
      complex(ki), dimension(0:2), intent(in) :: param
      complex(ki), dimension(4) :: vecA0, vecA1, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
	     vecA0(1:4) = + a0(0:3)
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
end module     p9_part21part21_part25part25part21_d121h0l132_qp
