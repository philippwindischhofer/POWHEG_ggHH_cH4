module     p2_part21part21_part25part25part21_d188h0l132
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d188h0l132.f90
   ! generator: buildfortran_tn3.py
   use p2_part21part21_part25part25part21_config, only: ki
   use p2_part21part21_part25part25part21_util, only: cond_t, d => metric_tensor
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
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_color
      use p2_part21part21_part25part25part21_abbrevd188h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(1) :: acd188
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
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_color
      use p2_part21part21_part25part25part21_abbrevd188h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(70) :: acd188
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd188(1)=dotproduct(ninjaA0,ninjaE3)
      acd188(2)=dotproduct(ninjaE3,spvae2e1)
      acd188(3)=abb188(82)
      acd188(4)=dotproduct(ninjaE3,spvae1e2)
      acd188(5)=abb188(52)
      acd188(6)=dotproduct(ninjaE3,spvae1k1)
      acd188(7)=abb188(81)
      acd188(8)=dotproduct(ninjaE3,spvak2e2)
      acd188(9)=abb188(89)
      acd188(10)=dotproduct(ninjaE3,spvae1k2)
      acd188(11)=abb188(70)
      acd188(12)=dotproduct(ninjaE3,spvak2e5)
      acd188(13)=abb188(84)
      acd188(14)=dotproduct(ninjaE3,spvae5k1)
      acd188(15)=abb188(86)
      acd188(16)=dotproduct(ninjaE3,spval3e2)
      acd188(17)=abb188(83)
      acd188(18)=dotproduct(ninjaE3,spvae1l4)
      acd188(19)=abb188(62)
      acd188(20)=dotproduct(ninjaE3,spvak1e1)
      acd188(21)=abb188(19)
      acd188(22)=dotproduct(ninjaE3,spvak2e1)
      acd188(23)=abb188(87)
      acd188(24)=dotproduct(ninjaE3,spvae5k2)
      acd188(25)=abb188(100)
      acd188(26)=dotproduct(ninjaE3,spvae2l3)
      acd188(27)=abb188(45)
      acd188(28)=dotproduct(ninjaE3,spvak1e5)
      acd188(29)=abb188(98)
      acd188(30)=dotproduct(ninjaE3,spval4e1)
      acd188(31)=abb188(65)
      acd188(32)=dotproduct(ninjaE3,spvae2k2)
      acd188(33)=abb188(107)
      acd188(34)=abb188(7)
      acd188(35)=abb188(80)
      acd188(36)=abb188(34)
      acd188(37)=abb188(64)
      acd188(38)=abb188(112)
      acd188(39)=abb188(12)
      acd188(40)=abb188(50)
      acd188(41)=abb188(33)
      acd188(42)=abb188(103)
      acd188(43)=abb188(32)
      acd188(44)=abb188(79)
      acd188(45)=dotproduct(ninjaE3,spvae1e5)
      acd188(46)=dotproduct(ninjaE3,spvae5e2)
      acd188(47)=abb188(88)
      acd188(48)=abb188(71)
      acd188(49)=abb188(31)
      acd188(50)=abb188(85)
      acd188(51)=abb188(22)
      acd188(52)=abb188(44)
      acd188(53)=abb188(77)
      acd188(54)=abb188(102)
      acd188(55)=abb188(46)
      acd188(56)=abb188(99)
      acd188(57)=abb188(72)
      acd188(58)=abb188(108)
      acd188(59)=dotproduct(ninjaE3,spvae5e1)
      acd188(60)=dotproduct(ninjaE3,spvae2e5)
      acd188(61)=acd188(22)*acd188(23)
      acd188(62)=acd188(32)*acd188(33)
      acd188(63)=acd188(30)*acd188(31)
      acd188(64)=acd188(28)*acd188(29)
      acd188(65)=acd188(20)*acd188(21)
      acd188(66)=acd188(26)*acd188(27)
      acd188(67)=acd188(24)*acd188(25)
      acd188(68)=2.0_ki*acd188(1)
      acd188(69)=acd188(5)*acd188(68)
      acd188(61)=acd188(69)+acd188(67)+acd188(66)+acd188(65)+acd188(64)+acd188(&
      &63)+acd188(61)+acd188(62)
      acd188(61)=acd188(61)*acd188(68)
      acd188(62)=acd188(22)*acd188(52)
      acd188(63)=acd188(30)*acd188(56)
      acd188(64)=acd188(28)*acd188(55)
      acd188(65)=acd188(20)*acd188(49)
      acd188(62)=acd188(65)+acd188(64)+acd188(62)+acd188(63)
      acd188(62)=acd188(26)*acd188(62)
      acd188(63)=acd188(22)*acd188(51)
      acd188(64)=acd188(30)*acd188(54)
      acd188(65)=acd188(28)*acd188(53)
      acd188(66)=acd188(20)*acd188(48)
      acd188(63)=acd188(66)+acd188(65)+acd188(63)+acd188(64)
      acd188(63)=acd188(24)*acd188(63)
      acd188(64)=acd188(30)*acd188(58)
      acd188(65)=acd188(28)*acd188(57)
      acd188(66)=acd188(20)*acd188(50)
      acd188(64)=acd188(66)+acd188(65)+acd188(64)
      acd188(64)=acd188(32)*acd188(64)
      acd188(65)=-acd188(47)*acd188(59)*acd188(60)
      acd188(61)=acd188(61)+acd188(63)+acd188(62)+acd188(65)+acd188(64)
      acd188(61)=acd188(4)*acd188(61)
      acd188(62)=acd188(10)*acd188(11)
      acd188(63)=acd188(18)*acd188(19)
      acd188(64)=acd188(14)*acd188(15)
      acd188(65)=acd188(8)*acd188(9)
      acd188(66)=acd188(6)*acd188(7)
      acd188(67)=acd188(16)*acd188(17)
      acd188(69)=acd188(12)*acd188(13)
      acd188(70)=acd188(3)*acd188(68)
      acd188(62)=acd188(70)+acd188(69)+acd188(67)+acd188(66)+acd188(65)+acd188(&
      &64)+acd188(62)+acd188(63)
      acd188(62)=acd188(62)*acd188(68)
      acd188(63)=acd188(10)*acd188(40)
      acd188(64)=acd188(18)*acd188(44)
      acd188(65)=acd188(14)*acd188(43)
      acd188(66)=acd188(6)*acd188(36)
      acd188(63)=acd188(66)+acd188(65)+acd188(63)+acd188(64)
      acd188(63)=acd188(16)*acd188(63)
      acd188(64)=acd188(10)*acd188(39)
      acd188(65)=acd188(18)*acd188(42)
      acd188(66)=acd188(14)*acd188(41)
      acd188(67)=acd188(6)*acd188(35)
      acd188(64)=acd188(67)+acd188(66)+acd188(64)+acd188(65)
      acd188(64)=acd188(12)*acd188(64)
      acd188(65)=acd188(18)*acd188(38)
      acd188(66)=acd188(14)*acd188(37)
      acd188(67)=acd188(6)*acd188(34)
      acd188(65)=acd188(67)+acd188(65)+acd188(66)
      acd188(65)=acd188(8)*acd188(65)
      acd188(66)=-acd188(47)*acd188(45)*acd188(46)
      acd188(62)=acd188(62)+acd188(64)+acd188(63)+acd188(66)+acd188(65)
      acd188(62)=acd188(2)*acd188(62)
      acd188(61)=acd188(61)+acd188(62)
      brack(ninjaidxt0x0mu0)=acd188(61)
      brack(ninjaidxt0x0mu2)=0.0_ki
      brack(ninjaidxt0x1mu0)=0.0_ki
      brack(ninjaidxt0x2mu0)=0.0_ki
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d188h0_ninja_t2")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1, only: epspow
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_abbrevd188h0
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
end module     p2_part21part21_part25part25part21_d188h0l132
