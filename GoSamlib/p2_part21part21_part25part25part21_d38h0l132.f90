module     p2_part21part21_part25part25part21_d38h0l132
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d38h0l132.f90
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
      use p2_part21part21_part25part25part21_abbrevd38h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(14) :: acd38
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd38(1)=dotproduct(ninjaE3,spvae2e5)
      acd38(2)=dotproduct(ninjaE3,spvae5e2)
      acd38(3)=abb38(4)
      acd38(4)=dotproduct(ninjaE3,spvae1e2)
      acd38(5)=dotproduct(ninjaE3,spvae2e1)
      acd38(6)=abb38(7)
      acd38(7)=dotproduct(ninjaE3,spvak1e2)
      acd38(8)=dotproduct(ninjaE3,spvae2k1)
      acd38(9)=abb38(10)
      acd38(10)=dotproduct(ninjaE3,spvak5e2)
      acd38(11)=dotproduct(ninjaE3,spvae2k5)
      acd38(12)=-acd38(8)*acd38(7)
      acd38(13)=-acd38(11)*acd38(10)
      acd38(12)=acd38(13)+acd38(12)
      acd38(12)=acd38(9)*acd38(12)
      acd38(13)=acd38(3)*acd38(2)*acd38(1)
      acd38(14)=acd38(6)*acd38(5)*acd38(4)
      acd38(12)=acd38(13)+acd38(14)+acd38(12)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd38(12)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_color
      use p2_part21part21_part25part25part21_abbrevd38h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(74) :: acd38
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd38(1)=dotproduct(ninjaA1,spvae2e5)
      acd38(2)=dotproduct(ninjaE3,spvae5e2)
      acd38(3)=abb38(4)
      acd38(4)=dotproduct(ninjaA1,spvae5e2)
      acd38(5)=dotproduct(ninjaE3,spvae2e5)
      acd38(6)=dotproduct(ninjaA1,spvak1e2)
      acd38(7)=dotproduct(ninjaE3,spvae2k1)
      acd38(8)=abb38(10)
      acd38(9)=dotproduct(ninjaA1,spvae1e2)
      acd38(10)=dotproduct(ninjaE3,spvae2e1)
      acd38(11)=abb38(7)
      acd38(12)=dotproduct(ninjaA1,spvae2e1)
      acd38(13)=dotproduct(ninjaE3,spvae1e2)
      acd38(14)=dotproduct(ninjaA1,spvae2k1)
      acd38(15)=dotproduct(ninjaE3,spvak1e2)
      acd38(16)=dotproduct(ninjaA1,spvak5e2)
      acd38(17)=dotproduct(ninjaE3,spvae2k5)
      acd38(18)=dotproduct(ninjaA1,spvae2k5)
      acd38(19)=dotproduct(ninjaE3,spvak5e2)
      acd38(20)=dotproduct(ninjaA0,ninjaE3)
      acd38(21)=abb38(6)
      acd38(22)=dotproduct(ninjaA0,spvae2e5)
      acd38(23)=dotproduct(ninjaA0,spvae5e2)
      acd38(24)=dotproduct(ninjaA0,spvak1e2)
      acd38(25)=dotproduct(ninjaA0,spvae1e2)
      acd38(26)=dotproduct(ninjaA0,spvae2e1)
      acd38(27)=dotproduct(ninjaA0,spvae2k1)
      acd38(28)=dotproduct(ninjaA0,spvak5e2)
      acd38(29)=dotproduct(ninjaA0,spvae2k5)
      acd38(30)=abb38(13)
      acd38(31)=abb38(14)
      acd38(32)=abb38(5)
      acd38(33)=abb38(17)
      acd38(34)=abb38(20)
      acd38(35)=dotproduct(ninjaE3,spvak1k2)
      acd38(36)=abb38(8)
      acd38(37)=abb38(18)
      acd38(38)=abb38(22)
      acd38(39)=abb38(16)
      acd38(40)=dotproduct(ninjaE3,spvak5k2)
      acd38(41)=abb38(11)
      acd38(42)=dotproduct(ninjaE3,spvak2k1)
      acd38(43)=abb38(15)
      acd38(44)=dotproduct(ninjaE3,spvae1k2)
      acd38(45)=abb38(19)
      acd38(46)=dotproduct(ninjaE3,spvak2e1)
      acd38(47)=abb38(21)
      acd38(48)=dotproduct(ninjaE3,spvak2k5)
      acd38(49)=abb38(23)
      acd38(50)=dotproduct(ninjaE3,spvae5k2)
      acd38(51)=abb38(24)
      acd38(52)=dotproduct(ninjaE3,spvak2e5)
      acd38(53)=abb38(26)
      acd38(54)=-acd38(19)*acd38(18)
      acd38(55)=-acd38(17)*acd38(16)
      acd38(56)=-acd38(15)*acd38(14)
      acd38(57)=-acd38(7)*acd38(6)
      acd38(54)=acd38(57)+acd38(56)+acd38(54)+acd38(55)
      acd38(54)=acd38(8)*acd38(54)
      acd38(55)=acd38(13)*acd38(12)
      acd38(56)=acd38(10)*acd38(9)
      acd38(55)=acd38(55)+acd38(56)
      acd38(55)=acd38(11)*acd38(55)
      acd38(56)=acd38(5)*acd38(4)
      acd38(57)=acd38(2)*acd38(1)
      acd38(56)=acd38(56)+acd38(57)
      acd38(56)=acd38(3)*acd38(56)
      acd38(54)=acd38(54)+acd38(55)+acd38(56)
      acd38(55)=-acd38(19)*acd38(29)
      acd38(56)=-acd38(17)*acd38(28)
      acd38(57)=-acd38(15)*acd38(27)
      acd38(58)=-acd38(7)*acd38(24)
      acd38(55)=acd38(58)+acd38(57)+acd38(55)+acd38(56)
      acd38(55)=acd38(8)*acd38(55)
      acd38(56)=acd38(13)*acd38(26)
      acd38(57)=acd38(10)*acd38(25)
      acd38(56)=acd38(56)+acd38(57)
      acd38(56)=acd38(11)*acd38(56)
      acd38(57)=acd38(5)*acd38(23)
      acd38(58)=acd38(2)*acd38(22)
      acd38(57)=acd38(57)+acd38(58)
      acd38(57)=acd38(3)*acd38(57)
      acd38(58)=acd38(52)*acd38(53)
      acd38(59)=acd38(50)*acd38(51)
      acd38(60)=acd38(48)*acd38(49)
      acd38(61)=acd38(46)*acd38(47)
      acd38(62)=acd38(44)*acd38(45)
      acd38(63)=acd38(42)*acd38(43)
      acd38(64)=acd38(40)*acd38(41)
      acd38(65)=acd38(35)*acd38(36)
      acd38(66)=acd38(20)*acd38(21)
      acd38(67)=acd38(19)*acd38(38)
      acd38(68)=acd38(17)*acd38(39)
      acd38(69)=acd38(15)*acd38(32)
      acd38(70)=acd38(13)*acd38(33)
      acd38(71)=acd38(10)*acd38(34)
      acd38(72)=acd38(7)*acd38(37)
      acd38(73)=acd38(5)*acd38(30)
      acd38(74)=acd38(2)*acd38(31)
      acd38(55)=acd38(55)+acd38(57)+acd38(56)+acd38(74)+acd38(73)+acd38(72)+acd&
      &38(71)+acd38(70)+acd38(69)+acd38(68)+acd38(67)+2.0_ki*acd38(66)+acd38(65&
      &)+acd38(64)+acd38(63)+acd38(62)+acd38(61)+acd38(60)+acd38(58)+acd38(59)
      brack(ninjaidxt0x0mu0)=acd38(55)
      brack(ninjaidxt0x0mu2)=0.0_ki
      brack(ninjaidxt0x1mu0)=acd38(54)
      brack(ninjaidxt0x2mu0)=0.0_ki
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d38h0_ninja_t2")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1, only: epspow
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_abbrevd38h0
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a0, a1, b, c
      complex(ki), dimension(0:2), intent(in) :: param
      complex(ki), dimension(4) :: vecA0, vecA1, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      qshift = -k2
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
end module     p2_part21part21_part25part25part21_d38h0l132
