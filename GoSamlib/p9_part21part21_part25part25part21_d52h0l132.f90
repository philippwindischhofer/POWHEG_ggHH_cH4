module     p9_part21part21_part25part25part21_d52h0l132
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d52h0l132.f90
   ! generator: buildfortran_tn3.py
   use p9_part21part21_part25part25part21_config, only: ki
   use p9_part21part21_part25part25part21_util, only: cond_t, d => metric_tensor
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
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd52h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(9) :: acd52
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd52(1)=dotproduct(ninjaE3,spvae1e2)
      acd52(2)=dotproduct(ninjaE3,spvae5e1)
      acd52(3)=dotproduct(ninjaE3,spvae2e5)
      acd52(4)=abb52(12)
      acd52(5)=dotproduct(ninjaE3,spvae2e1)
      acd52(6)=dotproduct(ninjaE3,spvae1e5)
      acd52(7)=dotproduct(ninjaE3,spvae5e2)
      acd52(8)=acd52(3)*acd52(2)*acd52(1)
      acd52(9)=acd52(7)*acd52(6)*acd52(5)
      acd52(8)=acd52(8)+acd52(9)
      acd52(8)=acd52(4)*acd52(8)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd52(8)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd52h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(57) :: acd52
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd52(1)=dotproduct(ninjaA1,spvae1e2)
      acd52(2)=dotproduct(ninjaE3,spvae2e5)
      acd52(3)=dotproduct(ninjaE3,spvae5e1)
      acd52(4)=abb52(12)
      acd52(5)=dotproduct(ninjaA1,spvae2e5)
      acd52(6)=dotproduct(ninjaE3,spvae1e2)
      acd52(7)=dotproduct(ninjaA1,spvae2e1)
      acd52(8)=dotproduct(ninjaE3,spvae5e2)
      acd52(9)=dotproduct(ninjaE3,spvae1e5)
      acd52(10)=dotproduct(ninjaA1,spvae5e2)
      acd52(11)=dotproduct(ninjaE3,spvae2e1)
      acd52(12)=dotproduct(ninjaA1,spvae5e1)
      acd52(13)=dotproduct(ninjaA1,spvae1e5)
      acd52(14)=dotproduct(ninjaA0,ninjaE3)
      acd52(15)=abb52(16)
      acd52(16)=abb52(24)
      acd52(17)=abb52(26)
      acd52(18)=abb52(30)
      acd52(19)=dotproduct(ninjaA0,spvae1e2)
      acd52(20)=dotproduct(ninjaA0,spvae2e5)
      acd52(21)=dotproduct(ninjaA0,spvae2e1)
      acd52(22)=dotproduct(ninjaA0,spvae5e2)
      acd52(23)=dotproduct(ninjaA0,spvae5e1)
      acd52(24)=dotproduct(ninjaA0,spvae1e5)
      acd52(25)=abb52(2)
      acd52(26)=abb52(18)
      acd52(27)=dotproduct(ninjaE3,spvae2k1)
      acd52(28)=abb52(13)
      acd52(29)=dotproduct(ninjaE3,spvae5k2)
      acd52(30)=abb52(37)
      acd52(31)=dotproduct(ninjaE3,spvae2k5)
      acd52(32)=abb52(33)
      acd52(33)=abb52(5)
      acd52(34)=dotproduct(ninjaE3,spvak2e5)
      acd52(35)=abb52(15)
      acd52(36)=abb52(31)
      acd52(37)=dotproduct(ninjaE3,spvak1e2)
      acd52(38)=abb52(19)
      acd52(39)=dotproduct(ninjaE3,spvak5e2)
      acd52(40)=abb52(35)
      acd52(41)=abb52(7)
      acd52(42)=abb52(23)
      acd52(43)=abb52(21)
      acd52(44)=abb52(28)
      acd52(45)=abb52(17)
      acd52(46)=abb52(27)
      acd52(47)=acd52(8)*acd52(13)
      acd52(48)=acd52(9)*acd52(10)
      acd52(47)=acd52(47)+acd52(48)
      acd52(47)=acd52(11)*acd52(47)
      acd52(48)=acd52(2)*acd52(12)
      acd52(49)=acd52(3)*acd52(5)
      acd52(48)=acd52(48)+acd52(49)
      acd52(48)=acd52(6)*acd52(48)
      acd52(49)=acd52(9)*acd52(8)
      acd52(50)=acd52(7)*acd52(49)
      acd52(51)=acd52(3)*acd52(2)
      acd52(52)=acd52(1)*acd52(51)
      acd52(47)=acd52(48)+acd52(47)+acd52(50)+acd52(52)
      acd52(47)=acd52(4)*acd52(47)
      acd52(48)=acd52(39)*acd52(40)
      acd52(50)=acd52(37)*acd52(38)
      acd52(52)=acd52(34)*acd52(35)
      acd52(53)=2.0_ki*acd52(14)
      acd52(54)=acd52(16)*acd52(53)
      acd52(55)=acd52(8)*acd52(33)
      acd52(56)=acd52(9)*acd52(36)
      acd52(48)=acd52(56)+acd52(55)+acd52(54)+acd52(52)+acd52(48)+acd52(50)
      acd52(48)=acd52(11)*acd52(48)
      acd52(50)=acd52(31)*acd52(32)
      acd52(52)=acd52(29)*acd52(30)
      acd52(54)=acd52(27)*acd52(28)
      acd52(55)=acd52(15)*acd52(53)
      acd52(56)=acd52(2)*acd52(25)
      acd52(57)=acd52(3)*acd52(26)
      acd52(50)=acd52(57)+acd52(56)+acd52(55)+acd52(54)+acd52(50)+acd52(52)
      acd52(50)=acd52(6)*acd52(50)
      acd52(52)=acd52(31)*acd52(46)
      acd52(54)=acd52(29)*acd52(45)
      acd52(55)=acd52(27)*acd52(44)
      acd52(56)=acd52(18)*acd52(53)
      acd52(52)=acd52(56)+acd52(55)+acd52(52)+acd52(54)
      acd52(52)=acd52(9)*acd52(52)
      acd52(54)=acd52(39)*acd52(43)
      acd52(55)=acd52(37)*acd52(42)
      acd52(56)=acd52(34)*acd52(41)
      acd52(53)=acd52(17)*acd52(53)
      acd52(53)=acd52(53)+acd52(56)+acd52(54)+acd52(55)
      acd52(53)=acd52(3)*acd52(53)
      acd52(54)=acd52(8)*acd52(24)
      acd52(55)=acd52(9)*acd52(22)
      acd52(54)=acd52(54)+acd52(55)
      acd52(54)=acd52(11)*acd52(54)
      acd52(55)=acd52(2)*acd52(23)
      acd52(56)=acd52(3)*acd52(20)
      acd52(55)=acd52(55)+acd52(56)
      acd52(55)=acd52(6)*acd52(55)
      acd52(49)=acd52(21)*acd52(49)
      acd52(51)=acd52(19)*acd52(51)
      acd52(49)=acd52(55)+acd52(54)+acd52(49)+acd52(51)
      acd52(49)=acd52(4)*acd52(49)
      acd52(48)=acd52(49)+acd52(50)+acd52(48)+acd52(52)+acd52(53)
      brack(ninjaidxt0x0mu0)=acd52(48)
      brack(ninjaidxt0x0mu2)=0.0_ki
      brack(ninjaidxt0x1mu0)=acd52(47)
      brack(ninjaidxt0x2mu0)=0.0_ki
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p9_part21part21_part25part25part21_d52h0_ninja_t2")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki => ki_nin
      use p9_part21part21_part25part25part21_globalsl1, only: epspow
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_abbrevd52h0
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
end module     p9_part21part21_part25part25part21_d52h0l132
