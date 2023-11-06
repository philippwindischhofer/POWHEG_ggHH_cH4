module     p9_part21part21_part25part25part21_d198h0l132
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d198h0l132.f90
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
      use p9_part21part21_part25part25part21_abbrevd198h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(1) :: acd198
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
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd198h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(57) :: acd198
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd198(1)=dotproduct(ninjaA0,ninjaE3)
      acd198(2)=dotproduct(ninjaE3,spvae5e2)
      acd198(3)=abb198(101)
      acd198(4)=dotproduct(ninjaE3,spvae2e5)
      acd198(5)=abb198(103)
      acd198(6)=dotproduct(ninjaE3,spvae1k1)
      acd198(7)=abb198(91)
      acd198(8)=dotproduct(ninjaE3,spvak2e1)
      acd198(9)=abb198(25)
      acd198(10)=abb198(33)
      acd198(11)=dotproduct(ninjaE3,spvae1k2)
      acd198(12)=abb198(43)
      acd198(13)=dotproduct(ninjaE3,spvae2k1)
      acd198(14)=abb198(52)
      acd198(15)=dotproduct(ninjaE3,spval3e1)
      acd198(16)=abb198(123)
      acd198(17)=dotproduct(ninjaE3,spvae1l4)
      acd198(18)=abb198(113)
      acd198(19)=abb198(28)
      acd198(20)=dotproduct(ninjaE3,spvak1e2)
      acd198(21)=abb198(41)
      acd198(22)=dotproduct(ninjaE3,spvae1l3)
      acd198(23)=abb198(116)
      acd198(24)=dotproduct(ninjaE3,spvak1e1)
      acd198(25)=abb198(68)
      acd198(26)=dotproduct(ninjaE3,spval4e1)
      acd198(27)=abb198(114)
      acd198(28)=abb198(9)
      acd198(29)=abb198(92)
      acd198(30)=abb198(15)
      acd198(31)=abb198(78)
      acd198(32)=abb198(79)
      acd198(33)=abb198(14)
      acd198(34)=abb198(56)
      acd198(35)=abb198(75)
      acd198(36)=abb198(87)
      acd198(37)=abb198(128)
      acd198(38)=dotproduct(ninjaE3,spvae2e1)
      acd198(39)=dotproduct(ninjaE3,spvae1e5)
      acd198(40)=abb198(108)
      acd198(41)=abb198(72)
      acd198(42)=abb198(57)
      acd198(43)=abb198(84)
      acd198(44)=abb198(89)
      acd198(45)=abb198(77)
      acd198(46)=abb198(121)
      acd198(47)=dotproduct(ninjaE3,spvae1e2)
      acd198(48)=dotproduct(ninjaE3,spvae5e1)
      acd198(49)=acd198(26)*acd198(27)
      acd198(50)=acd198(24)*acd198(25)
      acd198(51)=acd198(20)*acd198(21)
      acd198(52)=acd198(22)*acd198(23)
      acd198(53)=acd198(11)*acd198(19)
      acd198(54)=acd198(8)*acd198(10)
      acd198(55)=2.0_ki*acd198(1)
      acd198(56)=acd198(5)*acd198(55)
      acd198(49)=acd198(56)+acd198(54)+acd198(53)+acd198(52)+acd198(51)+acd198(&
      &49)+acd198(50)
      acd198(49)=acd198(49)*acd198(55)
      acd198(50)=acd198(26)*acd198(46)
      acd198(51)=acd198(24)*acd198(45)
      acd198(52)=acd198(20)*acd198(44)
      acd198(50)=acd198(52)+acd198(50)+acd198(51)
      acd198(50)=acd198(22)*acd198(50)
      acd198(51)=acd198(26)*acd198(43)
      acd198(52)=acd198(24)*acd198(42)
      acd198(53)=acd198(20)*acd198(41)
      acd198(51)=acd198(53)+acd198(51)+acd198(52)
      acd198(51)=acd198(11)*acd198(51)
      acd198(52)=acd198(22)*acd198(34)
      acd198(53)=acd198(11)*acd198(33)
      acd198(52)=acd198(52)+acd198(53)
      acd198(52)=acd198(8)*acd198(52)
      acd198(53)=acd198(40)*acd198(47)*acd198(48)
      acd198(49)=acd198(49)+acd198(52)+acd198(51)+acd198(53)+acd198(50)
      acd198(49)=acd198(4)*acd198(49)
      acd198(50)=acd198(17)*acd198(18)
      acd198(51)=acd198(13)*acd198(14)
      acd198(52)=acd198(6)*acd198(7)
      acd198(53)=acd198(15)*acd198(16)
      acd198(54)=acd198(11)*acd198(12)
      acd198(56)=acd198(8)*acd198(9)
      acd198(57)=acd198(3)*acd198(55)
      acd198(50)=acd198(57)+acd198(56)+acd198(54)+acd198(53)+acd198(52)+acd198(&
      &50)+acd198(51)
      acd198(50)=acd198(50)*acd198(55)
      acd198(51)=acd198(17)*acd198(37)
      acd198(52)=acd198(13)*acd198(36)
      acd198(53)=acd198(6)*acd198(29)
      acd198(54)=acd198(11)*acd198(35)
      acd198(51)=acd198(54)+acd198(53)+acd198(51)+acd198(52)
      acd198(51)=acd198(15)*acd198(51)
      acd198(52)=acd198(17)*acd198(32)
      acd198(53)=acd198(13)*acd198(31)
      acd198(54)=acd198(6)*acd198(28)
      acd198(55)=acd198(11)*acd198(30)
      acd198(52)=acd198(55)+acd198(54)+acd198(52)+acd198(53)
      acd198(52)=acd198(8)*acd198(52)
      acd198(53)=acd198(40)*acd198(38)*acd198(39)
      acd198(50)=acd198(50)+acd198(52)+acd198(53)+acd198(51)
      acd198(50)=acd198(2)*acd198(50)
      acd198(49)=acd198(49)+acd198(50)
      brack(ninjaidxt0x0mu0)=acd198(49)
      brack(ninjaidxt0x0mu2)=0.0_ki
      brack(ninjaidxt0x1mu0)=0.0_ki
      brack(ninjaidxt0x2mu0)=0.0_ki
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p9_part21part21_part25part25part21_d198h0_ninja_t2")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki => ki_nin
      use p9_part21part21_part25part25part21_globalsl1, only: epspow
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_abbrevd198h0
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a0, a1, b, c
      complex(ki), dimension(0:2), intent(in) :: param
      complex(ki), dimension(4) :: vecA0, vecA1, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      qshift = -k3+k2-k5
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
end module     p9_part21part21_part25part25part21_d198h0l132
