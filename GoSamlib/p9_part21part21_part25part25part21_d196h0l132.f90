module     p9_part21part21_part25part25part21_d196h0l132
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d196h0l132.f90
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
      use p9_part21part21_part25part25part21_abbrevd196h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(1) :: acd196
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
      use p9_part21part21_part25part25part21_abbrevd196h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(57) :: acd196
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd196(1)=dotproduct(ninjaA0,ninjaE3)
      acd196(2)=dotproduct(ninjaE3,spvae2e5)
      acd196(3)=abb196(95)
      acd196(4)=dotproduct(ninjaE3,spvae5e2)
      acd196(5)=abb196(89)
      acd196(6)=dotproduct(ninjaE3,spvae5k1)
      acd196(7)=abb196(10)
      acd196(8)=dotproduct(ninjaE3,spvak2e1)
      acd196(9)=abb196(32)
      acd196(10)=dotproduct(ninjaE3,spvae1k2)
      acd196(11)=abb196(39)
      acd196(12)=dotproduct(ninjaE3,spval3e1)
      acd196(13)=abb196(72)
      acd196(14)=dotproduct(ninjaE3,spvae1l4)
      acd196(15)=abb196(110)
      acd196(16)=dotproduct(ninjaE3,spvae1k1)
      acd196(17)=abb196(86)
      acd196(18)=abb196(75)
      acd196(19)=dotproduct(ninjaE3,spvae1l3)
      acd196(20)=abb196(93)
      acd196(21)=abb196(23)
      acd196(22)=dotproduct(ninjaE3,spvak1e5)
      acd196(23)=abb196(26)
      acd196(24)=dotproduct(ninjaE3,spvak1e1)
      acd196(25)=abb196(88)
      acd196(26)=dotproduct(ninjaE3,spval4e1)
      acd196(27)=abb196(116)
      acd196(28)=abb196(84)
      acd196(29)=abb196(82)
      acd196(30)=abb196(21)
      acd196(31)=abb196(49)
      acd196(32)=abb196(56)
      acd196(33)=abb196(38)
      acd196(34)=abb196(119)
      acd196(35)=abb196(55)
      acd196(36)=dotproduct(ninjaE3,spvae1e2)
      acd196(37)=dotproduct(ninjaE3,spvae5e1)
      acd196(38)=abb196(99)
      acd196(39)=abb196(15)
      acd196(40)=abb196(20)
      acd196(41)=abb196(91)
      acd196(42)=abb196(65)
      acd196(43)=abb196(117)
      acd196(44)=abb196(101)
      acd196(45)=abb196(76)
      acd196(46)=abb196(73)
      acd196(47)=dotproduct(ninjaE3,spvae2e1)
      acd196(48)=dotproduct(ninjaE3,spvae1e5)
      acd196(49)=acd196(26)*acd196(27)
      acd196(50)=acd196(24)*acd196(25)
      acd196(51)=acd196(22)*acd196(23)
      acd196(52)=acd196(19)*acd196(20)
      acd196(53)=acd196(10)*acd196(21)
      acd196(54)=acd196(8)*acd196(18)
      acd196(55)=2.0_ki*acd196(1)
      acd196(56)=acd196(5)*acd196(55)
      acd196(49)=acd196(56)+acd196(54)+acd196(53)+acd196(52)+acd196(51)+acd196(&
      &49)+acd196(50)
      acd196(49)=acd196(49)*acd196(55)
      acd196(50)=acd196(26)*acd196(43)
      acd196(51)=acd196(24)*acd196(42)
      acd196(52)=acd196(22)*acd196(41)
      acd196(50)=acd196(52)+acd196(50)+acd196(51)
      acd196(50)=acd196(19)*acd196(50)
      acd196(51)=acd196(26)*acd196(46)
      acd196(52)=acd196(24)*acd196(45)
      acd196(53)=acd196(22)*acd196(44)
      acd196(51)=acd196(53)+acd196(51)+acd196(52)
      acd196(51)=acd196(10)*acd196(51)
      acd196(52)=acd196(19)*acd196(39)
      acd196(53)=acd196(10)*acd196(40)
      acd196(52)=acd196(52)+acd196(53)
      acd196(52)=acd196(8)*acd196(52)
      acd196(53)=-acd196(38)*acd196(47)*acd196(48)
      acd196(49)=acd196(49)+acd196(52)+acd196(51)+acd196(53)+acd196(50)
      acd196(49)=acd196(4)*acd196(49)
      acd196(50)=acd196(16)*acd196(17)
      acd196(51)=acd196(14)*acd196(15)
      acd196(52)=acd196(6)*acd196(7)
      acd196(53)=acd196(12)*acd196(13)
      acd196(54)=acd196(10)*acd196(11)
      acd196(56)=acd196(8)*acd196(9)
      acd196(57)=acd196(3)*acd196(55)
      acd196(50)=acd196(57)+acd196(56)+acd196(54)+acd196(53)+acd196(52)+acd196(&
      &50)+acd196(51)
      acd196(50)=acd196(50)*acd196(55)
      acd196(51)=acd196(16)*acd196(35)
      acd196(52)=acd196(14)*acd196(34)
      acd196(53)=acd196(6)*acd196(29)
      acd196(54)=acd196(10)*acd196(33)
      acd196(51)=acd196(54)+acd196(53)+acd196(51)+acd196(52)
      acd196(51)=acd196(12)*acd196(51)
      acd196(52)=acd196(16)*acd196(32)
      acd196(53)=acd196(14)*acd196(31)
      acd196(54)=acd196(6)*acd196(28)
      acd196(55)=-acd196(10)*acd196(30)
      acd196(52)=acd196(55)+acd196(54)+acd196(52)+acd196(53)
      acd196(52)=acd196(8)*acd196(52)
      acd196(53)=-acd196(38)*acd196(36)*acd196(37)
      acd196(50)=acd196(50)+acd196(52)+acd196(53)+acd196(51)
      acd196(50)=acd196(2)*acd196(50)
      acd196(49)=acd196(49)+acd196(50)
      brack(ninjaidxt0x0mu0)=acd196(49)
      brack(ninjaidxt0x0mu2)=0.0_ki
      brack(ninjaidxt0x1mu0)=0.0_ki
      brack(ninjaidxt0x2mu0)=0.0_ki
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p9_part21part21_part25part25part21_d196h0_ninja_t2")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki => ki_nin
      use p9_part21part21_part25part25part21_globalsl1, only: epspow
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_abbrevd196h0
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a0, a1, b, c
      complex(ki), dimension(0:2), intent(in) :: param
      complex(ki), dimension(4) :: vecA0, vecA1, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      qshift = -k3
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
end module     p9_part21part21_part25part25part21_d196h0l132
