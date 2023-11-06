module     p2_part21part21_part25part25part21_d145h0l132
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d145h0l132.f90
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
      use p2_part21part21_part25part25part21_abbrevd145h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(9) :: acd145
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd145(1)=dotproduct(ninjaE3,spvae1e2)
      acd145(2)=dotproduct(ninjaE3,spvae5e1)
      acd145(3)=dotproduct(ninjaE3,spvae2e5)
      acd145(4)=abb145(7)
      acd145(5)=dotproduct(ninjaE3,spvae2e1)
      acd145(6)=dotproduct(ninjaE3,spvae1e5)
      acd145(7)=dotproduct(ninjaE3,spvae5e2)
      acd145(8)=-acd145(3)*acd145(2)*acd145(1)
      acd145(9)=-acd145(7)*acd145(6)*acd145(5)
      acd145(8)=acd145(8)+acd145(9)
      acd145(8)=acd145(4)*acd145(8)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd145(8)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_color
      use p2_part21part21_part25part25part21_abbrevd145h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(57) :: acd145
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd145(1)=dotproduct(ninjaA1,spvae2e1)
      acd145(2)=dotproduct(ninjaE3,spvae1e5)
      acd145(3)=dotproduct(ninjaE3,spvae5e2)
      acd145(4)=abb145(7)
      acd145(5)=dotproduct(ninjaA1,spvae1e2)
      acd145(6)=dotproduct(ninjaE3,spvae5e1)
      acd145(7)=dotproduct(ninjaE3,spvae2e5)
      acd145(8)=dotproduct(ninjaA1,spvae5e1)
      acd145(9)=dotproduct(ninjaE3,spvae1e2)
      acd145(10)=dotproduct(ninjaA1,spvae2e5)
      acd145(11)=dotproduct(ninjaA1,spvae1e5)
      acd145(12)=dotproduct(ninjaE3,spvae2e1)
      acd145(13)=dotproduct(ninjaA1,spvae5e2)
      acd145(14)=dotproduct(ninjaA0,ninjaE3)
      acd145(15)=abb145(19)
      acd145(16)=abb145(25)
      acd145(17)=abb145(22)
      acd145(18)=abb145(27)
      acd145(19)=dotproduct(ninjaA0,spvae2e1)
      acd145(20)=dotproduct(ninjaA0,spvae1e2)
      acd145(21)=dotproduct(ninjaA0,spvae5e1)
      acd145(22)=dotproduct(ninjaA0,spvae2e5)
      acd145(23)=dotproduct(ninjaA0,spvae1e5)
      acd145(24)=dotproduct(ninjaA0,spvae5e2)
      acd145(25)=dotproduct(ninjaE3,spvak1e2)
      acd145(26)=abb145(4)
      acd145(27)=abb145(21)
      acd145(28)=abb145(16)
      acd145(29)=abb145(15)
      acd145(30)=dotproduct(ninjaE3,spvak5e2)
      acd145(31)=abb145(34)
      acd145(32)=dotproduct(ninjaE3,spvak2e5)
      acd145(33)=abb145(37)
      acd145(34)=abb145(6)
      acd145(35)=abb145(23)
      acd145(36)=dotproduct(ninjaE3,spvae5k2)
      acd145(37)=abb145(20)
      acd145(38)=dotproduct(ninjaE3,spvae2k1)
      acd145(39)=abb145(24)
      acd145(40)=dotproduct(ninjaE3,spvae2k5)
      acd145(41)=abb145(31)
      acd145(42)=abb145(33)
      acd145(43)=abb145(36)
      acd145(44)=abb145(8)
      acd145(45)=abb145(26)
      acd145(46)=abb145(30)
      acd145(47)=-acd145(3)*acd145(11)
      acd145(48)=-acd145(2)*acd145(13)
      acd145(47)=acd145(47)+acd145(48)
      acd145(47)=acd145(12)*acd145(47)
      acd145(48)=-acd145(7)*acd145(8)
      acd145(49)=-acd145(6)*acd145(10)
      acd145(48)=acd145(48)+acd145(49)
      acd145(48)=acd145(9)*acd145(48)
      acd145(49)=acd145(6)*acd145(7)
      acd145(50)=-acd145(5)*acd145(49)
      acd145(51)=acd145(2)*acd145(3)
      acd145(52)=-acd145(1)*acd145(51)
      acd145(47)=acd145(48)+acd145(47)+acd145(50)+acd145(52)
      acd145(47)=acd145(4)*acd145(47)
      acd145(48)=acd145(32)*acd145(33)
      acd145(50)=acd145(30)*acd145(31)
      acd145(52)=acd145(25)*acd145(26)
      acd145(53)=2.0_ki*acd145(14)
      acd145(54)=acd145(15)*acd145(53)
      acd145(55)=-acd145(3)*acd145(29)
      acd145(56)=acd145(2)*acd145(28)
      acd145(48)=acd145(56)+acd145(55)+acd145(54)+acd145(52)+acd145(48)+acd145(&
      &50)
      acd145(48)=acd145(12)*acd145(48)
      acd145(50)=acd145(40)*acd145(41)
      acd145(52)=acd145(38)*acd145(39)
      acd145(54)=acd145(36)*acd145(37)
      acd145(55)=acd145(16)*acd145(53)
      acd145(56)=-acd145(7)*acd145(35)
      acd145(57)=acd145(6)*acd145(34)
      acd145(50)=acd145(57)+acd145(56)+acd145(55)+acd145(54)+acd145(50)+acd145(&
      &52)
      acd145(50)=acd145(9)*acd145(50)
      acd145(52)=acd145(32)*acd145(43)
      acd145(54)=acd145(30)*acd145(42)
      acd145(55)=acd145(25)*acd145(27)
      acd145(56)=acd145(17)*acd145(53)
      acd145(52)=acd145(56)+acd145(55)+acd145(52)+acd145(54)
      acd145(52)=acd145(6)*acd145(52)
      acd145(54)=-acd145(40)*acd145(46)
      acd145(55)=acd145(38)*acd145(45)
      acd145(56)=acd145(36)*acd145(44)
      acd145(53)=acd145(18)*acd145(53)
      acd145(53)=acd145(53)+acd145(56)+acd145(54)+acd145(55)
      acd145(53)=acd145(2)*acd145(53)
      acd145(54)=-acd145(3)*acd145(23)
      acd145(55)=-acd145(2)*acd145(24)
      acd145(54)=acd145(54)+acd145(55)
      acd145(54)=acd145(12)*acd145(54)
      acd145(55)=-acd145(7)*acd145(21)
      acd145(56)=-acd145(6)*acd145(22)
      acd145(55)=acd145(55)+acd145(56)
      acd145(55)=acd145(9)*acd145(55)
      acd145(49)=-acd145(20)*acd145(49)
      acd145(51)=-acd145(19)*acd145(51)
      acd145(49)=acd145(55)+acd145(54)+acd145(49)+acd145(51)
      acd145(49)=acd145(4)*acd145(49)
      acd145(48)=acd145(49)+acd145(50)+acd145(48)+acd145(52)+acd145(53)
      brack(ninjaidxt0x0mu0)=acd145(48)
      brack(ninjaidxt0x0mu2)=0.0_ki
      brack(ninjaidxt0x1mu0)=acd145(47)
      brack(ninjaidxt0x2mu0)=0.0_ki
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d145h0_ninja_t2")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1, only: epspow
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_abbrevd145h0
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
end module     p2_part21part21_part25part25part21_d145h0l132
