module     p2_part21part21_part25part25part21_d204h0l132
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d204h0l132.f90
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
      use p2_part21part21_part25part25part21_abbrevd204h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(1) :: acd204
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
      use p2_part21part21_part25part25part21_abbrevd204h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(56) :: acd204
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd204(1)=dotproduct(ninjaA0,ninjaE3)
      acd204(2)=dotproduct(ninjaE3,spvae1e5)
      acd204(3)=abb204(104)
      acd204(4)=dotproduct(ninjaE3,spvae5e1)
      acd204(5)=abb204(98)
      acd204(6)=dotproduct(ninjaE3,spvae5k2)
      acd204(7)=abb204(100)
      acd204(8)=dotproduct(ninjaE3,spval4e1)
      acd204(9)=abb204(55)
      acd204(10)=dotproduct(ninjaE3,spvak2e1)
      acd204(11)=abb204(24)
      acd204(12)=dotproduct(ninjaE3,spvae2k1)
      acd204(13)=abb204(74)
      acd204(14)=dotproduct(ninjaE3,spvae5l4)
      acd204(15)=abb204(59)
      acd204(16)=dotproduct(ninjaE3,spvae2k5)
      acd204(17)=abb204(109)
      acd204(18)=dotproduct(ninjaE3,spvak2e5)
      acd204(19)=abb204(84)
      acd204(20)=dotproduct(ninjaE3,spvae1l4)
      acd204(21)=abb204(47)
      acd204(22)=dotproduct(ninjaE3,spvak1e2)
      acd204(23)=abb204(58)
      acd204(24)=dotproduct(ninjaE3,spvae1k2)
      acd204(25)=abb204(70)
      acd204(26)=dotproduct(ninjaE3,spvak5e2)
      acd204(27)=abb204(111)
      acd204(28)=dotproduct(ninjaE3,spval4e5)
      acd204(29)=abb204(91)
      acd204(30)=abb204(17)
      acd204(31)=abb204(44)
      acd204(32)=abb204(87)
      acd204(33)=abb204(63)
      acd204(34)=abb204(79)
      acd204(35)=abb204(106)
      acd204(36)=abb204(101)
      acd204(37)=dotproduct(ninjaE3,spvae2e1)
      acd204(38)=dotproduct(ninjaE3,spvae5e2)
      acd204(39)=abb204(62)
      acd204(40)=abb204(35)
      acd204(41)=abb204(40)
      acd204(42)=abb204(86)
      acd204(43)=abb204(102)
      acd204(44)=dotproduct(ninjaE3,spvae1e2)
      acd204(45)=dotproduct(ninjaE3,spvae2e5)
      acd204(46)=abb204(88)
      acd204(47)=abb204(99)
      acd204(48)=-acd204(18)*acd204(19)
      acd204(49)=-acd204(28)*acd204(29)
      acd204(50)=acd204(26)*acd204(27)
      acd204(51)=acd204(22)*acd204(23)
      acd204(52)=acd204(24)*acd204(25)
      acd204(53)=acd204(20)*acd204(21)
      acd204(54)=2.0_ki*acd204(1)
      acd204(55)=acd204(5)*acd204(54)
      acd204(48)=acd204(55)+acd204(53)+acd204(52)+acd204(51)+acd204(50)+acd204(&
      &48)+acd204(49)
      acd204(48)=acd204(48)*acd204(54)
      acd204(49)=-acd204(18)*acd204(40)
      acd204(50)=-acd204(28)*acd204(32)
      acd204(51)=acd204(26)*acd204(42)
      acd204(52)=acd204(22)*acd204(41)
      acd204(49)=acd204(52)+acd204(51)+acd204(49)+acd204(50)
      acd204(49)=acd204(20)*acd204(49)
      acd204(50)=acd204(28)*acd204(47)
      acd204(51)=acd204(26)*acd204(46)
      acd204(52)=acd204(22)*acd204(43)
      acd204(50)=acd204(52)+acd204(50)+acd204(51)
      acd204(50)=acd204(24)*acd204(50)
      acd204(51)=-acd204(39)*acd204(44)*acd204(45)
      acd204(48)=acd204(48)+acd204(49)+acd204(51)+acd204(50)
      acd204(48)=acd204(4)*acd204(48)
      acd204(49)=-acd204(6)*acd204(7)
      acd204(50)=acd204(16)*acd204(17)
      acd204(51)=-acd204(14)*acd204(15)
      acd204(52)=acd204(12)*acd204(13)
      acd204(53)=acd204(10)*acd204(11)
      acd204(55)=acd204(8)*acd204(9)
      acd204(56)=acd204(3)*acd204(54)
      acd204(49)=acd204(56)+acd204(55)+acd204(53)+acd204(52)+acd204(51)+acd204(&
      &49)+acd204(50)
      acd204(49)=acd204(49)*acd204(54)
      acd204(50)=-acd204(6)*acd204(30)
      acd204(51)=acd204(16)*acd204(33)
      acd204(52)=-acd204(14)*acd204(32)
      acd204(53)=acd204(12)*acd204(31)
      acd204(50)=acd204(53)+acd204(52)+acd204(50)+acd204(51)
      acd204(50)=acd204(8)*acd204(50)
      acd204(51)=acd204(16)*acd204(36)
      acd204(52)=acd204(14)*acd204(35)
      acd204(53)=acd204(12)*acd204(34)
      acd204(51)=acd204(53)+acd204(51)+acd204(52)
      acd204(51)=acd204(10)*acd204(51)
      acd204(52)=-acd204(39)*acd204(37)*acd204(38)
      acd204(49)=acd204(49)+acd204(50)+acd204(52)+acd204(51)
      acd204(49)=acd204(2)*acd204(49)
      acd204(48)=acd204(48)+acd204(49)
      brack(ninjaidxt0x0mu0)=acd204(48)
      brack(ninjaidxt0x0mu2)=0.0_ki
      brack(ninjaidxt0x1mu0)=0.0_ki
      brack(ninjaidxt0x2mu0)=0.0_ki
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d204h0_ninja_t2")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1, only: epspow
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_abbrevd204h0
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a0, a1, b, c
      complex(ki), dimension(0:2), intent(in) :: param
      complex(ki), dimension(4) :: vecA0, vecA1, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      qshift = k2-k4
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
end module     p2_part21part21_part25part25part21_d204h0l132
