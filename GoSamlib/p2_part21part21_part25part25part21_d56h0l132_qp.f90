module     p2_part21part21_part25part25part21_d56h0l132_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d56h0l132_qp.f90
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
      use p2_part21part21_part25part25part21_abbrevd56h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(9) :: acd56
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd56(1)=dotproduct(ninjaE3,spvae1e2)
      acd56(2)=dotproduct(ninjaE3,spvae5e1)
      acd56(3)=dotproduct(ninjaE3,spvae2e5)
      acd56(4)=abb56(7)
      acd56(5)=dotproduct(ninjaE3,spvae2e1)
      acd56(6)=dotproduct(ninjaE3,spvae1e5)
      acd56(7)=dotproduct(ninjaE3,spvae5e2)
      acd56(8)=acd56(3)*acd56(2)*acd56(1)
      acd56(9)=acd56(7)*acd56(6)*acd56(5)
      acd56(8)=acd56(8)+acd56(9)
      acd56(8)=acd56(4)*acd56(8)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd56(8)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd56h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(57) :: acd56
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd56(1)=dotproduct(ninjaA1,spvae5e1)
      acd56(2)=dotproduct(ninjaE3,spvae2e5)
      acd56(3)=dotproduct(ninjaE3,spvae1e2)
      acd56(4)=abb56(7)
      acd56(5)=dotproduct(ninjaA1,spvae2e5)
      acd56(6)=dotproduct(ninjaE3,spvae5e1)
      acd56(7)=dotproduct(ninjaA1,spvae1e5)
      acd56(8)=dotproduct(ninjaE3,spvae5e2)
      acd56(9)=dotproduct(ninjaE3,spvae2e1)
      acd56(10)=dotproduct(ninjaA1,spvae5e2)
      acd56(11)=dotproduct(ninjaE3,spvae1e5)
      acd56(12)=dotproduct(ninjaA1,spvae1e2)
      acd56(13)=dotproduct(ninjaA1,spvae2e1)
      acd56(14)=dotproduct(ninjaA0,ninjaE3)
      acd56(15)=abb56(30)
      acd56(16)=abb56(27)
      acd56(17)=abb56(25)
      acd56(18)=abb56(9)
      acd56(19)=dotproduct(ninjaA0,spvae5e1)
      acd56(20)=dotproduct(ninjaA0,spvae2e5)
      acd56(21)=dotproduct(ninjaA0,spvae1e5)
      acd56(22)=dotproduct(ninjaA0,spvae5e2)
      acd56(23)=dotproduct(ninjaA0,spvae1e2)
      acd56(24)=dotproduct(ninjaA0,spvae2e1)
      acd56(25)=abb56(4)
      acd56(26)=abb56(15)
      acd56(27)=dotproduct(ninjaE3,spvae2k1)
      acd56(28)=abb56(29)
      acd56(29)=dotproduct(ninjaE3,spvae1k2)
      acd56(30)=abb56(34)
      acd56(31)=dotproduct(ninjaE3,spvae1k5)
      acd56(32)=abb56(31)
      acd56(33)=dotproduct(ninjaE3,spvak1e2)
      acd56(34)=abb56(16)
      acd56(35)=dotproduct(ninjaE3,spvak2e1)
      acd56(36)=abb56(19)
      acd56(37)=dotproduct(ninjaE3,spvak5e1)
      acd56(38)=abb56(26)
      acd56(39)=abb56(5)
      acd56(40)=abb56(6)
      acd56(41)=abb56(18)
      acd56(42)=abb56(35)
      acd56(43)=abb56(32)
      acd56(44)=abb56(14)
      acd56(45)=abb56(23)
      acd56(46)=abb56(28)
      acd56(47)=acd56(9)*acd56(10)
      acd56(48)=acd56(8)*acd56(13)
      acd56(47)=acd56(47)+acd56(48)
      acd56(47)=acd56(11)*acd56(47)
      acd56(48)=acd56(3)*acd56(5)
      acd56(49)=acd56(2)*acd56(12)
      acd56(48)=acd56(48)+acd56(49)
      acd56(48)=acd56(6)*acd56(48)
      acd56(49)=acd56(8)*acd56(9)
      acd56(50)=acd56(7)*acd56(49)
      acd56(51)=acd56(2)*acd56(3)
      acd56(52)=acd56(1)*acd56(51)
      acd56(47)=acd56(48)+acd56(47)+acd56(50)+acd56(52)
      acd56(47)=acd56(4)*acd56(47)
      acd56(48)=acd56(37)*acd56(43)
      acd56(50)=acd56(35)*acd56(42)
      acd56(52)=acd56(33)*acd56(39)
      acd56(53)=2.0_ki*acd56(14)
      acd56(54)=acd56(17)*acd56(53)
      acd56(55)=acd56(9)*acd56(41)
      acd56(56)=acd56(8)*acd56(40)
      acd56(48)=acd56(56)+acd56(55)+acd56(54)+acd56(52)+acd56(48)+acd56(50)
      acd56(48)=acd56(11)*acd56(48)
      acd56(50)=acd56(31)*acd56(32)
      acd56(52)=acd56(29)*acd56(30)
      acd56(54)=acd56(27)*acd56(28)
      acd56(55)=acd56(15)*acd56(53)
      acd56(56)=acd56(3)*acd56(26)
      acd56(57)=acd56(2)*acd56(25)
      acd56(50)=acd56(57)+acd56(56)+acd56(55)+acd56(54)+acd56(50)+acd56(52)
      acd56(50)=acd56(6)*acd56(50)
      acd56(52)=acd56(31)*acd56(46)
      acd56(54)=acd56(29)*acd56(45)
      acd56(55)=acd56(27)*acd56(44)
      acd56(56)=acd56(18)*acd56(53)
      acd56(52)=acd56(56)+acd56(55)+acd56(52)+acd56(54)
      acd56(52)=acd56(8)*acd56(52)
      acd56(54)=acd56(37)*acd56(38)
      acd56(55)=acd56(35)*acd56(36)
      acd56(56)=acd56(33)*acd56(34)
      acd56(53)=acd56(16)*acd56(53)
      acd56(53)=acd56(53)+acd56(56)+acd56(54)+acd56(55)
      acd56(53)=acd56(2)*acd56(53)
      acd56(54)=acd56(9)*acd56(22)
      acd56(55)=acd56(8)*acd56(24)
      acd56(54)=acd56(54)+acd56(55)
      acd56(54)=acd56(11)*acd56(54)
      acd56(55)=acd56(3)*acd56(20)
      acd56(56)=acd56(2)*acd56(23)
      acd56(55)=acd56(55)+acd56(56)
      acd56(55)=acd56(6)*acd56(55)
      acd56(49)=acd56(21)*acd56(49)
      acd56(51)=acd56(19)*acd56(51)
      acd56(49)=acd56(55)+acd56(54)+acd56(49)+acd56(51)
      acd56(49)=acd56(4)*acd56(49)
      acd56(48)=acd56(49)+acd56(50)+acd56(48)+acd56(52)+acd56(53)
      brack(ninjaidxt0x0mu0)=acd56(48)
      brack(ninjaidxt0x0mu2)=0.0_ki
      brack(ninjaidxt0x1mu0)=acd56(47)
      brack(ninjaidxt0x2mu0)=0.0_ki
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d56h0_qp_ninja_t2")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd56h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a0, a1, b, c
      complex(ki), dimension(0:2), intent(in) :: param
      complex(ki), dimension(4) :: vecA0, vecA1, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      qshift = -k2+k5
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
end module     p2_part21part21_part25part25part21_d56h0l132_qp
