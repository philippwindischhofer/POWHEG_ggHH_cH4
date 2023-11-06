module     p2_part21part21_part25part25part21_d149h0l132_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d149h0l132_qp.f90
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
      use p2_part21part21_part25part25part21_abbrevd149h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(9) :: acd149
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd149(1)=dotproduct(ninjaE3,spvae1e2)
      acd149(2)=dotproduct(ninjaE3,spvae5e1)
      acd149(3)=dotproduct(ninjaE3,spvae2e5)
      acd149(4)=abb149(7)
      acd149(5)=dotproduct(ninjaE3,spvae2e1)
      acd149(6)=dotproduct(ninjaE3,spvae1e5)
      acd149(7)=dotproduct(ninjaE3,spvae5e2)
      acd149(8)=-acd149(3)*acd149(2)*acd149(1)
      acd149(9)=-acd149(7)*acd149(6)*acd149(5)
      acd149(8)=acd149(8)+acd149(9)
      acd149(8)=acd149(4)*acd149(8)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd149(8)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd149h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(57) :: acd149
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd149(1)=dotproduct(ninjaA1,spvae1e5)
      acd149(2)=dotproduct(ninjaE3,spvae5e2)
      acd149(3)=dotproduct(ninjaE3,spvae2e1)
      acd149(4)=abb149(7)
      acd149(5)=dotproduct(ninjaA1,spvae5e2)
      acd149(6)=dotproduct(ninjaE3,spvae1e5)
      acd149(7)=dotproduct(ninjaA1,spvae2e5)
      acd149(8)=dotproduct(ninjaE3,spvae5e1)
      acd149(9)=dotproduct(ninjaE3,spvae1e2)
      acd149(10)=dotproduct(ninjaA1,spvae5e1)
      acd149(11)=dotproduct(ninjaE3,spvae2e5)
      acd149(12)=dotproduct(ninjaA1,spvae1e2)
      acd149(13)=dotproduct(ninjaA1,spvae2e1)
      acd149(14)=dotproduct(ninjaA0,ninjaE3)
      acd149(15)=abb149(32)
      acd149(16)=abb149(24)
      acd149(17)=abb149(28)
      acd149(18)=abb149(30)
      acd149(19)=dotproduct(ninjaA0,spvae1e5)
      acd149(20)=dotproduct(ninjaA0,spvae5e2)
      acd149(21)=dotproduct(ninjaA0,spvae2e5)
      acd149(22)=dotproduct(ninjaA0,spvae5e1)
      acd149(23)=dotproduct(ninjaA0,spvae1e2)
      acd149(24)=dotproduct(ninjaA0,spvae2e1)
      acd149(25)=abb149(3)
      acd149(26)=dotproduct(ninjaE3,spvak1e2)
      acd149(27)=abb149(16)
      acd149(28)=abb149(21)
      acd149(29)=dotproduct(ninjaE3,spvak2e1)
      acd149(30)=abb149(26)
      acd149(31)=dotproduct(ninjaE3,spvak5e1)
      acd149(32)=abb149(34)
      acd149(33)=dotproduct(ninjaE3,spvae2k1)
      acd149(34)=abb149(8)
      acd149(35)=dotproduct(ninjaE3,spvae1k2)
      acd149(36)=abb149(23)
      acd149(37)=dotproduct(ninjaE3,spvae1k5)
      acd149(38)=abb149(29)
      acd149(39)=abb149(5)
      acd149(40)=abb149(6)
      acd149(41)=abb149(18)
      acd149(42)=abb149(27)
      acd149(43)=abb149(17)
      acd149(44)=abb149(22)
      acd149(45)=abb149(36)
      acd149(46)=abb149(33)
      acd149(47)=-acd149(9)*acd149(7)
      acd149(48)=-acd149(11)*acd149(12)
      acd149(47)=acd149(47)+acd149(48)
      acd149(47)=acd149(8)*acd149(47)
      acd149(48)=-acd149(3)*acd149(5)
      acd149(49)=-acd149(2)*acd149(13)
      acd149(48)=acd149(48)+acd149(49)
      acd149(48)=acd149(6)*acd149(48)
      acd149(49)=acd149(11)*acd149(9)
      acd149(50)=-acd149(10)*acd149(49)
      acd149(51)=acd149(2)*acd149(3)
      acd149(52)=-acd149(1)*acd149(51)
      acd149(47)=acd149(48)+acd149(47)+acd149(50)+acd149(52)
      acd149(47)=acd149(4)*acd149(47)
      acd149(48)=acd149(37)*acd149(46)
      acd149(50)=acd149(35)*acd149(45)
      acd149(52)=acd149(33)*acd149(44)
      acd149(53)=2.0_ki*acd149(14)
      acd149(54)=acd149(18)*acd149(53)
      acd149(55)=-acd149(9)*acd149(43)
      acd149(56)=acd149(11)*acd149(40)
      acd149(48)=acd149(56)+acd149(55)+acd149(54)+acd149(52)+acd149(48)+acd149(&
      &50)
      acd149(48)=acd149(8)*acd149(48)
      acd149(50)=acd149(31)*acd149(32)
      acd149(52)=acd149(29)*acd149(30)
      acd149(54)=acd149(26)*acd149(27)
      acd149(55)=acd149(15)*acd149(53)
      acd149(56)=-acd149(3)*acd149(28)
      acd149(57)=acd149(2)*acd149(25)
      acd149(50)=acd149(57)+acd149(56)+acd149(55)+acd149(54)+acd149(50)+acd149(&
      &52)
      acd149(50)=acd149(6)*acd149(50)
      acd149(52)=acd149(31)*acd149(42)
      acd149(54)=-acd149(29)*acd149(41)
      acd149(55)=-acd149(26)*acd149(39)
      acd149(56)=acd149(17)*acd149(53)
      acd149(52)=acd149(56)+acd149(55)+acd149(52)+acd149(54)
      acd149(52)=acd149(11)*acd149(52)
      acd149(54)=acd149(37)*acd149(38)
      acd149(55)=-acd149(35)*acd149(36)
      acd149(56)=-acd149(33)*acd149(34)
      acd149(53)=acd149(16)*acd149(53)
      acd149(53)=acd149(53)+acd149(56)+acd149(54)+acd149(55)
      acd149(53)=acd149(2)*acd149(53)
      acd149(54)=-acd149(9)*acd149(21)
      acd149(55)=-acd149(11)*acd149(23)
      acd149(54)=acd149(54)+acd149(55)
      acd149(54)=acd149(8)*acd149(54)
      acd149(55)=-acd149(3)*acd149(20)
      acd149(56)=-acd149(2)*acd149(24)
      acd149(55)=acd149(55)+acd149(56)
      acd149(55)=acd149(6)*acd149(55)
      acd149(49)=-acd149(22)*acd149(49)
      acd149(51)=-acd149(19)*acd149(51)
      acd149(49)=acd149(55)+acd149(54)+acd149(49)+acd149(51)
      acd149(49)=acd149(4)*acd149(49)
      acd149(48)=acd149(49)+acd149(50)+acd149(48)+acd149(52)+acd149(53)
      brack(ninjaidxt0x0mu0)=acd149(48)
      brack(ninjaidxt0x0mu2)=0.0_ki
      brack(ninjaidxt0x1mu0)=acd149(47)
      brack(ninjaidxt0x2mu0)=0.0_ki
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d149h0_qp_ninja_t2")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd149h0_qp
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
end module     p2_part21part21_part25part25part21_d149h0l132_qp
