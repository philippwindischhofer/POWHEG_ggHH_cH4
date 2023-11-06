module     p2_part21part21_part25part25part21_d147h0l132_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d147h0l132_qp.f90
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
      use p2_part21part21_part25part25part21_abbrevd147h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(9) :: acd147
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd147(1)=dotproduct(ninjaE3,spvae1e2)
      acd147(2)=dotproduct(ninjaE3,spvae5e1)
      acd147(3)=dotproduct(ninjaE3,spvae2e5)
      acd147(4)=abb147(7)
      acd147(5)=dotproduct(ninjaE3,spvae2e1)
      acd147(6)=dotproduct(ninjaE3,spvae1e5)
      acd147(7)=dotproduct(ninjaE3,spvae5e2)
      acd147(8)=acd147(3)*acd147(2)*acd147(1)
      acd147(9)=acd147(7)*acd147(6)*acd147(5)
      acd147(8)=acd147(8)+acd147(9)
      acd147(8)=acd147(4)*acd147(8)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd147(8)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd147h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(57) :: acd147
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd147(1)=dotproduct(ninjaA1,spvae1e2)
      acd147(2)=dotproduct(ninjaE3,spvae2e5)
      acd147(3)=dotproduct(ninjaE3,spvae5e1)
      acd147(4)=abb147(7)
      acd147(5)=dotproduct(ninjaA1,spvae2e5)
      acd147(6)=dotproduct(ninjaE3,spvae1e2)
      acd147(7)=dotproduct(ninjaA1,spvae5e2)
      acd147(8)=dotproduct(ninjaE3,spvae2e1)
      acd147(9)=dotproduct(ninjaE3,spvae1e5)
      acd147(10)=dotproduct(ninjaA1,spvae2e1)
      acd147(11)=dotproduct(ninjaE3,spvae5e2)
      acd147(12)=dotproduct(ninjaA1,spvae5e1)
      acd147(13)=dotproduct(ninjaA1,spvae1e5)
      acd147(14)=dotproduct(ninjaA0,ninjaE3)
      acd147(15)=abb147(32)
      acd147(16)=abb147(28)
      acd147(17)=abb147(25)
      acd147(18)=abb147(29)
      acd147(19)=dotproduct(ninjaA0,spvae1e2)
      acd147(20)=dotproduct(ninjaA0,spvae2e5)
      acd147(21)=dotproduct(ninjaA0,spvae5e2)
      acd147(22)=dotproduct(ninjaA0,spvae2e1)
      acd147(23)=dotproduct(ninjaA0,spvae5e1)
      acd147(24)=dotproduct(ninjaA0,spvae1e5)
      acd147(25)=abb147(3)
      acd147(26)=dotproduct(ninjaE3,spvak1e5)
      acd147(27)=abb147(16)
      acd147(28)=abb147(21)
      acd147(29)=dotproduct(ninjaE3,spvak5e1)
      acd147(30)=abb147(23)
      acd147(31)=dotproduct(ninjaE3,spvak2e1)
      acd147(32)=abb147(36)
      acd147(33)=dotproduct(ninjaE3,spvae5k1)
      acd147(34)=abb147(8)
      acd147(35)=dotproduct(ninjaE3,spvae1k5)
      acd147(36)=abb147(14)
      acd147(37)=dotproduct(ninjaE3,spvae1k2)
      acd147(38)=abb147(27)
      acd147(39)=abb147(5)
      acd147(40)=abb147(6)
      acd147(41)=abb147(18)
      acd147(42)=abb147(24)
      acd147(43)=abb147(17)
      acd147(44)=abb147(26)
      acd147(45)=abb147(34)
      acd147(46)=abb147(35)
      acd147(47)=acd147(9)*acd147(7)
      acd147(48)=acd147(11)*acd147(13)
      acd147(47)=acd147(47)+acd147(48)
      acd147(47)=acd147(8)*acd147(47)
      acd147(48)=acd147(3)*acd147(5)
      acd147(49)=acd147(2)*acd147(12)
      acd147(48)=acd147(48)+acd147(49)
      acd147(48)=acd147(6)*acd147(48)
      acd147(49)=acd147(11)*acd147(9)
      acd147(50)=acd147(10)*acd147(49)
      acd147(51)=acd147(2)*acd147(3)
      acd147(52)=acd147(1)*acd147(51)
      acd147(47)=acd147(48)+acd147(47)+acd147(50)+acd147(52)
      acd147(47)=acd147(4)*acd147(47)
      acd147(48)=acd147(37)*acd147(46)
      acd147(50)=acd147(35)*acd147(45)
      acd147(52)=acd147(33)*acd147(44)
      acd147(53)=2.0_ki*acd147(14)
      acd147(54)=acd147(18)*acd147(53)
      acd147(55)=acd147(9)*acd147(43)
      acd147(56)=acd147(11)*acd147(40)
      acd147(48)=acd147(56)+acd147(55)+acd147(54)+acd147(52)+acd147(48)+acd147(&
      &50)
      acd147(48)=acd147(8)*acd147(48)
      acd147(50)=acd147(31)*acd147(32)
      acd147(52)=acd147(29)*acd147(30)
      acd147(54)=acd147(26)*acd147(27)
      acd147(55)=acd147(15)*acd147(53)
      acd147(56)=acd147(3)*acd147(28)
      acd147(57)=acd147(2)*acd147(25)
      acd147(50)=acd147(57)+acd147(56)+acd147(55)+acd147(54)+acd147(50)+acd147(&
      &52)
      acd147(50)=acd147(6)*acd147(50)
      acd147(52)=acd147(31)*acd147(42)
      acd147(54)=-acd147(29)*acd147(41)
      acd147(55)=acd147(26)*acd147(39)
      acd147(56)=acd147(17)*acd147(53)
      acd147(52)=acd147(56)+acd147(55)+acd147(52)+acd147(54)
      acd147(52)=acd147(11)*acd147(52)
      acd147(54)=acd147(37)*acd147(38)
      acd147(55)=-acd147(35)*acd147(36)
      acd147(56)=acd147(33)*acd147(34)
      acd147(53)=acd147(16)*acd147(53)
      acd147(53)=acd147(53)+acd147(56)+acd147(54)+acd147(55)
      acd147(53)=acd147(2)*acd147(53)
      acd147(54)=acd147(9)*acd147(21)
      acd147(55)=acd147(11)*acd147(24)
      acd147(54)=acd147(54)+acd147(55)
      acd147(54)=acd147(8)*acd147(54)
      acd147(55)=acd147(3)*acd147(20)
      acd147(56)=acd147(2)*acd147(23)
      acd147(55)=acd147(55)+acd147(56)
      acd147(55)=acd147(6)*acd147(55)
      acd147(49)=acd147(22)*acd147(49)
      acd147(51)=acd147(19)*acd147(51)
      acd147(49)=acd147(55)+acd147(54)+acd147(49)+acd147(51)
      acd147(49)=acd147(4)*acd147(49)
      acd147(48)=acd147(49)+acd147(50)+acd147(48)+acd147(52)+acd147(53)
      brack(ninjaidxt0x0mu0)=acd147(48)
      brack(ninjaidxt0x0mu2)=0.0_ki
      brack(ninjaidxt0x1mu0)=acd147(47)
      brack(ninjaidxt0x2mu0)=0.0_ki
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d147h0_qp_ninja_t2")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd147h0_qp
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
end module     p2_part21part21_part25part25part21_d147h0l132_qp
