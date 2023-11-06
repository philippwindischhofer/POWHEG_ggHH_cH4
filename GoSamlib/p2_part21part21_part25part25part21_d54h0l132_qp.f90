module     p2_part21part21_part25part25part21_d54h0l132_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d54h0l132_qp.f90
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
      use p2_part21part21_part25part25part21_abbrevd54h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(9) :: acd54
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd54(1)=dotproduct(ninjaE3,spvae1e2)
      acd54(2)=dotproduct(ninjaE3,spvae5e1)
      acd54(3)=dotproduct(ninjaE3,spvae2e5)
      acd54(4)=abb54(7)
      acd54(5)=dotproduct(ninjaE3,spvae2e1)
      acd54(6)=dotproduct(ninjaE3,spvae1e5)
      acd54(7)=dotproduct(ninjaE3,spvae5e2)
      acd54(8)=-acd54(3)*acd54(2)*acd54(1)
      acd54(9)=-acd54(7)*acd54(6)*acd54(5)
      acd54(8)=acd54(8)+acd54(9)
      acd54(8)=acd54(4)*acd54(8)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd54(8)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd54h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(57) :: acd54
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd54(1)=dotproduct(ninjaA1,spvae2e1)
      acd54(2)=dotproduct(ninjaE3,spvae5e2)
      acd54(3)=dotproduct(ninjaE3,spvae1e5)
      acd54(4)=abb54(7)
      acd54(5)=dotproduct(ninjaA1,spvae5e2)
      acd54(6)=dotproduct(ninjaE3,spvae2e1)
      acd54(7)=dotproduct(ninjaA1,spvae1e2)
      acd54(8)=dotproduct(ninjaE3,spvae2e5)
      acd54(9)=dotproduct(ninjaE3,spvae5e1)
      acd54(10)=dotproduct(ninjaA1,spvae2e5)
      acd54(11)=dotproduct(ninjaE3,spvae1e2)
      acd54(12)=dotproduct(ninjaA1,spvae5e1)
      acd54(13)=dotproduct(ninjaA1,spvae1e5)
      acd54(14)=dotproduct(ninjaA0,ninjaE3)
      acd54(15)=abb54(30)
      acd54(16)=abb54(25)
      acd54(17)=abb54(23)
      acd54(18)=abb54(9)
      acd54(19)=dotproduct(ninjaA0,spvae2e1)
      acd54(20)=dotproduct(ninjaA0,spvae5e2)
      acd54(21)=dotproduct(ninjaA0,spvae1e2)
      acd54(22)=dotproduct(ninjaA0,spvae2e5)
      acd54(23)=dotproduct(ninjaA0,spvae5e1)
      acd54(24)=dotproduct(ninjaA0,spvae1e5)
      acd54(25)=abb54(4)
      acd54(26)=abb54(29)
      acd54(27)=dotproduct(ninjaE3,spvae5k1)
      acd54(28)=abb54(26)
      acd54(29)=dotproduct(ninjaE3,spvae1k2)
      acd54(30)=abb54(34)
      acd54(31)=dotproduct(ninjaE3,spvae1k5)
      acd54(32)=abb54(32)
      acd54(33)=dotproduct(ninjaE3,spvak1e5)
      acd54(34)=abb54(16)
      acd54(35)=dotproduct(ninjaE3,spvak5e1)
      acd54(36)=abb54(19)
      acd54(37)=dotproduct(ninjaE3,spvak2e1)
      acd54(38)=abb54(24)
      acd54(39)=abb54(5)
      acd54(40)=abb54(6)
      acd54(41)=abb54(18)
      acd54(42)=abb54(33)
      acd54(43)=abb54(35)
      acd54(44)=abb54(14)
      acd54(45)=abb54(15)
      acd54(46)=abb54(21)
      acd54(47)=-acd54(9)*acd54(10)
      acd54(48)=-acd54(8)*acd54(12)
      acd54(47)=acd54(47)+acd54(48)
      acd54(47)=acd54(11)*acd54(47)
      acd54(48)=-acd54(3)*acd54(5)
      acd54(49)=-acd54(2)*acd54(13)
      acd54(48)=acd54(48)+acd54(49)
      acd54(48)=acd54(6)*acd54(48)
      acd54(49)=acd54(8)*acd54(9)
      acd54(50)=-acd54(7)*acd54(49)
      acd54(51)=acd54(2)*acd54(3)
      acd54(52)=-acd54(1)*acd54(51)
      acd54(47)=acd54(48)+acd54(47)+acd54(50)+acd54(52)
      acd54(47)=acd54(4)*acd54(47)
      acd54(48)=acd54(37)*acd54(43)
      acd54(50)=acd54(35)*acd54(42)
      acd54(52)=acd54(33)*acd54(39)
      acd54(53)=2.0_ki*acd54(14)
      acd54(54)=acd54(17)*acd54(53)
      acd54(55)=-acd54(9)*acd54(41)
      acd54(56)=acd54(8)*acd54(40)
      acd54(48)=acd54(56)+acd54(55)+acd54(54)+acd54(52)+acd54(48)+acd54(50)
      acd54(48)=acd54(11)*acd54(48)
      acd54(50)=acd54(31)*acd54(32)
      acd54(52)=acd54(29)*acd54(30)
      acd54(54)=acd54(27)*acd54(28)
      acd54(55)=acd54(15)*acd54(53)
      acd54(56)=-acd54(3)*acd54(26)
      acd54(57)=acd54(2)*acd54(25)
      acd54(50)=acd54(57)+acd54(56)+acd54(55)+acd54(54)+acd54(50)+acd54(52)
      acd54(50)=acd54(6)*acd54(50)
      acd54(52)=acd54(31)*acd54(46)
      acd54(54)=acd54(29)*acd54(45)
      acd54(55)=-acd54(27)*acd54(44)
      acd54(56)=acd54(18)*acd54(53)
      acd54(52)=acd54(56)+acd54(55)+acd54(52)+acd54(54)
      acd54(52)=acd54(8)*acd54(52)
      acd54(54)=acd54(37)*acd54(38)
      acd54(55)=acd54(35)*acd54(36)
      acd54(56)=-acd54(33)*acd54(34)
      acd54(53)=acd54(16)*acd54(53)
      acd54(53)=acd54(53)+acd54(56)+acd54(54)+acd54(55)
      acd54(53)=acd54(2)*acd54(53)
      acd54(54)=-acd54(9)*acd54(22)
      acd54(55)=-acd54(8)*acd54(23)
      acd54(54)=acd54(54)+acd54(55)
      acd54(54)=acd54(11)*acd54(54)
      acd54(55)=-acd54(3)*acd54(20)
      acd54(56)=-acd54(2)*acd54(24)
      acd54(55)=acd54(55)+acd54(56)
      acd54(55)=acd54(6)*acd54(55)
      acd54(49)=-acd54(21)*acd54(49)
      acd54(51)=-acd54(19)*acd54(51)
      acd54(49)=acd54(55)+acd54(54)+acd54(49)+acd54(51)
      acd54(49)=acd54(4)*acd54(49)
      acd54(48)=acd54(49)+acd54(50)+acd54(48)+acd54(52)+acd54(53)
      brack(ninjaidxt0x0mu0)=acd54(48)
      brack(ninjaidxt0x0mu2)=0.0_ki
      brack(ninjaidxt0x1mu0)=acd54(47)
      brack(ninjaidxt0x2mu0)=0.0_ki
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d54h0_qp_ninja_t2")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd54h0_qp
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
end module     p2_part21part21_part25part25part21_d54h0l132_qp
