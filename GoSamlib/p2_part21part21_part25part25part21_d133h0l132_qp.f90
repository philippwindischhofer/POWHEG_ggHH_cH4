module     p2_part21part21_part25part25part21_d133h0l132_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d133h0l132_qp.f90
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
      use p2_part21part21_part25part25part21_abbrevd133h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(1) :: acd133
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
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd133h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(60) :: acd133
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd133(1)=dotproduct(ninjaA0,ninjaE3)
      acd133(2)=abb133(12)
      acd133(3)=dotproduct(ninjaE3,spvak2k1)
      acd133(4)=abb133(51)
      acd133(5)=dotproduct(ninjaE3,spval3e1)
      acd133(6)=abb133(45)
      acd133(7)=dotproduct(ninjaE3,spvak5k1)
      acd133(8)=abb133(8)
      acd133(9)=dotproduct(ninjaE3,spvak2e1)
      acd133(10)=abb133(17)
      acd133(11)=dotproduct(ninjaE3,spvae1k2)
      acd133(12)=abb133(15)
      acd133(13)=dotproduct(ninjaE3,spvae1l4)
      acd133(14)=abb133(28)
      acd133(15)=dotproduct(ninjaE3,spvae1l3)
      acd133(16)=abb133(40)
      acd133(17)=dotproduct(ninjaE3,spvae1k1)
      acd133(18)=abb133(26)
      acd133(19)=dotproduct(ninjaE3,spvak1k2)
      acd133(20)=abb133(27)
      acd133(21)=dotproduct(ninjaE3,spval4e1)
      acd133(22)=abb133(32)
      acd133(23)=dotproduct(ninjaE3,spvak1k5)
      acd133(24)=abb133(52)
      acd133(25)=dotproduct(ninjaE3,spvak1e1)
      acd133(26)=abb133(42)
      acd133(27)=abb133(7)
      acd133(28)=abb133(50)
      acd133(29)=abb133(36)
      acd133(30)=abb133(35)
      acd133(31)=abb133(34)
      acd133(32)=abb133(20)
      acd133(33)=dotproduct(ninjaE3,spvak5e1)
      acd133(34)=dotproduct(ninjaE3,spvae1k5)
      acd133(35)=abb133(9)
      acd133(36)=abb133(10)
      acd133(37)=abb133(23)
      acd133(38)=abb133(24)
      acd133(39)=abb133(37)
      acd133(40)=abb133(41)
      acd133(41)=abb133(39)
      acd133(42)=abb133(47)
      acd133(43)=abb133(30)
      acd133(44)=abb133(33)
      acd133(45)=abb133(55)
      acd133(46)=abb133(46)
      acd133(47)=acd133(19)*acd133(20)
      acd133(48)=acd133(3)*acd133(4)
      acd133(49)=acd133(25)*acd133(26)
      acd133(50)=acd133(23)*acd133(24)
      acd133(51)=acd133(21)*acd133(22)
      acd133(52)=acd133(17)*acd133(18)
      acd133(53)=acd133(13)*acd133(14)
      acd133(54)=acd133(7)*acd133(8)
      acd133(55)=acd133(15)*acd133(16)
      acd133(56)=acd133(11)*acd133(12)
      acd133(57)=acd133(9)*acd133(10)
      acd133(58)=acd133(5)*acd133(6)
      acd133(59)=2.0_ki*acd133(1)
      acd133(60)=acd133(2)*acd133(59)
      acd133(47)=acd133(60)+acd133(58)+acd133(57)+acd133(56)+acd133(55)+acd133(&
      &54)+acd133(53)+acd133(52)+acd133(51)+acd133(50)+acd133(49)+acd133(47)+ac&
      &d133(48)
      acd133(47)=acd133(47)*acd133(59)
      acd133(48)=acd133(17)*acd133(39)
      acd133(49)=acd133(13)*acd133(37)
      acd133(50)=acd133(7)*acd133(32)
      acd133(51)=acd133(15)*acd133(38)
      acd133(52)=acd133(11)*acd133(36)
      acd133(48)=acd133(52)+acd133(51)+acd133(50)+acd133(48)+acd133(49)
      acd133(48)=acd133(9)*acd133(48)
      acd133(49)=acd133(3)*acd133(27)
      acd133(50)=acd133(17)*acd133(31)
      acd133(51)=acd133(13)*acd133(30)
      acd133(52)=acd133(7)*acd133(28)
      acd133(53)=acd133(11)*acd133(29)
      acd133(49)=acd133(53)+acd133(52)+acd133(51)+acd133(49)+acd133(50)
      acd133(49)=acd133(5)*acd133(49)
      acd133(50)=acd133(19)*acd133(43)
      acd133(51)=acd133(25)*acd133(46)
      acd133(52)=acd133(23)*acd133(45)
      acd133(53)=acd133(21)*acd133(44)
      acd133(50)=acd133(53)+acd133(52)+acd133(50)+acd133(51)
      acd133(50)=acd133(15)*acd133(50)
      acd133(51)=acd133(25)*acd133(42)
      acd133(52)=acd133(23)*acd133(41)
      acd133(53)=acd133(21)*acd133(40)
      acd133(51)=acd133(53)+acd133(51)+acd133(52)
      acd133(51)=acd133(11)*acd133(51)
      acd133(52)=acd133(33)*acd133(34)*acd133(35)
      acd133(47)=acd133(47)+acd133(49)+acd133(48)+acd133(51)+acd133(52)+acd133(&
      &50)
      brack(ninjaidxt0x0mu0)=acd133(47)
      brack(ninjaidxt0x0mu2)=0.0_ki
      brack(ninjaidxt0x1mu0)=0.0_ki
      brack(ninjaidxt0x2mu0)=0.0_ki
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d133h0_qp_ninja_t2")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd133h0_qp
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
end module     p2_part21part21_part25part25part21_d133h0l132_qp
