module     p3_part1part21_part25part25part1_d88h2l132_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p3_part1part2 &
   ! &1_part25part25part1/helicity2d88h2l132_qp.f90
   ! generator: buildfortran_tn3.py
   use p3_part1part21_part25part25part1_config, only: ki => ki_qp
   use p3_part1part21_part25part25part1_util_qp, only: cond_t, d => metric_tens&
   &or
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
      use p3_part1part21_part25part25part1_model_qp
      use p3_part1part21_part25part25part1_kinematics_qp
      use p3_part1part21_part25part25part1_color_qp
      use p3_part1part21_part25part25part1_abbrevd88h2_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(17) :: acd88
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd88(1)=dotproduct(ninjaA0,ninjaE3)
      acd88(2)=dotproduct(ninjaE3,spvae2k5)
      acd88(3)=dotproduct(ninjaE3,spvak1e2)
      acd88(4)=abb88(9)
      acd88(5)=dotproduct(ninjaE3,spvak1l4)
      acd88(6)=dotproduct(ninjaE3,spval4e2)
      acd88(7)=abb88(7)
      acd88(8)=dotproduct(ninjaE3,spvak1k2)
      acd88(9)=dotproduct(ninjaE3,spvak2e2)
      acd88(10)=abb88(18)
      acd88(11)=dotproduct(ninjaE3,spval4k5)
      acd88(12)=dotproduct(ninjaE3,spvae2l4)
      acd88(13)=dotproduct(ninjaE3,spvak2k5)
      acd88(14)=dotproduct(ninjaE3,spvae2k2)
      acd88(15)=acd88(4)*acd88(1)*acd88(2)
      acd88(16)=acd88(12)*acd88(11)*acd88(7)
      acd88(17)=acd88(14)*acd88(13)*acd88(10)
      acd88(15)=acd88(17)+acd88(16)+2.0_ki*acd88(15)
      acd88(15)=acd88(3)*acd88(15)
      acd88(16)=acd88(6)*acd88(5)*acd88(7)
      acd88(17)=acd88(9)*acd88(8)*acd88(10)
      acd88(16)=acd88(17)+acd88(16)
      acd88(16)=acd88(2)*acd88(16)
      acd88(15)=acd88(15)+acd88(16)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd88(15)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p3_part1part21_part25part25part1_model_qp
      use p3_part1part21_part25part25part1_kinematics_qp
      use p3_part1part21_part25part25part1_color_qp
      use p3_part1part21_part25part25part1_abbrevd88h2_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(67) :: acd88
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd88(1)=dotproduct(ninjaE3,spvae2k5)
      acd88(2)=dotproduct(ninjaE3,spvak1e2)
      acd88(3)=abb88(9)
      acd88(4)=dotproduct(ninjaA0,ninjaA1)
      acd88(5)=dotproduct(ninjaA0,ninjaE3)
      acd88(6)=dotproduct(ninjaA1,spvae2k5)
      acd88(7)=dotproduct(ninjaA1,spvak1e2)
      acd88(8)=dotproduct(ninjaA1,spvak1l4)
      acd88(9)=dotproduct(ninjaE3,spval4e2)
      acd88(10)=abb88(7)
      acd88(11)=dotproduct(ninjaA1,spval4e2)
      acd88(12)=dotproduct(ninjaE3,spvak1l4)
      acd88(13)=dotproduct(ninjaE3,spvak2e2)
      acd88(14)=dotproduct(ninjaE3,spvak1k2)
      acd88(15)=abb88(18)
      acd88(16)=dotproduct(ninjaA1,spval4k5)
      acd88(17)=dotproduct(ninjaE3,spvae2l4)
      acd88(18)=dotproduct(ninjaE3,spval4k5)
      acd88(19)=dotproduct(ninjaE3,spvae2k2)
      acd88(20)=dotproduct(ninjaE3,spvak2k5)
      acd88(21)=dotproduct(ninjaA1,spvae2l4)
      acd88(22)=dotproduct(ninjaA1,spvak2e2)
      acd88(23)=dotproduct(ninjaA1,spvae2k2)
      acd88(24)=dotproduct(ninjaA1,spvak1k2)
      acd88(25)=dotproduct(ninjaA1,spvak2k5)
      acd88(26)=dotproduct(ninjaA1,ninjaA1)
      acd88(27)=dotproduct(ninjaA0,ninjaA0)
      acd88(28)=dotproduct(ninjaA0,spvae2k5)
      acd88(29)=dotproduct(ninjaA0,spvak1e2)
      acd88(30)=abb88(25)
      acd88(31)=abb88(23)
      acd88(32)=dotproduct(ninjaA0,spvak1l4)
      acd88(33)=dotproduct(ninjaA0,spval4e2)
      acd88(34)=dotproduct(ninjaA0,spval4k5)
      acd88(35)=dotproduct(ninjaA0,spvae2l4)
      acd88(36)=dotproduct(ninjaA0,spvak2e2)
      acd88(37)=dotproduct(ninjaA0,spvae2k2)
      acd88(38)=dotproduct(ninjaA0,spvak1k2)
      acd88(39)=dotproduct(ninjaA0,spvak2k5)
      acd88(40)=abb88(28)
      acd88(41)=abb88(14)
      acd88(42)=abb88(30)
      acd88(43)=abb88(13)
      acd88(44)=abb88(21)
      acd88(45)=abb88(20)
      acd88(46)=abb88(16)
      acd88(47)=abb88(32)
      acd88(48)=abb88(22)
      acd88(49)=acd88(3)*acd88(2)
      acd88(50)=acd88(49)*acd88(1)
      acd88(51)=ninjaP1+2.0_ki*acd88(4)
      acd88(51)=acd88(49)*acd88(51)
      acd88(52)=acd88(9)*acd88(10)
      acd88(53)=acd88(8)*acd88(52)
      acd88(54)=acd88(12)*acd88(10)
      acd88(55)=acd88(11)*acd88(54)
      acd88(56)=acd88(14)*acd88(15)
      acd88(57)=acd88(22)*acd88(56)
      acd88(58)=acd88(13)*acd88(15)
      acd88(59)=acd88(24)*acd88(58)
      acd88(51)=acd88(59)+acd88(57)+acd88(55)+acd88(53)+acd88(51)
      acd88(51)=acd88(1)*acd88(51)
      acd88(53)=acd88(17)*acd88(10)
      acd88(55)=acd88(16)*acd88(53)
      acd88(57)=acd88(18)*acd88(10)
      acd88(59)=acd88(21)*acd88(57)
      acd88(60)=acd88(20)*acd88(15)
      acd88(61)=acd88(23)*acd88(60)
      acd88(62)=acd88(19)*acd88(15)
      acd88(63)=acd88(25)*acd88(62)
      acd88(55)=acd88(63)+acd88(61)+acd88(59)+acd88(55)
      acd88(55)=acd88(2)*acd88(55)
      acd88(59)=2.0_ki*acd88(5)
      acd88(61)=acd88(59)*acd88(49)
      acd88(63)=acd88(56)*acd88(13)
      acd88(64)=acd88(54)*acd88(9)
      acd88(61)=acd88(64)+acd88(61)+acd88(63)
      acd88(63)=acd88(6)*acd88(61)
      acd88(64)=acd88(59)*acd88(3)*acd88(1)
      acd88(65)=acd88(60)*acd88(19)
      acd88(66)=acd88(57)*acd88(17)
      acd88(64)=acd88(66)+acd88(64)+acd88(65)
      acd88(65)=acd88(7)*acd88(64)
      acd88(51)=acd88(63)+acd88(65)+acd88(51)+acd88(55)
      acd88(55)=acd88(26)+ninjaP2
      acd88(55)=acd88(50)*acd88(55)
      acd88(63)=acd88(31)*acd88(59)
      acd88(53)=acd88(34)*acd88(53)
      acd88(57)=acd88(35)*acd88(57)
      acd88(60)=acd88(37)*acd88(60)
      acd88(62)=acd88(39)*acd88(62)
      acd88(65)=acd88(42)*acd88(1)
      acd88(66)=acd88(45)*acd88(17)
      acd88(67)=acd88(46)*acd88(19)
      acd88(53)=acd88(65)+acd88(67)+acd88(66)+acd88(62)+acd88(60)+acd88(57)+acd&
      &88(53)+acd88(63)
      acd88(53)=acd88(2)*acd88(53)
      acd88(57)=ninjaP0+acd88(27)
      acd88(49)=acd88(49)*acd88(57)
      acd88(57)=acd88(30)*acd88(59)
      acd88(52)=acd88(32)*acd88(52)
      acd88(54)=acd88(33)*acd88(54)
      acd88(56)=acd88(36)*acd88(56)
      acd88(58)=acd88(38)*acd88(58)
      acd88(59)=acd88(41)*acd88(9)
      acd88(60)=acd88(43)*acd88(13)
      acd88(49)=acd88(60)+acd88(59)+acd88(58)+acd88(56)+acd88(54)+acd88(52)+acd&
      &88(57)+acd88(49)
      acd88(49)=acd88(1)*acd88(49)
      acd88(52)=acd88(28)*acd88(61)
      acd88(54)=acd88(29)*acd88(64)
      acd88(56)=acd88(40)*acd88(12)*acd88(9)
      acd88(57)=acd88(44)*acd88(18)*acd88(17)
      acd88(58)=acd88(47)*acd88(14)*acd88(13)
      acd88(59)=acd88(48)*acd88(20)*acd88(19)
      acd88(49)=acd88(59)+acd88(58)+acd88(57)+acd88(56)+acd88(52)+acd88(54)+acd&
      &88(49)+acd88(53)
      brack(ninjaidxt0x0mu0)=acd88(49)
      brack(ninjaidxt0x0mu2)=acd88(50)
      brack(ninjaidxt0x1mu0)=acd88(51)
      brack(ninjaidxt0x2mu0)=acd88(55)
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p3_part1part21_part25part25part1_d88h2_qp_ninja_t2")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p3_part1part21_part25part25part1_globalsl1_qp, only: epspow
      use p3_part1part21_part25part25part1_kinematics_qp
      use p3_part1part21_part25part25part1_abbrevd88h2_qp
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a0, a1, b, c
      complex(ki), dimension(0:2), intent(in) :: param
      complex(ki), dimension(4) :: vecA0, vecA1, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
	     vecA0(1:4) = - a0(0:3)
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
end module     p3_part1part21_part25part25part1_d88h2l132_qp
