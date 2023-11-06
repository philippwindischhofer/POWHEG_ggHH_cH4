module     p10_part1part21_part25part25part1_d86h2l132_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p10_part1part &
   ! &21_part25part25part1/helicity2d86h2l132_qp.f90
   ! generator: buildfortran_tn3.py
   use p10_part1part21_part25part25part1_config, only: ki => ki_qp
   use p10_part1part21_part25part25part1_util_qp, only: cond_t, d => metric_ten&
   &sor
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
      use p10_part1part21_part25part25part1_model_qp
      use p10_part1part21_part25part25part1_kinematics_qp
      use p10_part1part21_part25part25part1_color_qp
      use p10_part1part21_part25part25part1_abbrevd86h2_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(17) :: acd86
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd86(1)=dotproduct(ninjaA0,ninjaE3)
      acd86(2)=dotproduct(ninjaE3,spvae2k5)
      acd86(3)=dotproduct(ninjaE3,spvak1e2)
      acd86(4)=abb86(9)
      acd86(5)=dotproduct(ninjaE3,spvak1l3)
      acd86(6)=dotproduct(ninjaE3,spval3e2)
      acd86(7)=abb86(7)
      acd86(8)=dotproduct(ninjaE3,spvak1k2)
      acd86(9)=dotproduct(ninjaE3,spvak2e2)
      acd86(10)=abb86(19)
      acd86(11)=dotproduct(ninjaE3,spval3k5)
      acd86(12)=dotproduct(ninjaE3,spvae2l3)
      acd86(13)=dotproduct(ninjaE3,spvak2k5)
      acd86(14)=dotproduct(ninjaE3,spvae2k2)
      acd86(15)=acd86(4)*acd86(1)*acd86(2)
      acd86(16)=acd86(12)*acd86(11)*acd86(7)
      acd86(17)=acd86(14)*acd86(13)*acd86(10)
      acd86(15)=acd86(17)+acd86(16)+2.0_ki*acd86(15)
      acd86(15)=acd86(3)*acd86(15)
      acd86(16)=acd86(6)*acd86(5)*acd86(7)
      acd86(17)=acd86(9)*acd86(8)*acd86(10)
      acd86(16)=acd86(17)+acd86(16)
      acd86(16)=acd86(2)*acd86(16)
      acd86(15)=acd86(15)+acd86(16)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd86(15)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p10_part1part21_part25part25part1_model_qp
      use p10_part1part21_part25part25part1_kinematics_qp
      use p10_part1part21_part25part25part1_color_qp
      use p10_part1part21_part25part25part1_abbrevd86h2_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(67) :: acd86
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd86(1)=dotproduct(ninjaE3,spvae2k5)
      acd86(2)=dotproduct(ninjaE3,spvak1e2)
      acd86(3)=abb86(9)
      acd86(4)=dotproduct(ninjaA0,ninjaA1)
      acd86(5)=dotproduct(ninjaA0,ninjaE3)
      acd86(6)=dotproduct(ninjaA1,spvae2k5)
      acd86(7)=dotproduct(ninjaA1,spvak1e2)
      acd86(8)=dotproduct(ninjaA1,spvak1l3)
      acd86(9)=dotproduct(ninjaE3,spval3e2)
      acd86(10)=abb86(7)
      acd86(11)=dotproduct(ninjaA1,spval3e2)
      acd86(12)=dotproduct(ninjaE3,spvak1l3)
      acd86(13)=dotproduct(ninjaE3,spvak2e2)
      acd86(14)=dotproduct(ninjaE3,spvak1k2)
      acd86(15)=abb86(19)
      acd86(16)=dotproduct(ninjaA1,spval3k5)
      acd86(17)=dotproduct(ninjaE3,spvae2l3)
      acd86(18)=dotproduct(ninjaE3,spval3k5)
      acd86(19)=dotproduct(ninjaE3,spvae2k2)
      acd86(20)=dotproduct(ninjaE3,spvak2k5)
      acd86(21)=dotproduct(ninjaA1,spvae2l3)
      acd86(22)=dotproduct(ninjaA1,spvae2k2)
      acd86(23)=dotproduct(ninjaA1,spvak2e2)
      acd86(24)=dotproduct(ninjaA1,spvak1k2)
      acd86(25)=dotproduct(ninjaA1,spvak2k5)
      acd86(26)=dotproduct(ninjaA1,ninjaA1)
      acd86(27)=dotproduct(ninjaA0,ninjaA0)
      acd86(28)=dotproduct(ninjaA0,spvae2k5)
      acd86(29)=dotproduct(ninjaA0,spvak1e2)
      acd86(30)=abb86(26)
      acd86(31)=abb86(24)
      acd86(32)=dotproduct(ninjaA0,spvak1l3)
      acd86(33)=dotproduct(ninjaA0,spval3e2)
      acd86(34)=dotproduct(ninjaA0,spval3k5)
      acd86(35)=dotproduct(ninjaA0,spvae2l3)
      acd86(36)=dotproduct(ninjaA0,spvae2k2)
      acd86(37)=dotproduct(ninjaA0,spvak2e2)
      acd86(38)=dotproduct(ninjaA0,spvak1k2)
      acd86(39)=dotproduct(ninjaA0,spvak2k5)
      acd86(40)=abb86(29)
      acd86(41)=abb86(23)
      acd86(42)=abb86(30)
      acd86(43)=abb86(16)
      acd86(44)=abb86(18)
      acd86(45)=abb86(17)
      acd86(46)=abb86(11)
      acd86(47)=abb86(22)
      acd86(48)=abb86(33)
      acd86(49)=acd86(3)*acd86(2)
      acd86(50)=acd86(49)*acd86(1)
      acd86(51)=ninjaP1+2.0_ki*acd86(4)
      acd86(51)=acd86(49)*acd86(51)
      acd86(52)=acd86(9)*acd86(10)
      acd86(53)=acd86(8)*acd86(52)
      acd86(54)=acd86(12)*acd86(10)
      acd86(55)=acd86(11)*acd86(54)
      acd86(56)=acd86(14)*acd86(15)
      acd86(57)=acd86(23)*acd86(56)
      acd86(58)=acd86(13)*acd86(15)
      acd86(59)=acd86(24)*acd86(58)
      acd86(51)=acd86(59)+acd86(57)+acd86(55)+acd86(53)+acd86(51)
      acd86(51)=acd86(1)*acd86(51)
      acd86(53)=acd86(17)*acd86(10)
      acd86(55)=acd86(16)*acd86(53)
      acd86(57)=acd86(18)*acd86(10)
      acd86(59)=acd86(21)*acd86(57)
      acd86(60)=acd86(20)*acd86(15)
      acd86(61)=acd86(22)*acd86(60)
      acd86(62)=acd86(19)*acd86(15)
      acd86(63)=acd86(25)*acd86(62)
      acd86(55)=acd86(63)+acd86(61)+acd86(59)+acd86(55)
      acd86(55)=acd86(2)*acd86(55)
      acd86(59)=2.0_ki*acd86(5)
      acd86(61)=acd86(59)*acd86(49)
      acd86(63)=acd86(56)*acd86(13)
      acd86(64)=acd86(54)*acd86(9)
      acd86(61)=acd86(64)+acd86(61)+acd86(63)
      acd86(63)=acd86(6)*acd86(61)
      acd86(64)=acd86(59)*acd86(3)*acd86(1)
      acd86(65)=acd86(60)*acd86(19)
      acd86(66)=acd86(57)*acd86(17)
      acd86(64)=acd86(66)+acd86(64)+acd86(65)
      acd86(65)=acd86(7)*acd86(64)
      acd86(51)=acd86(63)+acd86(65)+acd86(51)+acd86(55)
      acd86(55)=acd86(26)+ninjaP2
      acd86(55)=acd86(50)*acd86(55)
      acd86(63)=acd86(31)*acd86(59)
      acd86(53)=acd86(34)*acd86(53)
      acd86(57)=acd86(35)*acd86(57)
      acd86(60)=acd86(36)*acd86(60)
      acd86(62)=acd86(39)*acd86(62)
      acd86(65)=acd86(42)*acd86(1)
      acd86(66)=acd86(45)*acd86(17)
      acd86(67)=acd86(46)*acd86(19)
      acd86(53)=acd86(65)+acd86(67)+acd86(66)+acd86(62)+acd86(60)+acd86(57)+acd&
      &86(53)+acd86(63)
      acd86(53)=acd86(2)*acd86(53)
      acd86(57)=ninjaP0+acd86(27)
      acd86(49)=acd86(49)*acd86(57)
      acd86(57)=acd86(30)*acd86(59)
      acd86(52)=acd86(32)*acd86(52)
      acd86(54)=acd86(33)*acd86(54)
      acd86(56)=acd86(37)*acd86(56)
      acd86(58)=acd86(38)*acd86(58)
      acd86(59)=acd86(41)*acd86(9)
      acd86(60)=acd86(43)*acd86(13)
      acd86(49)=acd86(60)+acd86(59)+acd86(58)+acd86(56)+acd86(54)+acd86(52)+acd&
      &86(57)+acd86(49)
      acd86(49)=acd86(1)*acd86(49)
      acd86(52)=acd86(28)*acd86(61)
      acd86(54)=acd86(29)*acd86(64)
      acd86(56)=acd86(40)*acd86(12)*acd86(9)
      acd86(57)=acd86(44)*acd86(18)*acd86(17)
      acd86(58)=acd86(47)*acd86(20)*acd86(19)
      acd86(59)=acd86(48)*acd86(14)*acd86(13)
      acd86(49)=acd86(59)+acd86(58)+acd86(57)+acd86(56)+acd86(52)+acd86(54)+acd&
      &86(49)+acd86(53)
      brack(ninjaidxt0x0mu0)=acd86(49)
      brack(ninjaidxt0x0mu2)=acd86(50)
      brack(ninjaidxt0x1mu0)=acd86(51)
      brack(ninjaidxt0x2mu0)=acd86(55)
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p10_part1part21_part25part25part1_d86h2_qp_ninja_t2")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p10_part1part21_part25part25part1_globalsl1_qp, only: epspow
      use p10_part1part21_part25part25part1_kinematics_qp
      use p10_part1part21_part25part25part1_abbrevd86h2_qp
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
end module     p10_part1part21_part25part25part1_d86h2l132_qp
