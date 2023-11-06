module     p3_part1part21_part25part25part1_d84h2l132
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p3_part1part2 &
   ! &1_part25part25part1/helicity2d84h2l132.f90
   ! generator: buildfortran_tn3.py
   use p3_part1part21_part25part25part1_config, only: ki
   use p3_part1part21_part25part25part1_util, only: cond_t, d => metric_tensor
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
      use p3_part1part21_part25part25part1_model
      use p3_part1part21_part25part25part1_kinematics
      use p3_part1part21_part25part25part1_color
      use p3_part1part21_part25part25part1_abbrevd84h2
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(1) :: acd84
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
      use p3_part1part21_part25part25part1_model
      use p3_part1part21_part25part25part1_kinematics
      use p3_part1part21_part25part25part1_color
      use p3_part1part21_part25part25part1_abbrevd84h2
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(42) :: acd84
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd84(1)=dotproduct(ninjaA0,ninjaE3)
      acd84(2)=abb84(10)
      acd84(3)=dotproduct(ninjaE3,spvak2e2)
      acd84(4)=abb84(20)
      acd84(5)=dotproduct(ninjaE3,spvae2k2)
      acd84(6)=abb84(28)
      acd84(7)=dotproduct(ninjaE3,spvak2k5)
      acd84(8)=abb84(13)
      acd84(9)=dotproduct(ninjaE3,spvae2l3)
      acd84(10)=abb84(37)
      acd84(11)=dotproduct(ninjaE3,spvak1k2)
      acd84(12)=abb84(26)
      acd84(13)=dotproduct(ninjaE3,spval3e2)
      acd84(14)=abb84(41)
      acd84(15)=dotproduct(ninjaE3,spvae2l4)
      acd84(16)=abb84(30)
      acd84(17)=dotproduct(ninjaE3,spval4e2)
      acd84(18)=abb84(34)
      acd84(19)=abb84(7)
      acd84(20)=abb84(18)
      acd84(21)=abb84(27)
      acd84(22)=abb84(46)
      acd84(23)=abb84(25)
      acd84(24)=abb84(39)
      acd84(25)=abb84(44)
      acd84(26)=dotproduct(ninjaE3,spvak1e2)
      acd84(27)=dotproduct(ninjaE3,spvae2k5)
      acd84(28)=abb84(8)
      acd84(29)=abb84(16)
      acd84(30)=abb84(38)
      acd84(31)=abb84(22)
      acd84(32)=abb84(23)
      acd84(33)=-acd84(17)*acd84(18)
      acd84(34)=acd84(15)*acd84(16)
      acd84(35)=acd84(11)*acd84(12)
      acd84(36)=acd84(7)*acd84(8)
      acd84(37)=acd84(13)*acd84(14)
      acd84(38)=acd84(9)*acd84(10)
      acd84(39)=acd84(5)*acd84(6)
      acd84(40)=acd84(3)*acd84(4)
      acd84(41)=2.0_ki*acd84(1)
      acd84(42)=acd84(2)*acd84(41)
      acd84(33)=acd84(42)+acd84(40)+acd84(39)+acd84(38)+acd84(37)+acd84(36)+acd&
      &84(35)+acd84(33)+acd84(34)
      acd84(33)=acd84(33)*acd84(41)
      acd84(34)=acd84(15)*acd84(22)
      acd84(35)=acd84(11)*acd84(21)
      acd84(36)=acd84(9)*acd84(20)
      acd84(37)=acd84(5)*acd84(19)
      acd84(34)=acd84(37)+acd84(36)+acd84(34)+acd84(35)
      acd84(34)=acd84(3)*acd84(34)
      acd84(35)=acd84(17)*acd84(25)
      acd84(36)=acd84(7)*acd84(23)
      acd84(37)=acd84(13)*acd84(24)
      acd84(35)=acd84(37)+acd84(35)+acd84(36)
      acd84(35)=acd84(5)*acd84(35)
      acd84(36)=acd84(15)*acd84(32)
      acd84(37)=acd84(11)*acd84(31)
      acd84(36)=acd84(36)+acd84(37)
      acd84(36)=acd84(13)*acd84(36)
      acd84(37)=acd84(17)*acd84(30)
      acd84(38)=acd84(7)*acd84(29)
      acd84(37)=acd84(37)+acd84(38)
      acd84(37)=acd84(9)*acd84(37)
      acd84(38)=acd84(26)*acd84(27)*acd84(28)
      acd84(33)=acd84(33)+acd84(34)+acd84(35)+acd84(37)+acd84(38)+acd84(36)
      brack(ninjaidxt0x0mu0)=acd84(33)
      brack(ninjaidxt0x0mu2)=0.0_ki
      brack(ninjaidxt0x1mu0)=0.0_ki
      brack(ninjaidxt0x2mu0)=0.0_ki
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p3_part1part21_part25part25part1_d84h2_ninja_t2")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki => ki_nin
      use p3_part1part21_part25part25part1_globalsl1, only: epspow
      use p3_part1part21_part25part25part1_kinematics
      use p3_part1part21_part25part25part1_abbrevd84h2
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a0, a1, b, c
      complex(ki), dimension(0:2), intent(in) :: param
      complex(ki), dimension(4) :: vecA0, vecA1, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
	     vecA0(1:4) = + a0(0:3)
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
end module     p3_part1part21_part25part25part1_d84h2l132
