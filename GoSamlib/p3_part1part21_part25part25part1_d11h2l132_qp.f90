module     p3_part1part21_part25part25part1_d11h2l132_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p3_part1part2 &
   ! &1_part25part25part1/helicity2d11h2l132_qp.f90
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
      use p3_part1part21_part25part25part1_abbrevd11h2_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(4) :: acd11
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd11(1)=dotproduct(ninjaE3,spvak1e2)
      acd11(2)=dotproduct(ninjaE3,spvae2k5)
      acd11(3)=abb11(4)
      acd11(4)=acd11(3)*acd11(2)*acd11(1)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd11(4)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p3_part1part21_part25part25part1_model_qp
      use p3_part1part21_part25part25part1_kinematics_qp
      use p3_part1part21_part25part25part1_color_qp
      use p3_part1part21_part25part25part1_abbrevd11h2_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(23) :: acd11
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd11(1)=dotproduct(ninjaA1,spvak1e2)
      acd11(2)=dotproduct(ninjaE3,spvae2k5)
      acd11(3)=abb11(4)
      acd11(4)=dotproduct(ninjaA1,spvae2k5)
      acd11(5)=dotproduct(ninjaE3,spvak1e2)
      acd11(6)=dotproduct(ninjaA0,ninjaE3)
      acd11(7)=abb11(5)
      acd11(8)=dotproduct(ninjaA0,spvak1e2)
      acd11(9)=dotproduct(ninjaA0,spvae2k5)
      acd11(10)=abb11(12)
      acd11(11)=abb11(9)
      acd11(12)=dotproduct(ninjaE3,spvak1k2)
      acd11(13)=abb11(7)
      acd11(14)=dotproduct(ninjaE3,spvak2k5)
      acd11(15)=abb11(8)
      acd11(16)=acd11(2)*acd11(3)
      acd11(17)=acd11(1)*acd11(16)
      acd11(18)=acd11(5)*acd11(3)
      acd11(19)=acd11(4)*acd11(18)
      acd11(17)=acd11(17)+acd11(19)
      acd11(19)=acd11(7)*acd11(6)
      acd11(16)=acd11(8)*acd11(16)
      acd11(18)=acd11(9)*acd11(18)
      acd11(20)=acd11(10)*acd11(5)
      acd11(21)=acd11(11)*acd11(2)
      acd11(22)=acd11(13)*acd11(12)
      acd11(23)=acd11(15)*acd11(14)
      acd11(16)=acd11(23)+acd11(22)+acd11(21)+acd11(20)+acd11(18)+2.0_ki*acd11(&
      &19)+acd11(16)
      brack(ninjaidxt0x0mu0)=acd11(16)
      brack(ninjaidxt0x0mu2)=0.0_ki
      brack(ninjaidxt0x1mu0)=acd11(17)
      brack(ninjaidxt0x2mu0)=0.0_ki
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p3_part1part21_part25part25part1_d11h2_qp_ninja_t2")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p3_part1part21_part25part25part1_globalsl1_qp, only: epspow
      use p3_part1part21_part25part25part1_kinematics_qp
      use p3_part1part21_part25part25part1_abbrevd11h2_qp
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
end module     p3_part1part21_part25part25part1_d11h2l132_qp
