module     p10_part1part21_part25part25part1_d11h2l131_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p10_part1part &
   ! &21_part25part25part1/helicity2d11h2l131_qp.f90
   ! generator: buildfortran_tn3.py
   use p10_part1part21_part25part25part1_config, only: ki => ki_qp
   use p10_part1part21_part25part25part1_util_qp, only: cond_t, d => metric_ten&
   &sor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, parameter :: ninjaidxt3mu0 = 0
   integer, parameter :: ninjaidxt2mu0 = 1
   integer, parameter :: ninjaidxt1mu0 = 2
   integer, parameter :: ninjaidxt1mu2 = 3
   integer, parameter :: ninjaidxt0mu0 = 4
   integer, parameter :: ninjaidxt0mu2 = 5
   public :: numerator_t3
contains
!---#[ subroutine brack_31:
   pure subroutine brack_31(ninjaA, ninjaE3, ninjaE4, ninjaP, brack)
      use p10_part1part21_part25part25part1_model_qp
      use p10_part1part21_part25part25part1_kinematics_qp
      use p10_part1part21_part25part25part1_color_qp
      use p10_part1part21_part25part25part1_abbrevd11h2_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(4) :: acd11
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd11(1)=dotproduct(ninjaE3,spvak1e2)
      acd11(2)=dotproduct(ninjaE3,spvae2k5)
      acd11(3)=abb11(4)
      acd11(4)=acd11(3)*acd11(2)*acd11(1)
      brack(ninjaidxt3mu0)=0.0_ki
      brack(ninjaidxt2mu0)=acd11(4)
   end subroutine brack_31
!---#] subroutine brack_31:
!---#[ subroutine brack_32:
   pure subroutine brack_32(ninjaA, ninjaE3, ninjaE4, ninjaP, brack)
      use p10_part1part21_part25part25part1_model_qp
      use p10_part1part21_part25part25part1_kinematics_qp
      use p10_part1part21_part25part25part1_color_qp
      use p10_part1part21_part25part25part1_abbrevd11h2_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(27) :: acd11
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd11(1)=dotproduct(ninjaE3,spvak1e2)
      acd11(2)=dotproduct(ninjaE4,spvae2k5)
      acd11(3)=abb11(4)
      acd11(4)=dotproduct(ninjaE3,spvae2k5)
      acd11(5)=dotproduct(ninjaE4,spvak1e2)
      acd11(6)=abb11(5)
      acd11(7)=dotproduct(ninjaA,ninjaE3)
      acd11(8)=dotproduct(ninjaA,spvak1e2)
      acd11(9)=dotproduct(ninjaA,spvae2k5)
      acd11(10)=abb11(12)
      acd11(11)=abb11(9)
      acd11(12)=dotproduct(ninjaE3,spvak1k2)
      acd11(13)=abb11(7)
      acd11(14)=dotproduct(ninjaE3,spvak2k5)
      acd11(15)=abb11(8)
      acd11(16)=dotproduct(ninjaA,ninjaA)
      acd11(17)=dotproduct(ninjaA,spvak1k2)
      acd11(18)=dotproduct(ninjaA,spvak2k5)
      acd11(19)=abb11(6)
      acd11(20)=acd11(4)*acd11(3)
      acd11(21)=acd11(20)*acd11(5)
      acd11(22)=acd11(1)*acd11(3)
      acd11(23)=acd11(22)*acd11(2)
      acd11(21)=acd11(6)+acd11(21)+acd11(23)
      acd11(20)=acd11(8)*acd11(20)
      acd11(22)=acd11(9)*acd11(22)
      acd11(23)=acd11(10)*acd11(1)
      acd11(24)=acd11(11)*acd11(4)
      acd11(25)=acd11(7)*acd11(6)
      acd11(26)=acd11(12)*acd11(13)
      acd11(27)=acd11(14)*acd11(15)
      acd11(20)=acd11(27)+acd11(26)+2.0_ki*acd11(25)+acd11(24)+acd11(23)+acd11(&
      &20)+acd11(22)
      acd11(22)=ninjaP*acd11(21)
      acd11(23)=acd11(9)*acd11(3)
      acd11(23)=acd11(10)+acd11(23)
      acd11(23)=acd11(8)*acd11(23)
      acd11(24)=acd11(11)*acd11(9)
      acd11(25)=acd11(16)*acd11(6)
      acd11(26)=acd11(17)*acd11(13)
      acd11(27)=acd11(18)*acd11(15)
      acd11(22)=acd11(19)+acd11(27)+acd11(26)+acd11(25)+acd11(24)+acd11(22)+acd&
      &11(23)
      brack(ninjaidxt1mu0)=acd11(20)
      brack(ninjaidxt1mu2)=0.0_ki
      brack(ninjaidxt0mu0)=acd11(22)
      brack(ninjaidxt0mu2)=acd11(21)
   end subroutine brack_32
!---#] subroutine brack_32:
!---#[ subroutine numerator_t3:
   subroutine numerator_t3(ncut, a, b, c, param, deg, coeffs) &
   & bind(c, name="p10_part1part21_part25part25part1_d11h2_qp_ninja_t3")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p10_part1part21_part25part25part1_globalsl1_qp, only: epspow
      use p10_part1part21_part25part25part1_kinematics_qp
      use p10_part1part21_part25part25part1_abbrevd11h2_qp
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a, b, c
      complex(ki), intent(in) :: param
      complex(ki), dimension(4) :: vecA, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      qshift = -k2
	     vecA(1:4) = - a(0:3) - qshift(1:4)
	     vecB(1:4) = - b(0:3)
	     vecC(1:4) = - c(0:3)
      if (deg.lt.0) return
      t1 = 0
      call cond_t(epspow.eq.t1,brack_31,vecA,vecB,vecC,param,coeffs)
      if (deg.le.(1+(0))) return
      call cond_t(epspow.eq.t1,brack_32,vecA,vecB,vecC,param,coeffs)
   end subroutine numerator_t3
!---#] subroutine numerator_t3:
end module     p10_part1part21_part25part25part1_d11h2l131_qp
