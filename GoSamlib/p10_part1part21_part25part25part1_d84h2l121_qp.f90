module     p10_part1part21_part25part25part1_d84h2l121_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p10_part1part &
   ! &21_part25part25part1/helicity2d84h2l121_qp.f90
   ! generator: buildfortran_n3.py
   use p10_part1part21_part25part25part1_config, only: ki => ki_qp
   use p10_part1part21_part25part25part1_util_qp, only: cond_t, d => metric_ten&
   &sor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, parameter :: ninjaidxt0 = 0
   public :: numerator_tmu
contains
!---#[ subroutine brack_0:
   pure subroutine brack_0(ninjaA0, brack)
      use p10_part1part21_part25part25part1_model_qp
      use p10_part1part21_part25part25part1_kinematics_qp
      use p10_part1part21_part25part25part1_color_qp
      use p10_part1part21_part25part25part1_abbrevd84h2_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0
      complex(ki), dimension(3) :: acd84
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd84(1)=dotproduct(ninjaA0,ninjaA0)
      acd84(2)=abb84(10)
      acd84(3)=acd84(2)*acd84(1)**2
      brack(ninjaidxt0)=acd84(3)
   end subroutine brack_0
!---#] subroutine brack_0:
!---#[ subroutine numerator_tmu:
   subroutine numerator_tmu(ncut, a, coeffs) &
   & bind(c, name="p10_part1part21_part25part25part1_d84h2_qp_ninja_tmu")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p10_part1part21_part25part25part1_globalsl1_qp, only: epspow
      use p10_part1part21_part25part25part1_kinematics_qp
      use p10_part1part21_part25part25part1_abbrevd84h2_qp
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki), dimension(0:3,0:*), intent(in) :: a
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: vecA0
	     vecA0(1:4) = + a(0:3,0)
      t1 = 0
      call cond_t(epspow.eq.t1,brack_0,vecA0, coeffs)
   end subroutine numerator_tmu
!---#] subroutine numerator_tmu:
end module     p10_part1part21_part25part25part1_d84h2l121_qp
