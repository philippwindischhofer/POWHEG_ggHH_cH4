module     p3_part1part21_part25part25part1_d11h0l121_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p3_part1part2 &
   ! &1_part25part25part1/helicity0d11h0l121_qp.f90
   ! generator: buildfortran_n3.py
   use p3_part1part21_part25part25part1_config, only: ki => ki_qp
   use p3_part1part21_part25part25part1_util_qp, only: cond_t, d => metric_tens&
   &or
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, parameter :: ninjaidxt0 = 0
   public :: numerator_tmu
contains
!---#[ subroutine numerator_tmu:
   subroutine numerator_tmu(ncut, a, coeffs) &
   & bind(c, name="p3_part1part21_part25part25part1_d11h0_qp_ninja_tmu")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p3_part1part21_part25part25part1_globalsl1_qp, only: epspow
      use p3_part1part21_part25part25part1_kinematics_qp
      use p3_part1part21_part25part25part1_abbrevd11h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki), dimension(0:3,0:*), intent(in) :: a
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      complex(ki), dimension(4) :: vecA0
      qshift = -k2
	     vecA0(1:4) = - a(0:3,0)
      t1 = 0
      coeffs(0) = 0.0_ki
   end subroutine numerator_tmu
!---#] subroutine numerator_tmu:
end module     p3_part1part21_part25part25part1_d11h0l121_qp
