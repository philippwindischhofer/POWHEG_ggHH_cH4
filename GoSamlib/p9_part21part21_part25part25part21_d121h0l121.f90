module     p9_part21part21_part25part25part21_d121h0l121
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d121h0l121.f90
   ! generator: buildfortran_n3.py
   use p9_part21part21_part25part25part21_config, only: ki
   use p9_part21part21_part25part25part21_util, only: cond_t, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, parameter :: ninjaidxt0 = 0
   public :: numerator_tmu
contains
!---#[ subroutine brack_0:
   pure subroutine brack_0(ninjaA0, brack)
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd121h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0
      complex(ki), dimension(3) :: acd121
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd121(1)=dotproduct(ninjaA0,ninjaA0)
      acd121(2)=abb121(9)
      acd121(3)=acd121(2)*acd121(1)**2
      brack(ninjaidxt0)=acd121(3)
   end subroutine brack_0
!---#] subroutine brack_0:
!---#[ subroutine numerator_tmu:
   subroutine numerator_tmu(ncut, a, coeffs) &
   & bind(c, name="p9_part21part21_part25part25part21_d121h0_ninja_tmu")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki => ki_nin
      use p9_part21part21_part25part25part21_globalsl1, only: epspow
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_abbrevd121h0
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
end module     p9_part21part21_part25part25part21_d121h0l121
