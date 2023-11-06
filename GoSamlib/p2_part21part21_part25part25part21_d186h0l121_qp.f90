module     p2_part21part21_part25part25part21_d186h0l121_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d186h0l121_qp.f90
   ! generator: buildfortran_n3.py
   use p2_part21part21_part25part25part21_config, only: ki => ki_qp
   use p2_part21part21_part25part25part21_util_qp, only: cond_t, d => metric_te&
   &nsor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, parameter :: ninjaidxt0 = 0
   public :: numerator_tmu
contains
!---#[ subroutine brack_0:
   pure subroutine brack_0(ninjaA0, brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd186h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0
      complex(ki), dimension(10) :: acd186
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd186(1)=dotproduct(ninjaA0,ninjaA0)
      acd186(2)=dotproduct(ninjaA0,spvae1e2)
      acd186(3)=dotproduct(ninjaA0,spvae5e1)
      acd186(4)=dotproduct(ninjaA0,spvae2e5)
      acd186(5)=abb186(90)
      acd186(6)=dotproduct(ninjaA0,spvae2e1)
      acd186(7)=dotproduct(ninjaA0,spvae1e5)
      acd186(8)=dotproduct(ninjaA0,spvae5e2)
      acd186(9)=acd186(4)*acd186(3)*acd186(2)
      acd186(10)=acd186(8)*acd186(7)*acd186(6)
      acd186(9)=acd186(9)+acd186(10)
      acd186(9)=acd186(9)*acd186(5)*acd186(1)
      brack(ninjaidxt0)=acd186(9)
   end subroutine brack_0
!---#] subroutine brack_0:
!---#[ subroutine numerator_tmu:
   subroutine numerator_tmu(ncut, a, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d186h0_qp_ninja_tmu")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd186h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki), dimension(0:3,0:*), intent(in) :: a
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      complex(ki), dimension(4) :: vecA0
      qshift = k2
	     vecA0(1:4) = + a(0:3,0)
      t1 = 0
      call cond_t(epspow.eq.t1,brack_0,vecA0, coeffs)
   end subroutine numerator_tmu
!---#] subroutine numerator_tmu:
end module     p2_part21part21_part25part25part21_d186h0l121_qp
