module     p2_part21part21_part25part25part21_d202h0l121_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d202h0l121_qp.f90
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
      use p2_part21part21_part25part25part21_abbrevd202h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0
      complex(ki), dimension(10) :: acd202
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd202(1)=dotproduct(ninjaA0,ninjaA0)
      acd202(2)=dotproduct(ninjaA0,spvae1e2)
      acd202(3)=dotproduct(ninjaA0,spvae5e1)
      acd202(4)=dotproduct(ninjaA0,spvae2e5)
      acd202(5)=abb202(81)
      acd202(6)=dotproduct(ninjaA0,spvae2e1)
      acd202(7)=dotproduct(ninjaA0,spvae1e5)
      acd202(8)=dotproduct(ninjaA0,spvae5e2)
      acd202(9)=-acd202(4)*acd202(3)*acd202(2)
      acd202(10)=-acd202(8)*acd202(7)*acd202(6)
      acd202(9)=acd202(9)+acd202(10)
      acd202(9)=acd202(9)*acd202(5)*acd202(1)
      brack(ninjaidxt0)=acd202(9)
   end subroutine brack_0
!---#] subroutine brack_0:
!---#[ subroutine numerator_tmu:
   subroutine numerator_tmu(ncut, a, coeffs) &
   & bind(c, name="p2_part21part21_part25part25part21_d202h0_qp_ninja_tmu")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd202h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki), dimension(0:3,0:*), intent(in) :: a
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      complex(ki), dimension(4) :: vecA0
      qshift = -k3-k4
	     vecA0(1:4) = - a(0:3,0)
      t1 = 0
      call cond_t(epspow.eq.t1,brack_0,vecA0, coeffs)
   end subroutine numerator_tmu
!---#] subroutine numerator_tmu:
end module     p2_part21part21_part25part25part21_d202h0l121_qp
