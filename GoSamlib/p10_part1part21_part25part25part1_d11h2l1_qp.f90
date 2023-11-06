module     p10_part1part21_part25part25part1_d11h2l1_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p10_part1part &
   ! &21_part25part25part1/helicity2d11h2l1_qp.f90
   ! generator: buildfortran.py
   use p10_part1part21_part25part25part1_config, only: ki => ki_qp
   use p10_part1part21_part25part25part1_util_qp, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_ninja
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p10_part1part21_part25part25part1_model_qp
      use p10_part1part21_part25part25part1_kinematics_qp
      use p10_part1part21_part25part25part1_color_qp
      use p10_part1part21_part25part25part1_abbrevd11h2_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc11(12)
      complex(ki) :: Qspvae2k5
      complex(ki) :: Qspvak1e2
      complex(ki) :: Qspvak2k5
      complex(ki) :: Qspvak1k2
      complex(ki) :: QspQ
      Qspvae2k5 = dotproduct(Q,spvae2k5)
      Qspvak1e2 = dotproduct(Q,spvak1e2)
      Qspvak2k5 = dotproduct(Q,spvak2k5)
      Qspvak1k2 = dotproduct(Q,spvak1k2)
      QspQ = dotproduct(Q,Q)
      acc11(1)=abb11(4)
      acc11(2)=abb11(5)
      acc11(3)=abb11(6)
      acc11(4)=abb11(7)
      acc11(5)=abb11(8)
      acc11(6)=abb11(9)
      acc11(7)=abb11(12)
      acc11(8)=acc11(1)*Qspvae2k5
      acc11(8)=acc11(7)+acc11(8)
      acc11(8)=Qspvak1e2*acc11(8)
      acc11(9)=acc11(6)*Qspvae2k5
      acc11(10)=Qspvak2k5*acc11(5)
      acc11(11)=Qspvak1k2*acc11(4)
      acc11(12)=QspQ*acc11(2)
      brack=acc11(3)+acc11(8)+acc11(9)+acc11(10)+acc11(11)+acc11(12)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p10_part1part21_part25part25part1_d11h2l1_qp_ninja")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki_nin
      use p10_part1part21_part25part25part1_globalsl1_qp, only: epspow
      use p10_part1part21_part25part25part1_kinematics_qp
      use p10_part1part21_part25part25part1_abbrevd11h2_qp
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d11
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = -k2
      Q(1:4)  =cmplx(real(-Q_ext(0:3)  -qshift(:),  ki_nin), aimag(-Q_ext(0:3))&
      &, ki)
      d11 = 0.0_ki
      d11 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d11, ki), aimag(d11), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p10_part1part21_part25part25part1_d11h2l1_qp
