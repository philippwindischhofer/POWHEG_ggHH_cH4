module     p10_part1part21_part25part25part1_d86h2l1_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p10_part1part &
   ! &21_part25part25part1/helicity2d86h2l1_qp.f90
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
      use p10_part1part21_part25part25part1_abbrevd86h2_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc86(41)
      complex(ki) :: QspQ
      complex(ki) :: Qspvak1k2
      complex(ki) :: Qspvak2e2
      complex(ki) :: Qspvak1l3
      complex(ki) :: Qspval3e2
      complex(ki) :: Qspvak1e2
      complex(ki) :: Qspvae2k5
      complex(ki) :: Qspvak2k5
      complex(ki) :: Qspvae2k2
      complex(ki) :: Qspval3k5
      complex(ki) :: Qspvae2l3
      complex(ki) :: Qspval4e2
      complex(ki) :: Qspvae2l4
      QspQ = dotproduct(Q,Q)
      Qspvak1k2 = dotproduct(Q,spvak1k2)
      Qspvak2e2 = dotproduct(Q,spvak2e2)
      Qspvak1l3 = dotproduct(Q,spvak1l3)
      Qspval3e2 = dotproduct(Q,spval3e2)
      Qspvak1e2 = dotproduct(Q,spvak1e2)
      Qspvae2k5 = dotproduct(Q,spvae2k5)
      Qspvak2k5 = dotproduct(Q,spvak2k5)
      Qspvae2k2 = dotproduct(Q,spvae2k2)
      Qspval3k5 = dotproduct(Q,spval3k5)
      Qspvae2l3 = dotproduct(Q,spvae2l3)
      Qspval4e2 = dotproduct(Q,spval4e2)
      Qspvae2l4 = dotproduct(Q,spvae2l4)
      acc86(1)=abb86(7)
      acc86(2)=abb86(8)
      acc86(3)=abb86(9)
      acc86(4)=abb86(11)
      acc86(5)=abb86(12)
      acc86(6)=abb86(13)
      acc86(7)=abb86(14)
      acc86(8)=abb86(15)
      acc86(9)=abb86(16)
      acc86(10)=abb86(17)
      acc86(11)=abb86(18)
      acc86(12)=abb86(19)
      acc86(13)=abb86(20)
      acc86(14)=abb86(21)
      acc86(15)=abb86(22)
      acc86(16)=abb86(23)
      acc86(17)=abb86(24)
      acc86(18)=abb86(25)
      acc86(19)=abb86(26)
      acc86(20)=abb86(27)
      acc86(21)=abb86(28)
      acc86(22)=abb86(29)
      acc86(23)=abb86(30)
      acc86(24)=abb86(31)
      acc86(25)=abb86(32)
      acc86(26)=abb86(33)
      acc86(27)=abb86(38)
      acc86(28)=abb86(40)
      acc86(29)=QspQ*acc86(19)
      acc86(30)=Qspvak1k2*acc86(12)
      acc86(30)=acc86(9)+acc86(30)
      acc86(30)=Qspvak2e2*acc86(30)
      acc86(31)=Qspvak1l3*acc86(1)
      acc86(31)=acc86(16)+acc86(31)
      acc86(31)=Qspval3e2*acc86(31)
      acc86(32)=QspQ*acc86(3)
      acc86(32)=acc86(23)+acc86(32)
      acc86(32)=Qspvak1e2*acc86(32)
      acc86(29)=acc86(32)+acc86(31)+acc86(30)+acc86(25)+acc86(29)
      acc86(29)=Qspvae2k5*acc86(29)
      acc86(30)=QspQ*acc86(17)
      acc86(31)=Qspvak2k5*acc86(12)
      acc86(31)=acc86(4)+acc86(31)
      acc86(31)=Qspvae2k2*acc86(31)
      acc86(32)=Qspval3k5*acc86(1)
      acc86(32)=acc86(10)+acc86(32)
      acc86(32)=Qspvae2l3*acc86(32)
      acc86(30)=acc86(32)+acc86(31)+acc86(20)+acc86(30)
      acc86(30)=Qspvak1e2*acc86(30)
      acc86(31)=acc86(27)*Qspval4e2
      acc86(32)=acc86(7)*Qspvae2l4
      acc86(33)=Qspvak1k2*acc86(24)
      acc86(34)=Qspvak1l3*acc86(21)
      acc86(35)=Qspvak2k5*acc86(5)
      acc86(36)=Qspval3k5*acc86(2)
      acc86(37)=QspQ*acc86(18)
      acc86(38)=Qspvak1k2*acc86(26)
      acc86(38)=acc86(8)+acc86(38)
      acc86(38)=Qspvak2e2*acc86(38)
      acc86(39)=Qspvak2k5*acc86(15)
      acc86(39)=acc86(14)+acc86(39)
      acc86(39)=Qspvae2k2*acc86(39)
      acc86(40)=Qspvak1l3*acc86(22)
      acc86(40)=acc86(28)+acc86(40)
      acc86(40)=Qspval3e2*acc86(40)
      acc86(41)=Qspval3k5*acc86(11)
      acc86(41)=acc86(13)+acc86(41)
      acc86(41)=Qspvae2l3*acc86(41)
      brack=acc86(6)+acc86(29)+acc86(30)+acc86(31)+acc86(32)+acc86(33)+acc86(34&
      &)+acc86(35)+acc86(36)+acc86(37)+acc86(38)+acc86(39)+acc86(40)+acc86(41)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p10_part1part21_part25part25part1_d86h2l1_qp_ninja")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki_nin
      use p10_part1part21_part25part25part1_globalsl1_qp, only: epspow
      use p10_part1part21_part25part25part1_kinematics_qp
      use p10_part1part21_part25part25part1_abbrevd86h2_qp
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d86
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      Q(1:4)  =cmplx(real(-Q_ext(0:3),  ki_nin), aimag(-Q_ext(0:3)), ki)
      d86 = 0.0_ki
      d86 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d86, ki), aimag(d86), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p10_part1part21_part25part25part1_d86h2l1_qp
