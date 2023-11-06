module     p10_part1part21_part25part25part1_d88h0l1
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p10_part1part &
   ! &21_part25part25part1/helicity0d88h0l1.f90
   ! generator: buildfortran.py
   use p10_part1part21_part25part25part1_config, only: ki
   use p10_part1part21_part25part25part1_util, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_ninja
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p10_part1part21_part25part25part1_model
      use p10_part1part21_part25part25part1_kinematics
      use p10_part1part21_part25part25part1_color
      use p10_part1part21_part25part25part1_abbrevd88h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc88(41)
      complex(ki) :: QspQ
      complex(ki) :: Qspvak2k1
      complex(ki) :: Qspvae2k2
      complex(ki) :: Qspval4k1
      complex(ki) :: Qspvae2l4
      complex(ki) :: Qspvae2k1
      complex(ki) :: Qspvak5e2
      complex(ki) :: Qspvak5k2
      complex(ki) :: Qspvak2e2
      complex(ki) :: Qspvak5l4
      complex(ki) :: Qspval4e2
      complex(ki) :: Qspvae2l3
      complex(ki) :: Qspval3e2
      QspQ = dotproduct(Q,Q)
      Qspvak2k1 = dotproduct(Q,spvak2k1)
      Qspvae2k2 = dotproduct(Q,spvae2k2)
      Qspval4k1 = dotproduct(Q,spval4k1)
      Qspvae2l4 = dotproduct(Q,spvae2l4)
      Qspvae2k1 = dotproduct(Q,spvae2k1)
      Qspvak5e2 = dotproduct(Q,spvak5e2)
      Qspvak5k2 = dotproduct(Q,spvak5k2)
      Qspvak2e2 = dotproduct(Q,spvak2e2)
      Qspvak5l4 = dotproduct(Q,spvak5l4)
      Qspval4e2 = dotproduct(Q,spval4e2)
      Qspvae2l3 = dotproduct(Q,spvae2l3)
      Qspval3e2 = dotproduct(Q,spval3e2)
      acc88(1)=abb88(7)
      acc88(2)=abb88(8)
      acc88(3)=abb88(9)
      acc88(4)=abb88(11)
      acc88(5)=abb88(12)
      acc88(6)=abb88(13)
      acc88(7)=abb88(14)
      acc88(8)=abb88(15)
      acc88(9)=abb88(16)
      acc88(10)=abb88(17)
      acc88(11)=abb88(18)
      acc88(12)=abb88(19)
      acc88(13)=abb88(20)
      acc88(14)=abb88(21)
      acc88(15)=abb88(22)
      acc88(16)=abb88(23)
      acc88(17)=abb88(24)
      acc88(18)=abb88(25)
      acc88(19)=abb88(26)
      acc88(20)=abb88(27)
      acc88(21)=abb88(28)
      acc88(22)=abb88(29)
      acc88(23)=abb88(30)
      acc88(24)=abb88(31)
      acc88(25)=abb88(32)
      acc88(26)=abb88(33)
      acc88(27)=abb88(34)
      acc88(28)=abb88(36)
      acc88(29)=QspQ*acc88(18)
      acc88(30)=Qspvak2k1*acc88(11)
      acc88(30)=acc88(26)+acc88(30)
      acc88(30)=Qspvae2k2*acc88(30)
      acc88(31)=Qspval4k1*acc88(1)
      acc88(31)=acc88(7)+acc88(31)
      acc88(31)=Qspvae2l4*acc88(31)
      acc88(32)=QspQ*acc88(3)
      acc88(32)=acc88(23)+acc88(32)
      acc88(32)=Qspvae2k1*acc88(32)
      acc88(29)=acc88(32)+acc88(31)+acc88(30)+acc88(24)+acc88(29)
      acc88(29)=Qspvak5e2*acc88(29)
      acc88(30)=QspQ*acc88(16)
      acc88(31)=Qspvak5k2*acc88(11)
      acc88(31)=acc88(9)+acc88(31)
      acc88(31)=Qspvak2e2*acc88(31)
      acc88(32)=Qspvak5l4*acc88(1)
      acc88(32)=acc88(13)+acc88(32)
      acc88(32)=Qspval4e2*acc88(32)
      acc88(30)=acc88(32)+acc88(31)+acc88(19)+acc88(30)
      acc88(30)=Qspvae2k1*acc88(30)
      acc88(31)=acc88(28)*Qspvae2l3
      acc88(32)=acc88(12)*Qspval3e2
      acc88(33)=Qspvak2k1*acc88(27)
      acc88(34)=Qspval4k1*acc88(20)
      acc88(35)=Qspvak5k2*acc88(10)
      acc88(36)=Qspvak5l4*acc88(2)
      acc88(37)=QspQ*acc88(17)
      acc88(38)=Qspvak5k2*acc88(15)
      acc88(38)=acc88(6)+acc88(38)
      acc88(38)=Qspvak2e2*acc88(38)
      acc88(39)=Qspvak2k1*acc88(25)
      acc88(39)=acc88(5)+acc88(39)
      acc88(39)=Qspvae2k2*acc88(39)
      acc88(40)=Qspvak5l4*acc88(14)
      acc88(40)=acc88(22)+acc88(40)
      acc88(40)=Qspval4e2*acc88(40)
      acc88(41)=Qspval4k1*acc88(21)
      acc88(41)=acc88(4)+acc88(41)
      acc88(41)=Qspvae2l4*acc88(41)
      brack=acc88(8)+acc88(29)+acc88(30)+acc88(31)+acc88(32)+acc88(33)+acc88(34&
      &)+acc88(35)+acc88(36)+acc88(37)+acc88(38)+acc88(39)+acc88(40)+acc88(41)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p10_part1part21_part25part25part1_d88h0l1_ninja")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki_nin
      use p10_part1part21_part25part25part1_globalsl1, only: epspow
      use p10_part1part21_part25part25part1_kinematics
      use p10_part1part21_part25part25part1_abbrevd88h0
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d88
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      Q(1:4)  =cmplx(real(-Q_ext(0:3),  ki_nin), aimag(-Q_ext(0:3)), ki)
      d88 = 0.0_ki
      d88 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d88, ki), aimag(d88), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p10_part1part21_part25part25part1_d88h0l1
