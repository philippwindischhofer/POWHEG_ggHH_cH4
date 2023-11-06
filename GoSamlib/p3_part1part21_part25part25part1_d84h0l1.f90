module     p3_part1part21_part25part25part1_d84h0l1
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p3_part1part2 &
   ! &1_part25part25part1/helicity0d84h0l1.f90
   ! generator: buildfortran.py
   use p3_part1part21_part25part25part1_config, only: ki
   use p3_part1part21_part25part25part1_util, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_ninja
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p3_part1part21_part25part25part1_model
      use p3_part1part21_part25part25part1_kinematics
      use p3_part1part21_part25part25part1_color
      use p3_part1part21_part25part25part1_abbrevd84h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc84(44)
      complex(ki) :: Qspvak2k1
      complex(ki) :: Qspvak5k2
      complex(ki) :: Qspval4e2
      complex(ki) :: Qspvae2l4
      complex(ki) :: Qspval3e2
      complex(ki) :: Qspvae2l3
      complex(ki) :: Qspvak2e2
      complex(ki) :: Qspvae2k2
      complex(ki) :: QspQ
      complex(ki) :: Qspvae2k1
      complex(ki) :: Qspvak5e2
      Qspvak2k1 = dotproduct(Q,spvak2k1)
      Qspvak5k2 = dotproduct(Q,spvak5k2)
      Qspval4e2 = dotproduct(Q,spval4e2)
      Qspvae2l4 = dotproduct(Q,spvae2l4)
      Qspval3e2 = dotproduct(Q,spval3e2)
      Qspvae2l3 = dotproduct(Q,spvae2l3)
      Qspvak2e2 = dotproduct(Q,spvak2e2)
      Qspvae2k2 = dotproduct(Q,spvae2k2)
      QspQ = dotproduct(Q,Q)
      Qspvae2k1 = dotproduct(Q,spvae2k1)
      Qspvak5e2 = dotproduct(Q,spvak5e2)
      acc84(1)=abb84(7)
      acc84(2)=abb84(8)
      acc84(3)=abb84(9)
      acc84(4)=abb84(10)
      acc84(5)=abb84(12)
      acc84(6)=abb84(13)
      acc84(7)=abb84(14)
      acc84(8)=abb84(15)
      acc84(9)=abb84(16)
      acc84(10)=abb84(17)
      acc84(11)=abb84(18)
      acc84(12)=abb84(19)
      acc84(13)=abb84(20)
      acc84(14)=abb84(21)
      acc84(15)=abb84(22)
      acc84(16)=abb84(23)
      acc84(17)=abb84(24)
      acc84(18)=abb84(25)
      acc84(19)=abb84(26)
      acc84(20)=abb84(30)
      acc84(21)=abb84(31)
      acc84(22)=abb84(32)
      acc84(23)=abb84(33)
      acc84(24)=abb84(34)
      acc84(25)=abb84(35)
      acc84(26)=abb84(36)
      acc84(27)=abb84(37)
      acc84(28)=abb84(39)
      acc84(29)=abb84(40)
      acc84(30)=abb84(41)
      acc84(31)=abb84(43)
      acc84(32)=abb84(44)
      acc84(33)=abb84(46)
      acc84(34)=Qspvak2k1*acc84(18)
      acc84(35)=Qspvak5k2*acc84(6)
      acc84(36)=Qspval4e2*acc84(24)
      acc84(37)=-Qspvae2l4*acc84(21)
      acc84(38)=Qspval3e2*acc84(29)
      acc84(39)=Qspvae2l3*acc84(27)
      acc84(40)=Qspvak2e2*acc84(13)
      acc84(41)=Qspvae2k2*acc84(31)
      acc84(42)=QspQ*acc84(4)
      acc84(34)=acc84(42)+acc84(41)+acc84(40)+acc84(39)+acc84(38)+acc84(37)+acc&
      &84(36)+acc84(35)+acc84(17)+acc84(34)
      acc84(34)=QspQ*acc84(34)
      acc84(35)=Qspvak2k1*acc84(19)
      acc84(36)=Qspval4e2*acc84(33)
      acc84(37)=Qspval3e2*acc84(11)
      acc84(38)=Qspvak2e2*acc84(1)
      acc84(35)=acc84(38)+acc84(37)+acc84(36)+acc84(12)+acc84(35)
      acc84(35)=Qspvae2k2*acc84(35)
      acc84(36)=Qspvak5k2*acc84(16)
      acc84(37)=Qspvae2l4*acc84(32)
      acc84(38)=Qspvae2l3*acc84(25)
      acc84(36)=acc84(38)+acc84(37)+acc84(8)+acc84(36)
      acc84(36)=Qspvak2e2*acc84(36)
      acc84(37)=Qspvak5k2*acc84(9)
      acc84(38)=Qspvae2l4*acc84(30)
      acc84(37)=acc84(38)+acc84(28)+acc84(37)
      acc84(37)=Qspval3e2*acc84(37)
      acc84(38)=Qspvak2k1*acc84(15)
      acc84(39)=Qspval4e2*acc84(22)
      acc84(38)=acc84(39)+acc84(26)+acc84(38)
      acc84(38)=Qspvae2l3*acc84(38)
      acc84(39)=Qspvae2k1*acc84(3)
      acc84(40)=Qspvae2k1*acc84(2)
      acc84(40)=acc84(10)+acc84(40)
      acc84(40)=Qspvak5e2*acc84(40)
      acc84(41)=Qspvak2k1*acc84(14)
      acc84(42)=Qspvak5k2*acc84(5)
      acc84(43)=Qspval4e2*acc84(23)
      acc84(44)=Qspvae2l4*acc84(20)
      brack=acc84(7)+acc84(34)+acc84(35)+acc84(36)+acc84(37)+acc84(38)+acc84(39&
      &)+acc84(40)+acc84(41)+acc84(42)+acc84(43)+acc84(44)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p3_part1part21_part25part25part1_d84h0l1_ninja")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki_nin
      use p3_part1part21_part25part25part1_globalsl1, only: epspow
      use p3_part1part21_part25part25part1_kinematics
      use p3_part1part21_part25part25part1_abbrevd84h0
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d84
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      Q(1:4)  =cmplx(real(+Q_ext(0:3),  ki_nin), aimag(+Q_ext(0:3)), ki)
      d84 = 0.0_ki
      d84 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d84, ki), aimag(d84), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p3_part1part21_part25part25part1_d84h0l1
