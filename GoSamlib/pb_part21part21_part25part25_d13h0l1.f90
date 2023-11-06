module     pb_part21part21_part25part25_d13h0l1
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/pb_part21part &
   ! &21_part25part25/helicity0d13h0l1.f90
   ! generator: buildfortran.py
   use pb_part21part21_part25part25_config, only: ki
   use pb_part21part21_part25part25_util, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_golem95
   public :: numerator_ninja
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use pb_part21part21_part25part25_model
      use pb_part21part21_part25part25_kinematics
      use pb_part21part21_part25part25_color
      use pb_part21part21_part25part25_abbrevd13h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc13(53)
      complex(ki) :: Qspvak1e1
      complex(ki) :: Qspvae1k1
      complex(ki) :: Qspvak1e2
      complex(ki) :: Qspvae2k1
      complex(ki) :: Qspval4e1
      complex(ki) :: Qspvae1l4
      complex(ki) :: Qspvak2e1
      complex(ki) :: Qspvae1k2
      complex(ki) :: Qspval3e1
      complex(ki) :: Qspvae1l3
      complex(ki) :: QspQ
      complex(ki) :: Qspvae1e2
      complex(ki) :: Qspvae2e1
      Qspvak1e1 = dotproduct(Q,spvak1e1)
      Qspvae1k1 = dotproduct(Q,spvae1k1)
      Qspvak1e2 = dotproduct(Q,spvak1e2)
      Qspvae2k1 = dotproduct(Q,spvae2k1)
      Qspval4e1 = dotproduct(Q,spval4e1)
      Qspvae1l4 = dotproduct(Q,spvae1l4)
      Qspvak2e1 = dotproduct(Q,spvak2e1)
      Qspvae1k2 = dotproduct(Q,spvae1k2)
      Qspval3e1 = dotproduct(Q,spval3e1)
      Qspvae1l3 = dotproduct(Q,spvae1l3)
      QspQ = dotproduct(Q,Q)
      Qspvae1e2 = dotproduct(Q,spvae1e2)
      Qspvae2e1 = dotproduct(Q,spvae2e1)
      acc13(1)=abb13(6)
      acc13(2)=abb13(7)
      acc13(3)=abb13(8)
      acc13(4)=abb13(10)
      acc13(5)=abb13(11)
      acc13(6)=abb13(12)
      acc13(7)=abb13(13)
      acc13(8)=abb13(14)
      acc13(9)=abb13(15)
      acc13(10)=abb13(16)
      acc13(11)=abb13(17)
      acc13(12)=abb13(18)
      acc13(13)=abb13(19)
      acc13(14)=abb13(20)
      acc13(15)=abb13(21)
      acc13(16)=abb13(22)
      acc13(17)=abb13(23)
      acc13(18)=abb13(24)
      acc13(19)=abb13(25)
      acc13(20)=abb13(26)
      acc13(21)=abb13(27)
      acc13(22)=abb13(28)
      acc13(23)=abb13(29)
      acc13(24)=abb13(30)
      acc13(25)=abb13(31)
      acc13(26)=abb13(32)
      acc13(27)=abb13(33)
      acc13(28)=abb13(34)
      acc13(29)=abb13(35)
      acc13(30)=abb13(36)
      acc13(31)=abb13(37)
      acc13(32)=abb13(38)
      acc13(33)=abb13(39)
      acc13(34)=abb13(40)
      acc13(35)=abb13(41)
      acc13(36)=abb13(42)
      acc13(37)=abb13(43)
      acc13(38)=abb13(44)
      acc13(39)=abb13(45)
      acc13(40)=abb13(46)
      acc13(41)=Qspvak1e1*acc13(29)
      acc13(42)=Qspvae1k1*acc13(2)
      acc13(43)=Qspvak1e2*acc13(20)
      acc13(44)=Qspvae2k1*acc13(7)
      acc13(45)=Qspval4e1*acc13(30)
      acc13(46)=Qspvae1l4*acc13(26)
      acc13(47)=Qspvak2e1*acc13(21)
      acc13(48)=Qspvae1k2*acc13(19)
      acc13(49)=Qspval3e1*acc13(36)
      acc13(50)=Qspvae1l3*acc13(32)
      acc13(51)=QspQ*acc13(3)
      acc13(41)=acc13(51)+acc13(50)+acc13(49)+acc13(48)+acc13(47)+acc13(46)+acc&
      &13(45)+acc13(44)+acc13(43)+acc13(42)+acc13(8)+acc13(41)
      acc13(41)=QspQ*acc13(41)
      acc13(42)=Qspvae1k1*acc13(35)
      acc13(43)=Qspvae2k1*acc13(24)
      acc13(44)=Qspvae1l4*acc13(37)
      acc13(45)=Qspvae1k2*acc13(5)
      acc13(42)=acc13(45)+acc13(44)+acc13(43)+acc13(6)+acc13(42)
      acc13(42)=Qspval3e1*acc13(42)
      acc13(43)=Qspvak1e1*acc13(33)
      acc13(44)=Qspvak1e2*acc13(39)
      acc13(45)=Qspval4e1*acc13(34)
      acc13(46)=Qspvak2e1*acc13(1)
      acc13(43)=acc13(46)+acc13(45)+acc13(44)+acc13(31)+acc13(43)
      acc13(43)=Qspvae1l3*acc13(43)
      acc13(44)=Qspvae1k1*acc13(18)
      acc13(45)=Qspvae2k1*acc13(38)
      acc13(46)=Qspvae1l4*acc13(22)
      acc13(44)=acc13(46)+acc13(45)+acc13(17)+acc13(44)
      acc13(44)=Qspvak2e1*acc13(44)
      acc13(45)=Qspvak1e1*acc13(13)
      acc13(46)=Qspvak1e2*acc13(40)
      acc13(47)=Qspval4e1*acc13(23)
      acc13(45)=acc13(47)+acc13(46)+acc13(15)+acc13(45)
      acc13(45)=Qspvae1k2*acc13(45)
      acc13(46)=Qspvae1e2*acc13(12)
      acc13(47)=Qspvae1e2*acc13(14)
      acc13(47)=acc13(10)+acc13(47)
      acc13(47)=Qspvae2e1*acc13(47)
      acc13(48)=Qspvak1e1*acc13(27)
      acc13(49)=Qspvae1k1*acc13(16)
      acc13(50)=Qspvak1e2*acc13(9)
      acc13(51)=Qspvae2k1*acc13(4)
      acc13(52)=Qspval4e1*acc13(28)
      acc13(53)=Qspvae1l4*acc13(25)
      brack=acc13(11)+acc13(41)+acc13(42)+acc13(43)+acc13(44)+acc13(45)+acc13(4&
      &6)+acc13(47)+acc13(48)+acc13(49)+acc13(50)+acc13(51)+acc13(52)+acc13(53)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use pb_part21part21_part25part25_globalsl1, only: epspow
      use pb_part21part21_part25part25_kinematics
      use pb_part21part21_part25part25_abbrevd13h0
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d13
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = -k3
      Q(:)  =cmplx(real(-Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d13 = 0.0_ki
      d13 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d13, ki), aimag(d13), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="pb_part21part21_part25part25_d13h0l1_ninja")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki_nin
      use pb_part21part21_part25part25_globalsl1, only: epspow
      use pb_part21part21_part25part25_kinematics
      use pb_part21part21_part25part25_abbrevd13h0
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d13
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = -k3
      Q(1:4)  =cmplx(real(-Q_ext(0:3)  -qshift(:),  ki_nin), aimag(-Q_ext(0:3))&
      &, ki)
      d13 = 0.0_ki
      d13 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d13, ki), aimag(d13), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module pb_part21part21_part25part25_d13h0l1
