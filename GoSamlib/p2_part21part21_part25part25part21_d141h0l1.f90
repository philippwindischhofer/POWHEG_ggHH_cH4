module     p2_part21part21_part25part25part21_d141h0l1
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d141h0l1.f90
   ! generator: buildfortran.py
   use p2_part21part21_part25part25part21_config, only: ki
   use p2_part21part21_part25part25part21_util, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_ninja
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_color
      use p2_part21part21_part25part25part21_abbrevd141h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc141(55)
      complex(ki) :: QspQ
      complex(ki) :: Qspk2
      complex(ki) :: Qspvak2e1
      complex(ki) :: Qspvak2k5
      complex(ki) :: Qspvak2l4
      complex(ki) :: Qspval4e1
      complex(ki) :: Qspvak5e1
      complex(ki) :: Qspvae1k2
      complex(ki) :: Qspvak5k2
      complex(ki) :: Qspval4k2
      complex(ki) :: Qspvae1l4
      complex(ki) :: Qspvae1k5
      complex(ki) :: Qspvak5l4
      complex(ki) :: Qspval4k5
      complex(ki) :: Qspvae1l3
      complex(ki) :: Qspval3e1
      QspQ = dotproduct(Q,Q)
      Qspk2 = dotproduct(Q,k2)
      Qspvak2e1 = dotproduct(Q,spvak2e1)
      Qspvak2k5 = dotproduct(Q,spvak2k5)
      Qspvak2l4 = dotproduct(Q,spvak2l4)
      Qspval4e1 = dotproduct(Q,spval4e1)
      Qspvak5e1 = dotproduct(Q,spvak5e1)
      Qspvae1k2 = dotproduct(Q,spvae1k2)
      Qspvak5k2 = dotproduct(Q,spvak5k2)
      Qspval4k2 = dotproduct(Q,spval4k2)
      Qspvae1l4 = dotproduct(Q,spvae1l4)
      Qspvae1k5 = dotproduct(Q,spvae1k5)
      Qspvak5l4 = dotproduct(Q,spvak5l4)
      Qspval4k5 = dotproduct(Q,spval4k5)
      Qspvae1l3 = dotproduct(Q,spvae1l3)
      Qspval3e1 = dotproduct(Q,spval3e1)
      acc141(1)=abb141(7)
      acc141(2)=abb141(8)
      acc141(3)=abb141(9)
      acc141(4)=abb141(10)
      acc141(5)=abb141(11)
      acc141(6)=abb141(12)
      acc141(7)=abb141(14)
      acc141(8)=abb141(15)
      acc141(9)=abb141(16)
      acc141(10)=abb141(17)
      acc141(11)=abb141(18)
      acc141(12)=abb141(19)
      acc141(13)=abb141(20)
      acc141(14)=abb141(21)
      acc141(15)=abb141(22)
      acc141(16)=abb141(23)
      acc141(17)=abb141(24)
      acc141(18)=abb141(25)
      acc141(19)=abb141(26)
      acc141(20)=abb141(27)
      acc141(21)=abb141(28)
      acc141(22)=abb141(29)
      acc141(23)=abb141(30)
      acc141(24)=abb141(31)
      acc141(25)=abb141(32)
      acc141(26)=abb141(33)
      acc141(27)=abb141(34)
      acc141(28)=abb141(35)
      acc141(29)=abb141(36)
      acc141(30)=abb141(37)
      acc141(31)=abb141(38)
      acc141(32)=abb141(39)
      acc141(33)=abb141(41)
      acc141(34)=abb141(43)
      acc141(35)=abb141(44)
      acc141(36)=abb141(45)
      acc141(37)=abb141(46)
      acc141(38)=abb141(47)
      acc141(39)=abb141(49)
      acc141(40)=abb141(53)
      acc141(41)=QspQ*acc141(5)
      acc141(42)=Qspk2*acc141(3)
      acc141(42)=acc141(41)+acc141(20)+acc141(42)
      acc141(42)=Qspvak2e1*acc141(42)
      acc141(43)=Qspk2*acc141(40)
      acc141(44)=Qspvak2k5*acc141(13)
      acc141(45)=QspQ*acc141(19)
      acc141(46)=acc141(2)*Qspvak2l4
      acc141(46)=acc141(30)+acc141(46)
      acc141(46)=Qspval4e1*acc141(46)
      acc141(47)=-Qspvak2k5*acc141(1)
      acc141(47)=acc141(10)+acc141(47)
      acc141(47)=Qspvak5e1*acc141(47)
      acc141(42)=acc141(42)+acc141(47)+acc141(46)+acc141(45)+acc141(44)+acc141(&
      &21)+acc141(43)
      acc141(42)=Qspvae1k2*acc141(42)
      acc141(43)=Qspk2*acc141(23)
      acc141(44)=Qspvak5k2*acc141(4)
      acc141(45)=QspQ*acc141(22)
      acc141(46)=acc141(2)*Qspval4k2
      acc141(46)=acc141(32)+acc141(46)
      acc141(46)=Qspvae1l4*acc141(46)
      acc141(47)=-Qspvak5k2*acc141(1)
      acc141(47)=acc141(7)+acc141(47)
      acc141(47)=Qspvae1k5*acc141(47)
      acc141(43)=acc141(47)+acc141(46)+acc141(45)+acc141(44)+acc141(18)+acc141(&
      &43)
      acc141(43)=Qspvak2e1*acc141(43)
      acc141(44)=QspQ*acc141(14)
      acc141(45)=acc141(2)*Qspvak5l4
      acc141(45)=acc141(25)+acc141(45)
      acc141(45)=Qspval4e1*acc141(45)
      acc141(41)=acc141(15)+acc141(41)
      acc141(41)=Qspvak5e1*acc141(41)
      acc141(41)=acc141(41)+acc141(45)+acc141(34)+acc141(44)
      acc141(41)=Qspvae1k5*acc141(41)
      acc141(44)=Qspvak2l4*acc141(31)
      acc141(45)=Qspvak5l4*acc141(35)
      acc141(44)=acc141(45)+acc141(39)+acc141(44)
      acc141(44)=Qspval4e1*acc141(44)
      acc141(45)=Qspval4k2*acc141(38)
      acc141(46)=Qspval4k5*acc141(27)
      acc141(45)=acc141(46)+acc141(24)+acc141(45)
      acc141(45)=Qspvae1l4*acc141(45)
      acc141(46)=QspQ*acc141(12)
      acc141(47)=acc141(2)*Qspval4k5
      acc141(47)=acc141(28)+acc141(47)
      acc141(47)=Qspvae1l4*acc141(47)
      acc141(46)=acc141(47)+acc141(33)+acc141(46)
      acc141(46)=Qspvak5e1*acc141(46)
      acc141(47)=acc141(36)*Qspvae1l3
      acc141(48)=acc141(29)*Qspval3e1
      acc141(49)=Qspvak2l4*acc141(26)
      acc141(50)=Qspvak2k5*acc141(11)
      acc141(51)=Qspval4k2*acc141(37)
      acc141(52)=Qspval4k5*acc141(16)
      acc141(53)=Qspvak5k2*acc141(8)
      acc141(54)=Qspvak5l4*acc141(17)
      acc141(55)=QspQ*acc141(6)
      brack=acc141(9)+acc141(41)+acc141(42)+acc141(43)+acc141(44)+acc141(45)+ac&
      &c141(46)+acc141(47)+acc141(48)+acc141(49)+acc141(50)+acc141(51)+acc141(5&
      &2)+acc141(53)+acc141(54)+acc141(55)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p2_part21part21_part25part25part21_d141h0l1_ninja")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki_nin
      use p2_part21part21_part25part25part21_globalsl1, only: epspow
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_abbrevd141h0
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d141
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k3+k5+k4
      Q(1:4)  =cmplx(real(+Q_ext(0:3)  -qshift(:),  ki_nin), aimag(+Q_ext(0:3))&
      &, ki)
      d141 = 0.0_ki
      d141 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d141, ki), aimag(d141), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p2_part21part21_part25part25part21_d141h0l1
