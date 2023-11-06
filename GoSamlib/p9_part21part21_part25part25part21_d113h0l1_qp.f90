module     p9_part21part21_part25part25part21_d113h0l1_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d113h0l1_qp.f90
   ! generator: buildfortran.py
   use p9_part21part21_part25part25part21_config, only: ki => ki_qp
   use p9_part21part21_part25part25part21_util_qp, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_ninja
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_color_qp
      use p9_part21part21_part25part25part21_abbrevd113h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc113(80)
      complex(ki) :: Qspvak1e5
      complex(ki) :: Qspvae5k1
      complex(ki) :: Qspvak2e5
      complex(ki) :: Qspvae5k2
      complex(ki) :: Qspvae1e5
      complex(ki) :: Qspvae5e1
      complex(ki) :: QspQ
      complex(ki) :: Qspk2
      complex(ki) :: Qspval4l3
      complex(ki) :: Qspval3l4
      complex(ki) :: Qspval3k2
      complex(ki) :: Qspvak2l3
      complex(ki) :: Qspval4k2
      complex(ki) :: Qspvak2l4
      complex(ki) :: Qspvae1l4
      complex(ki) :: Qspval4e1
      complex(ki) :: Qspvae5l3
      complex(ki) :: Qspval3e5
      complex(ki) :: Qspvae1k2
      complex(ki) :: Qspvak2e1
      complex(ki) :: Qspval4k1
      complex(ki) :: Qspvak2k1
      complex(ki) :: Qspvak1l4
      complex(ki) :: Qspvak1k2
      Qspvak1e5 = dotproduct(Q,spvak1e5)
      Qspvae5k1 = dotproduct(Q,spvae5k1)
      Qspvak2e5 = dotproduct(Q,spvak2e5)
      Qspvae5k2 = dotproduct(Q,spvae5k2)
      Qspvae1e5 = dotproduct(Q,spvae1e5)
      Qspvae5e1 = dotproduct(Q,spvae5e1)
      QspQ = dotproduct(Q,Q)
      Qspk2 = dotproduct(Q,k2)
      Qspval4l3 = dotproduct(Q,spval4l3)
      Qspval3l4 = dotproduct(Q,spval3l4)
      Qspval3k2 = dotproduct(Q,spval3k2)
      Qspvak2l3 = dotproduct(Q,spvak2l3)
      Qspval4k2 = dotproduct(Q,spval4k2)
      Qspvak2l4 = dotproduct(Q,spvak2l4)
      Qspvae1l4 = dotproduct(Q,spvae1l4)
      Qspval4e1 = dotproduct(Q,spval4e1)
      Qspvae5l3 = dotproduct(Q,spvae5l3)
      Qspval3e5 = dotproduct(Q,spval3e5)
      Qspvae1k2 = dotproduct(Q,spvae1k2)
      Qspvak2e1 = dotproduct(Q,spvak2e1)
      Qspval4k1 = dotproduct(Q,spval4k1)
      Qspvak2k1 = dotproduct(Q,spvak2k1)
      Qspvak1l4 = dotproduct(Q,spvak1l4)
      Qspvak1k2 = dotproduct(Q,spvak1k2)
      acc113(1)=abb113(7)
      acc113(2)=abb113(8)
      acc113(3)=abb113(9)
      acc113(4)=abb113(11)
      acc113(5)=abb113(12)
      acc113(6)=abb113(13)
      acc113(7)=abb113(14)
      acc113(8)=abb113(15)
      acc113(9)=abb113(16)
      acc113(10)=abb113(17)
      acc113(11)=abb113(18)
      acc113(12)=abb113(19)
      acc113(13)=abb113(20)
      acc113(14)=abb113(21)
      acc113(15)=abb113(22)
      acc113(16)=abb113(23)
      acc113(17)=abb113(24)
      acc113(18)=abb113(25)
      acc113(19)=abb113(26)
      acc113(20)=abb113(27)
      acc113(21)=abb113(28)
      acc113(22)=abb113(29)
      acc113(23)=abb113(31)
      acc113(24)=abb113(32)
      acc113(25)=abb113(33)
      acc113(26)=abb113(34)
      acc113(27)=abb113(37)
      acc113(28)=abb113(42)
      acc113(29)=abb113(43)
      acc113(30)=abb113(44)
      acc113(31)=abb113(45)
      acc113(32)=abb113(46)
      acc113(33)=abb113(47)
      acc113(34)=abb113(48)
      acc113(35)=abb113(49)
      acc113(36)=abb113(50)
      acc113(37)=abb113(51)
      acc113(38)=abb113(54)
      acc113(39)=abb113(56)
      acc113(40)=abb113(57)
      acc113(41)=abb113(58)
      acc113(42)=abb113(59)
      acc113(43)=abb113(60)
      acc113(44)=abb113(61)
      acc113(45)=abb113(65)
      acc113(46)=abb113(66)
      acc113(47)=abb113(67)
      acc113(48)=abb113(69)
      acc113(49)=abb113(71)
      acc113(50)=abb113(73)
      acc113(51)=abb113(75)
      acc113(52)=abb113(79)
      acc113(53)=abb113(82)
      acc113(54)=abb113(83)
      acc113(55)=abb113(87)
      acc113(56)=Qspvak1e5*Qspvae5k1
      acc113(57)=Qspvak2e5*Qspvae5k2
      acc113(56)=acc113(56)-acc113(57)
      acc113(57)=acc113(52)*acc113(56)
      acc113(58)=acc113(9)*Qspvae1e5
      acc113(59)=acc113(13)*Qspvae5e1
      acc113(60)=acc113(27)*Qspvak1e5
      acc113(61)=acc113(30)*Qspvae5k1
      acc113(62)=acc113(40)*Qspvae5k2
      acc113(63)=acc113(44)*Qspvak2e5
      acc113(57)=acc113(63)+acc113(62)+acc113(61)+acc113(60)+acc113(59)+acc113(&
      &58)+acc113(6)+acc113(57)
      acc113(57)=QspQ*acc113(57)
      acc113(58)=acc113(31)*Qspvak1e5
      acc113(59)=acc113(35)*Qspvae5k1
      acc113(60)=acc113(37)*Qspvae5e1
      acc113(61)=acc113(38)*Qspvae1e5
      acc113(58)=acc113(47)+acc113(61)+acc113(60)+acc113(59)+acc113(58)
      acc113(58)=Qspk2*acc113(58)
      acc113(59)=acc113(39)*Qspvae5k1
      acc113(60)=acc113(46)*Qspvae5e1
      acc113(61)=acc113(54)*Qspvae5k2
      acc113(59)=acc113(61)+acc113(60)+acc113(59)+acc113(32)
      acc113(59)=Qspval4l3*acc113(59)
      acc113(60)=acc113(3)*Qspvae1e5
      acc113(61)=acc113(29)*Qspvak1e5
      acc113(62)=acc113(42)*Qspvak2e5
      acc113(60)=acc113(62)+acc113(61)+acc113(60)+acc113(2)
      acc113(60)=Qspval3l4*acc113(60)
      acc113(61)=acc113(11)*Qspvae1e5
      acc113(62)=acc113(23)*Qspvak1e5
      acc113(63)=acc113(26)*Qspvak2e5
      acc113(61)=acc113(63)+acc113(62)+acc113(61)+acc113(10)
      acc113(61)=Qspval3k2*acc113(61)
      acc113(62)=acc113(1)*Qspvae5k1
      acc113(63)=acc113(8)*Qspvae5e1
      acc113(64)=acc113(15)*Qspvae5k2
      acc113(62)=acc113(64)+acc113(63)+acc113(4)+acc113(62)
      acc113(62)=Qspvak2l3*acc113(62)
      acc113(63)=acc113(21)*Qspvae5e1
      acc113(64)=acc113(51)*Qspvae5k1
      acc113(63)=acc113(64)+acc113(63)+acc113(18)
      acc113(63)=Qspval4k2*acc113(63)
      acc113(64)=acc113(7)*Qspvae1e5
      acc113(65)=acc113(12)*Qspvak1e5
      acc113(64)=acc113(65)+acc113(64)+acc113(5)
      acc113(64)=Qspvak2l4*acc113(64)
      acc113(56)=-acc113(48)*acc113(56)
      acc113(65)=acc113(33)*QspQ
      acc113(65)=acc113(65)+acc113(22)
      acc113(65)=acc113(65)*Qspvae5e1
      acc113(65)=acc113(16)+acc113(65)
      acc113(65)=Qspvae1e5*acc113(65)
      acc113(66)=acc113(17)*Qspvae5e1
      acc113(67)=acc113(41)*Qspvak1e5
      acc113(68)=acc113(45)*Qspvae5k2
      acc113(69)=acc113(50)*Qspvak2e5
      acc113(70)=acc113(53)*Qspvae5k1
      acc113(71)=Qspvae1l4*acc113(34)
      acc113(72)=Qspval4e1*acc113(43)
      acc113(73)=Qspvae5l3*acc113(36)
      acc113(74)=Qspval3e5*acc113(25)
      acc113(75)=Qspvae1k2*acc113(49)
      acc113(76)=Qspvak2e1*acc113(28)
      acc113(77)=Qspval4k1*acc113(55)
      acc113(78)=Qspvak2k1*acc113(14)
      acc113(79)=Qspvak1l4*acc113(24)
      acc113(80)=Qspvak1k2*acc113(19)
      brack=acc113(20)+acc113(56)+acc113(57)+acc113(58)+acc113(59)+acc113(60)+a&
      &cc113(61)+acc113(62)+acc113(63)+acc113(64)+acc113(65)+acc113(66)+acc113(&
      &67)+acc113(68)+acc113(69)+acc113(70)+acc113(71)+acc113(72)+acc113(73)+ac&
      &c113(74)+acc113(75)+acc113(76)+acc113(77)+acc113(78)+acc113(79)+acc113(8&
      &0)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p9_part21part21_part25part25part21_d113h0l1_qp_ninja")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki_nin
      use p9_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_abbrevd113h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d113
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k5+k4
      Q(1:4)  =cmplx(real(+Q_ext(0:3)  -qshift(:),  ki_nin), aimag(+Q_ext(0:3))&
      &, ki)
      d113 = 0.0_ki
      d113 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d113, ki), aimag(d113), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p9_part21part21_part25part25part21_d113h0l1_qp
