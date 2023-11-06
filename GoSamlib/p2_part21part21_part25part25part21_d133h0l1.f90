module     p2_part21part21_part25part25part21_d133h0l1
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d133h0l1.f90
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
      use p2_part21part21_part25part25part21_abbrevd133h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc133(62)
      complex(ki) :: Qspvak1k2
      complex(ki) :: Qspvak2k1
      complex(ki) :: Qspvak1k5
      complex(ki) :: Qspvak5k1
      complex(ki) :: Qspvak1e1
      complex(ki) :: Qspvae1k1
      complex(ki) :: Qspval4e1
      complex(ki) :: Qspvae1l4
      complex(ki) :: Qspvak2e1
      complex(ki) :: Qspvae1k2
      complex(ki) :: Qspval3e1
      complex(ki) :: Qspvae1l3
      complex(ki) :: QspQ
      complex(ki) :: Qspvak5e1
      complex(ki) :: Qspvae1k5
      Qspvak1k2 = dotproduct(Q,spvak1k2)
      Qspvak2k1 = dotproduct(Q,spvak2k1)
      Qspvak1k5 = dotproduct(Q,spvak1k5)
      Qspvak5k1 = dotproduct(Q,spvak5k1)
      Qspvak1e1 = dotproduct(Q,spvak1e1)
      Qspvae1k1 = dotproduct(Q,spvae1k1)
      Qspval4e1 = dotproduct(Q,spval4e1)
      Qspvae1l4 = dotproduct(Q,spvae1l4)
      Qspvak2e1 = dotproduct(Q,spvak2e1)
      Qspvae1k2 = dotproduct(Q,spvae1k2)
      Qspval3e1 = dotproduct(Q,spval3e1)
      Qspvae1l3 = dotproduct(Q,spvae1l3)
      QspQ = dotproduct(Q,Q)
      Qspvak5e1 = dotproduct(Q,spvak5e1)
      Qspvae1k5 = dotproduct(Q,spvae1k5)
      acc133(1)=abb133(7)
      acc133(2)=abb133(8)
      acc133(3)=abb133(9)
      acc133(4)=abb133(10)
      acc133(5)=abb133(11)
      acc133(6)=abb133(12)
      acc133(7)=abb133(14)
      acc133(8)=abb133(15)
      acc133(9)=abb133(16)
      acc133(10)=abb133(17)
      acc133(11)=abb133(18)
      acc133(12)=abb133(19)
      acc133(13)=abb133(20)
      acc133(14)=abb133(21)
      acc133(15)=abb133(22)
      acc133(16)=abb133(23)
      acc133(17)=abb133(24)
      acc133(18)=abb133(25)
      acc133(19)=abb133(26)
      acc133(20)=abb133(27)
      acc133(21)=abb133(28)
      acc133(22)=abb133(29)
      acc133(23)=abb133(30)
      acc133(24)=abb133(31)
      acc133(25)=abb133(32)
      acc133(26)=abb133(33)
      acc133(27)=abb133(34)
      acc133(28)=abb133(35)
      acc133(29)=abb133(36)
      acc133(30)=abb133(37)
      acc133(31)=abb133(38)
      acc133(32)=abb133(39)
      acc133(33)=abb133(40)
      acc133(34)=abb133(41)
      acc133(35)=abb133(42)
      acc133(36)=abb133(43)
      acc133(37)=abb133(44)
      acc133(38)=abb133(45)
      acc133(39)=abb133(46)
      acc133(40)=abb133(47)
      acc133(41)=abb133(48)
      acc133(42)=abb133(49)
      acc133(43)=abb133(50)
      acc133(44)=abb133(51)
      acc133(45)=abb133(52)
      acc133(46)=abb133(53)
      acc133(47)=abb133(55)
      acc133(48)=Qspvak1k2*acc133(20)
      acc133(49)=Qspvak2k1*acc133(44)
      acc133(50)=Qspvak1k5*acc133(45)
      acc133(51)=Qspvak5k1*acc133(2)
      acc133(52)=Qspvak1e1*acc133(35)
      acc133(53)=Qspvae1k1*acc133(19)
      acc133(54)=Qspval4e1*acc133(25)
      acc133(55)=Qspvae1l4*acc133(21)
      acc133(56)=Qspvak2e1*acc133(10)
      acc133(57)=Qspvae1k2*acc133(8)
      acc133(58)=Qspval3e1*acc133(38)
      acc133(59)=Qspvae1l3*acc133(33)
      acc133(60)=QspQ*acc133(6)
      acc133(48)=acc133(60)+acc133(59)+acc133(58)+acc133(57)+acc133(56)+acc133(&
      &55)+acc133(54)+acc133(53)+acc133(52)+acc133(51)+acc133(50)+acc133(49)+ac&
      &c133(5)+acc133(48)
      acc133(48)=QspQ*acc133(48)
      acc133(49)=Qspvak2k1*acc133(1)
      acc133(50)=Qspvak5k1*acc133(43)
      acc133(51)=Qspvae1k1*acc133(27)
      acc133(52)=Qspvae1l4*acc133(28)
      acc133(53)=Qspvae1k2*acc133(29)
      acc133(49)=acc133(53)+acc133(52)+acc133(51)+acc133(50)+acc133(36)+acc133(&
      &49)
      acc133(49)=Qspval3e1*acc133(49)
      acc133(50)=Qspvak1k2*acc133(23)
      acc133(51)=Qspvak1k5*acc133(47)
      acc133(52)=Qspvak1e1*acc133(39)
      acc133(53)=Qspval4e1*acc133(26)
      acc133(54)=Qspvak2e1*acc133(17)
      acc133(50)=acc133(54)+acc133(53)+acc133(52)+acc133(51)+acc133(31)+acc133(&
      &50)
      acc133(50)=Qspvae1l3*acc133(50)
      acc133(51)=Qspvak1k5*acc133(32)
      acc133(52)=Qspvak1e1*acc133(40)
      acc133(53)=Qspval4e1*acc133(34)
      acc133(54)=Qspvak2e1*acc133(4)
      acc133(51)=acc133(54)+acc133(53)+acc133(52)+acc133(11)+acc133(51)
      acc133(51)=Qspvae1k2*acc133(51)
      acc133(52)=Qspvak5k1*acc133(13)
      acc133(53)=Qspvae1k1*acc133(30)
      acc133(54)=Qspvae1l4*acc133(16)
      acc133(52)=acc133(54)+acc133(53)+acc133(7)+acc133(52)
      acc133(52)=Qspvak2e1*acc133(52)
      acc133(53)=Qspvak5e1*acc133(12)
      acc133(54)=Qspvak5e1*acc133(3)
      acc133(54)=acc133(9)+acc133(54)
      acc133(54)=Qspvae1k5*acc133(54)
      acc133(55)=Qspvak1k2*acc133(15)
      acc133(56)=Qspvak2k1*acc133(41)
      acc133(57)=Qspvak1k5*acc133(46)
      acc133(58)=Qspvak5k1*acc133(42)
      acc133(59)=Qspvak1e1*acc133(37)
      acc133(60)=Qspvae1k1*acc133(22)
      acc133(61)=Qspval4e1*acc133(24)
      acc133(62)=Qspvae1l4*acc133(18)
      brack=acc133(14)+acc133(48)+acc133(49)+acc133(50)+acc133(51)+acc133(52)+a&
      &cc133(53)+acc133(54)+acc133(55)+acc133(56)+acc133(57)+acc133(58)+acc133(&
      &59)+acc133(60)+acc133(61)+acc133(62)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p2_part21part21_part25part25part21_d133h0l1_ninja")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki_nin
      use p2_part21part21_part25part25part21_globalsl1, only: epspow
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_abbrevd133h0
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d133
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = -k3+k2-k5
      Q(1:4)  =cmplx(real(+Q_ext(0:3)  -qshift(:),  ki_nin), aimag(+Q_ext(0:3))&
      &, ki)
      d133 = 0.0_ki
      d133 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d133, ki), aimag(d133), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p2_part21part21_part25part25part21_d133h0l1
