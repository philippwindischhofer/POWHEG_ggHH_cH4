module     p9_part21part21_part25part25part21_d137h0l1_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d137h0l1_qp.f90
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
      use p9_part21part21_part25part25part21_abbrevd137h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc137(55)
      complex(ki) :: QspQ
      complex(ki) :: Qspk2
      complex(ki) :: Qspvak2e1
      complex(ki) :: Qspvak2k5
      complex(ki) :: Qspvak2l3
      complex(ki) :: Qspval3e1
      complex(ki) :: Qspvak5e1
      complex(ki) :: Qspvae1k2
      complex(ki) :: Qspvak5k2
      complex(ki) :: Qspval3k2
      complex(ki) :: Qspvae1l3
      complex(ki) :: Qspvae1k5
      complex(ki) :: Qspvak5l3
      complex(ki) :: Qspval3k5
      complex(ki) :: Qspval4e1
      complex(ki) :: Qspvae1l4
      QspQ = dotproduct(Q,Q)
      Qspk2 = dotproduct(Q,k2)
      Qspvak2e1 = dotproduct(Q,spvak2e1)
      Qspvak2k5 = dotproduct(Q,spvak2k5)
      Qspvak2l3 = dotproduct(Q,spvak2l3)
      Qspval3e1 = dotproduct(Q,spval3e1)
      Qspvak5e1 = dotproduct(Q,spvak5e1)
      Qspvae1k2 = dotproduct(Q,spvae1k2)
      Qspvak5k2 = dotproduct(Q,spvak5k2)
      Qspval3k2 = dotproduct(Q,spval3k2)
      Qspvae1l3 = dotproduct(Q,spvae1l3)
      Qspvae1k5 = dotproduct(Q,spvae1k5)
      Qspvak5l3 = dotproduct(Q,spvak5l3)
      Qspval3k5 = dotproduct(Q,spval3k5)
      Qspval4e1 = dotproduct(Q,spval4e1)
      Qspvae1l4 = dotproduct(Q,spvae1l4)
      acc137(1)=abb137(7)
      acc137(2)=abb137(8)
      acc137(3)=abb137(9)
      acc137(4)=abb137(10)
      acc137(5)=abb137(11)
      acc137(6)=abb137(12)
      acc137(7)=abb137(14)
      acc137(8)=abb137(15)
      acc137(9)=abb137(16)
      acc137(10)=abb137(17)
      acc137(11)=abb137(18)
      acc137(12)=abb137(19)
      acc137(13)=abb137(20)
      acc137(14)=abb137(21)
      acc137(15)=abb137(22)
      acc137(16)=abb137(23)
      acc137(17)=abb137(24)
      acc137(18)=abb137(25)
      acc137(19)=abb137(26)
      acc137(20)=abb137(27)
      acc137(21)=abb137(28)
      acc137(22)=abb137(29)
      acc137(23)=abb137(30)
      acc137(24)=abb137(31)
      acc137(25)=abb137(32)
      acc137(26)=abb137(33)
      acc137(27)=abb137(34)
      acc137(28)=abb137(35)
      acc137(29)=abb137(36)
      acc137(30)=abb137(37)
      acc137(31)=abb137(38)
      acc137(32)=abb137(39)
      acc137(33)=abb137(40)
      acc137(34)=abb137(41)
      acc137(35)=abb137(43)
      acc137(36)=abb137(45)
      acc137(37)=abb137(46)
      acc137(38)=abb137(47)
      acc137(39)=abb137(49)
      acc137(40)=abb137(53)
      acc137(41)=QspQ*acc137(5)
      acc137(42)=Qspk2*acc137(3)
      acc137(42)=acc137(41)+acc137(20)+acc137(42)
      acc137(42)=Qspvak2e1*acc137(42)
      acc137(43)=Qspk2*acc137(40)
      acc137(44)=Qspvak2k5*acc137(13)
      acc137(45)=QspQ*acc137(19)
      acc137(46)=acc137(2)*Qspvak2l3
      acc137(46)=acc137(30)+acc137(46)
      acc137(46)=Qspval3e1*acc137(46)
      acc137(47)=-Qspvak2k5*acc137(1)
      acc137(47)=acc137(10)+acc137(47)
      acc137(47)=Qspvak5e1*acc137(47)
      acc137(42)=acc137(42)+acc137(47)+acc137(46)+acc137(45)+acc137(44)+acc137(&
      &21)+acc137(43)
      acc137(42)=Qspvae1k2*acc137(42)
      acc137(43)=Qspk2*acc137(23)
      acc137(44)=Qspvak5k2*acc137(4)
      acc137(45)=QspQ*acc137(22)
      acc137(46)=acc137(2)*Qspval3k2
      acc137(46)=acc137(32)+acc137(46)
      acc137(46)=Qspvae1l3*acc137(46)
      acc137(47)=-Qspvak5k2*acc137(1)
      acc137(47)=acc137(7)+acc137(47)
      acc137(47)=Qspvae1k5*acc137(47)
      acc137(43)=acc137(47)+acc137(46)+acc137(45)+acc137(44)+acc137(18)+acc137(&
      &43)
      acc137(43)=Qspvak2e1*acc137(43)
      acc137(44)=QspQ*acc137(14)
      acc137(45)=acc137(2)*Qspvak5l3
      acc137(45)=acc137(25)+acc137(45)
      acc137(45)=Qspval3e1*acc137(45)
      acc137(41)=acc137(15)+acc137(41)
      acc137(41)=Qspvak5e1*acc137(41)
      acc137(41)=acc137(41)+acc137(45)+acc137(35)+acc137(44)
      acc137(41)=Qspvae1k5*acc137(41)
      acc137(44)=Qspvak2l3*acc137(31)
      acc137(45)=Qspvak5l3*acc137(36)
      acc137(44)=acc137(45)+acc137(39)+acc137(44)
      acc137(44)=Qspval3e1*acc137(44)
      acc137(45)=Qspval3k2*acc137(38)
      acc137(46)=Qspval3k5*acc137(27)
      acc137(45)=acc137(46)+acc137(33)+acc137(45)
      acc137(45)=Qspvae1l3*acc137(45)
      acc137(46)=QspQ*acc137(12)
      acc137(47)=acc137(2)*Qspval3k5
      acc137(47)=acc137(28)+acc137(47)
      acc137(47)=Qspvae1l3*acc137(47)
      acc137(46)=acc137(47)+acc137(34)+acc137(46)
      acc137(46)=Qspvak5e1*acc137(46)
      acc137(47)=acc137(29)*Qspval4e1
      acc137(48)=acc137(24)*Qspvae1l4
      acc137(49)=Qspvak2l3*acc137(26)
      acc137(50)=Qspvak2k5*acc137(11)
      acc137(51)=Qspval3k2*acc137(37)
      acc137(52)=Qspval3k5*acc137(16)
      acc137(53)=Qspvak5k2*acc137(8)
      acc137(54)=Qspvak5l3*acc137(17)
      acc137(55)=QspQ*acc137(6)
      brack=acc137(9)+acc137(41)+acc137(42)+acc137(43)+acc137(44)+acc137(45)+ac&
      &c137(46)+acc137(47)+acc137(48)+acc137(49)+acc137(50)+acc137(51)+acc137(5&
      &2)+acc137(53)+acc137(54)+acc137(55)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p9_part21part21_part25part25part21_d137h0l1_qp_ninja")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki_nin
      use p9_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_abbrevd137h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d137
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k3+k4
      Q(1:4)  =cmplx(real(+Q_ext(0:3)  -qshift(:),  ki_nin), aimag(+Q_ext(0:3))&
      &, ki)
      d137 = 0.0_ki
      d137 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d137, ki), aimag(d137), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p9_part21part21_part25part25part21_d137h0l1_qp
