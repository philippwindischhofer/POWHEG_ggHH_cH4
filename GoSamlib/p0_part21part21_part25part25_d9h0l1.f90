module     p0_part21part21_part25part25_d9h0l1
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p0_part21part &
   ! &21_part25part25/helicity0d9h0l1.f90
   ! generator: buildfortran.py
   use p0_part21part21_part25part25_config, only: ki
   use p0_part21part21_part25part25_util, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_golem95
   public :: numerator_ninja
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p0_part21part21_part25part25_model
      use p0_part21part21_part25part25_kinematics
      use p0_part21part21_part25part25_color
      use p0_part21part21_part25part25_abbrevd9h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc9(47)
      complex(ki) :: Qspvae1k1
      complex(ki) :: Qspvak1e2
      complex(ki) :: Qspval3e2
      complex(ki) :: Qspvae1l3
      complex(ki) :: Qspvak2k1
      complex(ki) :: Qspvae1k2
      complex(ki) :: Qspvak2e2
      complex(ki) :: QspQ
      complex(ki) :: Qspvae1e2
      complex(ki) :: Qspvae2e1
      complex(ki) :: Qspvak1e1
      complex(ki) :: Qspvae2k1
      complex(ki) :: Qspvae2l3
      complex(ki) :: Qspval3e1
      complex(ki) :: Qspvak1k2
      complex(ki) :: Qspvak2e1
      complex(ki) :: Qspvae2k2
      Qspvae1k1 = dotproduct(Q,spvae1k1)
      Qspvak1e2 = dotproduct(Q,spvak1e2)
      Qspval3e2 = dotproduct(Q,spval3e2)
      Qspvae1l3 = dotproduct(Q,spvae1l3)
      Qspvak2k1 = dotproduct(Q,spvak2k1)
      Qspvae1k2 = dotproduct(Q,spvae1k2)
      Qspvak2e2 = dotproduct(Q,spvak2e2)
      QspQ = dotproduct(Q,Q)
      Qspvae1e2 = dotproduct(Q,spvae1e2)
      Qspvae2e1 = dotproduct(Q,spvae2e1)
      Qspvak1e1 = dotproduct(Q,spvak1e1)
      Qspvae2k1 = dotproduct(Q,spvae2k1)
      Qspvae2l3 = dotproduct(Q,spvae2l3)
      Qspval3e1 = dotproduct(Q,spval3e1)
      Qspvak1k2 = dotproduct(Q,spvak1k2)
      Qspvak2e1 = dotproduct(Q,spvak2e1)
      Qspvae2k2 = dotproduct(Q,spvae2k2)
      acc9(1)=abb9(4)
      acc9(2)=abb9(5)
      acc9(3)=abb9(7)
      acc9(4)=abb9(8)
      acc9(5)=abb9(9)
      acc9(6)=abb9(10)
      acc9(7)=abb9(11)
      acc9(8)=abb9(12)
      acc9(9)=abb9(13)
      acc9(10)=abb9(14)
      acc9(11)=abb9(15)
      acc9(12)=abb9(16)
      acc9(13)=abb9(17)
      acc9(14)=abb9(18)
      acc9(15)=abb9(19)
      acc9(16)=abb9(20)
      acc9(17)=abb9(21)
      acc9(18)=abb9(22)
      acc9(19)=abb9(23)
      acc9(20)=abb9(24)
      acc9(21)=abb9(25)
      acc9(22)=abb9(26)
      acc9(23)=abb9(27)
      acc9(24)=abb9(28)
      acc9(25)=abb9(29)
      acc9(26)=abb9(30)
      acc9(27)=abb9(31)
      acc9(28)=abb9(32)
      acc9(29)=abb9(33)
      acc9(30)=abb9(34)
      acc9(31)=abb9(35)
      acc9(32)=abb9(36)
      acc9(33)=Qspvae1k1*Qspvak1e2
      acc9(34)=-Qspval3e2*Qspvae1l3
      acc9(33)=acc9(33)+acc9(34)
      acc9(33)=acc9(1)*acc9(33)
      acc9(34)=acc9(21)*Qspvak2k1
      acc9(35)=Qspvae1k1*acc9(31)
      acc9(36)=-acc9(10)*Qspvae1k2
      acc9(36)=acc9(13)+acc9(36)
      acc9(36)=Qspvak2e2*acc9(36)
      acc9(37)=Qspvae1l3*acc9(22)
      acc9(38)=Qspval3e2*acc9(8)
      acc9(39)=QspQ*acc9(18)
      acc9(40)=QspQ*acc9(20)
      acc9(40)=acc9(5)+acc9(40)
      acc9(40)=Qspvae1e2*acc9(40)
      acc9(33)=acc9(40)+acc9(33)+acc9(39)+acc9(38)+acc9(37)+acc9(36)+acc9(35)+a&
      &cc9(34)+acc9(14)
      acc9(33)=Qspvae2e1*acc9(33)
      acc9(34)=Qspvak1e1*Qspvae2k1
      acc9(35)=-Qspvae2l3*Qspval3e1
      acc9(34)=acc9(34)+acc9(35)
      acc9(34)=acc9(1)*acc9(34)
      acc9(35)=acc9(23)*Qspvak1k2
      acc9(36)=Qspvak1e1*acc9(32)
      acc9(37)=-acc9(10)*Qspvak2e1
      acc9(37)=acc9(17)+acc9(37)
      acc9(37)=Qspvae2k2*acc9(37)
      acc9(38)=Qspval3e1*acc9(24)
      acc9(39)=Qspvae2l3*acc9(4)
      acc9(40)=QspQ*acc9(19)
      acc9(34)=acc9(34)+acc9(40)+acc9(39)+acc9(38)+acc9(37)+acc9(36)+acc9(35)+a&
      &cc9(7)
      acc9(34)=Qspvae1e2*acc9(34)
      acc9(35)=Qspvak1e2*acc9(30)
      acc9(36)=Qspvae2k1*acc9(29)
      acc9(37)=Qspvak2e1*acc9(16)
      acc9(38)=Qspvae1k2*acc9(12)
      acc9(39)=Qspvak1e1*acc9(6)
      acc9(40)=Qspvae1k1*acc9(9)
      acc9(41)=Qspvak2e2*acc9(15)
      acc9(42)=Qspvae2k2*acc9(11)
      acc9(43)=Qspval3e1*acc9(28)
      acc9(44)=Qspvae1l3*acc9(27)
      acc9(45)=Qspval3e2*acc9(26)
      acc9(46)=Qspvae2l3*acc9(25)
      acc9(47)=QspQ*acc9(3)
      brack=acc9(2)+acc9(33)+acc9(34)+acc9(35)+acc9(36)+acc9(37)+acc9(38)+acc9(&
      &39)+acc9(40)+acc9(41)+acc9(42)+acc9(43)+acc9(44)+acc9(45)+acc9(46)+acc9(&
      &47)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p0_part21part21_part25part25_globalsl1, only: epspow
      use p0_part21part21_part25part25_kinematics
      use p0_part21part21_part25part25_abbrevd9h0
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d9
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = k2
      Q(:)  =cmplx(real(+Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d9 = 0.0_ki
      d9 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d9, ki), aimag(d9), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p0_part21part21_part25part25_d9h0l1_ninja")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki_nin
      use p0_part21part21_part25part25_globalsl1, only: epspow
      use p0_part21part21_part25part25_kinematics
      use p0_part21part21_part25part25_abbrevd9h0
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d9
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k2
      Q(1:4)  =cmplx(real(+Q_ext(0:3)  -qshift(:),  ki_nin), aimag(+Q_ext(0:3))&
      &, ki)
      d9 = 0.0_ki
      d9 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d9, ki), aimag(d9), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p0_part21part21_part25part25_d9h0l1
