module     pb_part21part21_part25part25_d11h0l1
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/pb_part21part &
   ! &21_part25part25/helicity0d11h0l1.f90
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
      use pb_part21part21_part25part25_abbrevd11h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc11(47)
      complex(ki) :: Qspvae1k1
      complex(ki) :: Qspvak1e2
      complex(ki) :: Qspval4e2
      complex(ki) :: Qspvae1l4
      complex(ki) :: Qspvak2k1
      complex(ki) :: Qspvae1k2
      complex(ki) :: Qspvak2e2
      complex(ki) :: QspQ
      complex(ki) :: Qspvae1e2
      complex(ki) :: Qspvae2e1
      complex(ki) :: Qspvak1e1
      complex(ki) :: Qspvae2k1
      complex(ki) :: Qspvae2l4
      complex(ki) :: Qspval4e1
      complex(ki) :: Qspvak1k2
      complex(ki) :: Qspvak2e1
      complex(ki) :: Qspvae2k2
      Qspvae1k1 = dotproduct(Q,spvae1k1)
      Qspvak1e2 = dotproduct(Q,spvak1e2)
      Qspval4e2 = dotproduct(Q,spval4e2)
      Qspvae1l4 = dotproduct(Q,spvae1l4)
      Qspvak2k1 = dotproduct(Q,spvak2k1)
      Qspvae1k2 = dotproduct(Q,spvae1k2)
      Qspvak2e2 = dotproduct(Q,spvak2e2)
      QspQ = dotproduct(Q,Q)
      Qspvae1e2 = dotproduct(Q,spvae1e2)
      Qspvae2e1 = dotproduct(Q,spvae2e1)
      Qspvak1e1 = dotproduct(Q,spvak1e1)
      Qspvae2k1 = dotproduct(Q,spvae2k1)
      Qspvae2l4 = dotproduct(Q,spvae2l4)
      Qspval4e1 = dotproduct(Q,spval4e1)
      Qspvak1k2 = dotproduct(Q,spvak1k2)
      Qspvak2e1 = dotproduct(Q,spvak2e1)
      Qspvae2k2 = dotproduct(Q,spvae2k2)
      acc11(1)=abb11(4)
      acc11(2)=abb11(5)
      acc11(3)=abb11(7)
      acc11(4)=abb11(8)
      acc11(5)=abb11(9)
      acc11(6)=abb11(10)
      acc11(7)=abb11(11)
      acc11(8)=abb11(12)
      acc11(9)=abb11(13)
      acc11(10)=abb11(14)
      acc11(11)=abb11(15)
      acc11(12)=abb11(16)
      acc11(13)=abb11(17)
      acc11(14)=abb11(18)
      acc11(15)=abb11(19)
      acc11(16)=abb11(20)
      acc11(17)=abb11(21)
      acc11(18)=abb11(22)
      acc11(19)=abb11(23)
      acc11(20)=abb11(24)
      acc11(21)=abb11(25)
      acc11(22)=abb11(26)
      acc11(23)=abb11(27)
      acc11(24)=abb11(28)
      acc11(25)=abb11(29)
      acc11(26)=abb11(30)
      acc11(27)=abb11(31)
      acc11(28)=abb11(32)
      acc11(29)=abb11(33)
      acc11(30)=abb11(34)
      acc11(31)=abb11(35)
      acc11(32)=abb11(36)
      acc11(33)=Qspvae1k1*Qspvak1e2
      acc11(34)=-Qspval4e2*Qspvae1l4
      acc11(33)=acc11(33)+acc11(34)
      acc11(33)=acc11(1)*acc11(33)
      acc11(34)=acc11(21)*Qspvak2k1
      acc11(35)=Qspvae1k1*acc11(31)
      acc11(36)=-acc11(10)*Qspvae1k2
      acc11(36)=acc11(13)+acc11(36)
      acc11(36)=Qspvak2e2*acc11(36)
      acc11(37)=Qspvae1l4*acc11(22)
      acc11(38)=Qspval4e2*acc11(8)
      acc11(39)=QspQ*acc11(18)
      acc11(40)=QspQ*acc11(20)
      acc11(40)=acc11(5)+acc11(40)
      acc11(40)=Qspvae1e2*acc11(40)
      acc11(33)=acc11(40)+acc11(33)+acc11(39)+acc11(38)+acc11(37)+acc11(36)+acc&
      &11(35)+acc11(34)+acc11(14)
      acc11(33)=Qspvae2e1*acc11(33)
      acc11(34)=Qspvak1e1*Qspvae2k1
      acc11(35)=-Qspvae2l4*Qspval4e1
      acc11(34)=acc11(34)+acc11(35)
      acc11(34)=acc11(1)*acc11(34)
      acc11(35)=acc11(23)*Qspvak1k2
      acc11(36)=Qspvak1e1*acc11(32)
      acc11(37)=-acc11(10)*Qspvak2e1
      acc11(37)=acc11(17)+acc11(37)
      acc11(37)=Qspvae2k2*acc11(37)
      acc11(38)=Qspval4e1*acc11(24)
      acc11(39)=Qspvae2l4*acc11(4)
      acc11(40)=QspQ*acc11(19)
      acc11(34)=acc11(34)+acc11(40)+acc11(39)+acc11(38)+acc11(37)+acc11(36)+acc&
      &11(35)+acc11(7)
      acc11(34)=Qspvae1e2*acc11(34)
      acc11(35)=Qspvak1e2*acc11(30)
      acc11(36)=Qspvae2k1*acc11(29)
      acc11(37)=Qspvak2e1*acc11(16)
      acc11(38)=Qspvae1k2*acc11(12)
      acc11(39)=Qspvak1e1*acc11(6)
      acc11(40)=Qspvae1k1*acc11(9)
      acc11(41)=Qspvak2e2*acc11(15)
      acc11(42)=Qspvae2k2*acc11(11)
      acc11(43)=Qspval4e1*acc11(28)
      acc11(44)=Qspvae1l4*acc11(27)
      acc11(45)=Qspval4e2*acc11(26)
      acc11(46)=Qspvae2l4*acc11(25)
      acc11(47)=QspQ*acc11(3)
      brack=acc11(2)+acc11(33)+acc11(34)+acc11(35)+acc11(36)+acc11(37)+acc11(38&
      &)+acc11(39)+acc11(40)+acc11(41)+acc11(42)+acc11(43)+acc11(44)+acc11(45)+&
      &acc11(46)+acc11(47)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use pb_part21part21_part25part25_globalsl1, only: epspow
      use pb_part21part21_part25part25_kinematics
      use pb_part21part21_part25part25_abbrevd11h0
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d11
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = k2
      Q(:)  =cmplx(real(+Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d11 = 0.0_ki
      d11 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d11, ki), aimag(d11), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="pb_part21part21_part25part25_d11h0l1_ninja")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki_nin
      use pb_part21part21_part25part25_globalsl1, only: epspow
      use pb_part21part21_part25part25_kinematics
      use pb_part21part21_part25part25_abbrevd11h0
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
      qshift = k2
      Q(1:4)  =cmplx(real(+Q_ext(0:3)  -qshift(:),  ki_nin), aimag(+Q_ext(0:3))&
      &, ki)
      d11 = 0.0_ki
      d11 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d11, ki), aimag(d11), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module pb_part21part21_part25part25_d11h0l1
