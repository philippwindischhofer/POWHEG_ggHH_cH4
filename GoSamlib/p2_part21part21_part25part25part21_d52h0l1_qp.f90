module     p2_part21part21_part25part25part21_d52h0l1_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d52h0l1_qp.f90
   ! generator: buildfortran.py
   use p2_part21part21_part25part25part21_config, only: ki => ki_qp
   use p2_part21part21_part25part25part21_util_qp, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_ninja
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_color_qp
      use p2_part21part21_part25part25part21_abbrevd52h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc52(48)
      complex(ki) :: Qspvae2k1
      complex(ki) :: Qspvae5k2
      complex(ki) :: Qspvae2k5
      complex(ki) :: Qspvae2e5
      complex(ki) :: QspQ
      complex(ki) :: Qspvae5e1
      complex(ki) :: Qspvae1e2
      complex(ki) :: Qspvak1e2
      complex(ki) :: Qspvak2e5
      complex(ki) :: Qspvak5e2
      complex(ki) :: Qspvae5e2
      complex(ki) :: Qspvae1e5
      complex(ki) :: Qspvae2e1
      Qspvae2k1 = dotproduct(Q,spvae2k1)
      Qspvae5k2 = dotproduct(Q,spvae5k2)
      Qspvae2k5 = dotproduct(Q,spvae2k5)
      Qspvae2e5 = dotproduct(Q,spvae2e5)
      QspQ = dotproduct(Q,Q)
      Qspvae5e1 = dotproduct(Q,spvae5e1)
      Qspvae1e2 = dotproduct(Q,spvae1e2)
      Qspvak1e2 = dotproduct(Q,spvak1e2)
      Qspvak2e5 = dotproduct(Q,spvak2e5)
      Qspvak5e2 = dotproduct(Q,spvak5e2)
      Qspvae5e2 = dotproduct(Q,spvae5e2)
      Qspvae1e5 = dotproduct(Q,spvae1e5)
      Qspvae2e1 = dotproduct(Q,spvae2e1)
      acc52(1)=abb52(2)
      acc52(2)=abb52(3)
      acc52(3)=abb52(4)
      acc52(4)=abb52(5)
      acc52(5)=abb52(6)
      acc52(6)=abb52(7)
      acc52(7)=abb52(8)
      acc52(8)=abb52(9)
      acc52(9)=abb52(10)
      acc52(10)=abb52(11)
      acc52(11)=abb52(12)
      acc52(12)=abb52(13)
      acc52(13)=abb52(14)
      acc52(14)=abb52(15)
      acc52(15)=abb52(16)
      acc52(16)=abb52(17)
      acc52(17)=abb52(18)
      acc52(18)=abb52(19)
      acc52(19)=abb52(20)
      acc52(20)=abb52(21)
      acc52(21)=abb52(22)
      acc52(22)=abb52(23)
      acc52(23)=abb52(24)
      acc52(24)=abb52(26)
      acc52(25)=abb52(27)
      acc52(26)=abb52(28)
      acc52(27)=abb52(30)
      acc52(28)=abb52(31)
      acc52(29)=abb52(32)
      acc52(30)=abb52(33)
      acc52(31)=abb52(34)
      acc52(32)=abb52(35)
      acc52(33)=abb52(36)
      acc52(34)=abb52(37)
      acc52(35)=abb52(38)
      acc52(36)=Qspvae2k1*acc52(12)
      acc52(37)=Qspvae5k2*acc52(34)
      acc52(38)=Qspvae2k5*acc52(30)
      acc52(39)=Qspvae2e5*acc52(1)
      acc52(40)=QspQ*acc52(15)
      acc52(41)=Qspvae2e5*acc52(11)
      acc52(41)=acc52(17)+acc52(41)
      acc52(41)=Qspvae5e1*acc52(41)
      acc52(36)=acc52(41)+acc52(40)+acc52(39)+acc52(38)+acc52(37)+acc52(9)+acc5&
      &2(36)
      acc52(36)=Qspvae1e2*acc52(36)
      acc52(37)=Qspvak1e2*acc52(18)
      acc52(38)=Qspvak2e5*acc52(14)
      acc52(39)=Qspvak5e2*acc52(32)
      acc52(40)=Qspvae5e2*acc52(4)
      acc52(41)=QspQ*acc52(23)
      acc52(42)=Qspvae5e2*acc52(11)
      acc52(42)=acc52(28)+acc52(42)
      acc52(42)=Qspvae1e5*acc52(42)
      acc52(37)=acc52(42)+acc52(41)+acc52(40)+acc52(39)+acc52(38)+acc52(8)+acc5&
      &2(37)
      acc52(37)=Qspvae2e1*acc52(37)
      acc52(38)=Qspvae2k1*acc52(26)
      acc52(39)=Qspvae5k2*acc52(16)
      acc52(40)=Qspvae2k5*acc52(25)
      acc52(41)=QspQ*acc52(27)
      acc52(38)=acc52(41)+acc52(40)+acc52(39)+acc52(13)+acc52(38)
      acc52(38)=Qspvae1e5*acc52(38)
      acc52(39)=Qspvak1e2*acc52(22)
      acc52(40)=Qspvak2e5*acc52(6)
      acc52(41)=Qspvak5e2*acc52(20)
      acc52(42)=QspQ*acc52(24)
      acc52(39)=acc52(42)+acc52(41)+acc52(40)+acc52(21)+acc52(39)
      acc52(39)=Qspvae5e1*acc52(39)
      acc52(40)=Qspvak1e2*acc52(2)
      acc52(41)=Qspvae2k1*acc52(5)
      acc52(42)=Qspvak2e5*acc52(35)
      acc52(43)=Qspvae5k2*acc52(33)
      acc52(44)=Qspvak5e2*acc52(31)
      acc52(45)=Qspvae2k5*acc52(29)
      acc52(46)=Qspvae2e5*acc52(10)
      acc52(47)=Qspvae5e2*acc52(19)
      acc52(48)=QspQ*acc52(3)
      brack=acc52(7)+acc52(36)+acc52(37)+acc52(38)+acc52(39)+acc52(40)+acc52(41&
      &)+acc52(42)+acc52(43)+acc52(44)+acc52(45)+acc52(46)+acc52(47)+acc52(48)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p2_part21part21_part25part25part21_d52h0l1_qp_ninja")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki_nin
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_abbrevd52h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d52
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k2
      Q(1:4)  =cmplx(real(+Q_ext(0:3)  -qshift(:),  ki_nin), aimag(+Q_ext(0:3))&
      &, ki)
      d52 = 0.0_ki
      d52 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d52, ki), aimag(d52), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p2_part21part21_part25part25part21_d52h0l1_qp
