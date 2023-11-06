module     p0_part21part21_part25part25_d7h0l1
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p0_part21part &
   ! &21_part25part25/helicity0d7h0l1.f90
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
      use p0_part21part21_part25part25_abbrevd7h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc7(16)
      complex(ki) :: Qspvae2e1
      complex(ki) :: Qspvae1e2
      complex(ki) :: Qspvae1k2
      complex(ki) :: Qspvak2e1
      complex(ki) :: Qspvae2k1
      complex(ki) :: Qspvak1e2
      complex(ki) :: QspQ
      Qspvae2e1 = dotproduct(Q,spvae2e1)
      Qspvae1e2 = dotproduct(Q,spvae1e2)
      Qspvae1k2 = dotproduct(Q,spvae1k2)
      Qspvak2e1 = dotproduct(Q,spvak2e1)
      Qspvae2k1 = dotproduct(Q,spvae2k1)
      Qspvak1e2 = dotproduct(Q,spvak1e2)
      QspQ = dotproduct(Q,Q)
      acc7(1)=abb7(3)
      acc7(2)=abb7(4)
      acc7(3)=abb7(5)
      acc7(4)=abb7(6)
      acc7(5)=abb7(8)
      acc7(6)=abb7(9)
      acc7(7)=abb7(10)
      acc7(8)=abb7(11)
      acc7(9)=abb7(12)
      acc7(10)=acc7(9)*Qspvae2e1
      acc7(10)=acc7(10)+acc7(8)
      acc7(10)=Qspvae1e2*acc7(10)
      acc7(11)=acc7(6)*Qspvae2e1
      acc7(12)=Qspvae1k2*acc7(7)
      acc7(13)=Qspvak2e1*acc7(5)
      acc7(14)=Qspvae2k1*acc7(4)
      acc7(15)=Qspvak1e2*acc7(3)
      acc7(16)=QspQ*acc7(2)
      brack=acc7(1)+acc7(10)+acc7(11)+acc7(12)+acc7(13)+acc7(14)+acc7(15)+acc7(&
      &16)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p0_part21part21_part25part25_globalsl1, only: epspow
      use p0_part21part21_part25part25_kinematics
      use p0_part21part21_part25part25_abbrevd7h0
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d7
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = k2
      Q(:)  =cmplx(real(+Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d7 = 0.0_ki
      d7 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d7, ki), aimag(d7), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p0_part21part21_part25part25_d7h0l1_ninja")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki_nin
      use p0_part21part21_part25part25_globalsl1, only: epspow
      use p0_part21part21_part25part25_kinematics
      use p0_part21part21_part25part25_abbrevd7h0
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d7
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k2
      Q(1:4)  =cmplx(real(+Q_ext(0:3)  -qshift(:),  ki_nin), aimag(+Q_ext(0:3))&
      &, ki)
      d7 = 0.0_ki
      d7 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d7, ki), aimag(d7), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p0_part21part21_part25part25_d7h0l1
