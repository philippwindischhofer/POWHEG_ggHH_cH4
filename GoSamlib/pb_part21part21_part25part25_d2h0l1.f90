module     pb_part21part21_part25part25_d2h0l1
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/pb_part21part &
   ! &21_part25part25/helicity0d2h0l1.f90
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
      use pb_part21part21_part25part25_abbrevd2h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc2(16)
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
      acc2(1)=abb2(2)
      acc2(2)=abb2(3)
      acc2(3)=abb2(5)
      acc2(4)=abb2(6)
      acc2(5)=abb2(7)
      acc2(6)=abb2(8)
      acc2(7)=abb2(9)
      acc2(8)=abb2(10)
      acc2(9)=abb2(11)
      acc2(10)=acc2(7)*Qspvae2e1
      acc2(10)=acc2(10)+acc2(6)
      acc2(10)=Qspvae1e2*acc2(10)
      acc2(11)=acc2(5)*Qspvae2e1
      acc2(12)=Qspvae1k2*acc2(8)
      acc2(13)=Qspvak2e1*acc2(9)
      acc2(14)=Qspvae2k1*acc2(2)
      acc2(15)=Qspvak1e2*acc2(1)
      acc2(16)=QspQ*acc2(4)
      brack=acc2(3)+acc2(10)+acc2(11)+acc2(12)+acc2(13)+acc2(14)+acc2(15)+acc2(&
      &16)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use pb_part21part21_part25part25_globalsl1, only: epspow
      use pb_part21part21_part25part25_kinematics
      use pb_part21part21_part25part25_abbrevd2h0
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d2
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = k2
      Q(:)  =cmplx(real(+Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d2 = 0.0_ki
      d2 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d2, ki), aimag(d2), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="pb_part21part21_part25part25_d2h0l1_ninja")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki_nin
      use pb_part21part21_part25part25_globalsl1, only: epspow
      use pb_part21part21_part25part25_kinematics
      use pb_part21part21_part25part25_abbrevd2h0
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d2
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k2
      Q(1:4)  =cmplx(real(+Q_ext(0:3)  -qshift(:),  ki_nin), aimag(+Q_ext(0:3))&
      &, ki)
      d2 = 0.0_ki
      d2 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d2, ki), aimag(d2), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module pb_part21part21_part25part25_d2h0l1
