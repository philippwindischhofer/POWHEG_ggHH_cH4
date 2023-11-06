module     p3_part1part21_part25part25part1_d11h0l1
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p3_part1part2 &
   ! &1_part25part25part1/helicity0d11h0l1.f90
   ! generator: buildfortran.py
   use p3_part1part21_part25part25part1_config, only: ki
   use p3_part1part21_part25part25part1_util, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_ninja
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p3_part1part21_part25part25part1_model
      use p3_part1part21_part25part25part1_kinematics
      use p3_part1part21_part25part25part1_color
      use p3_part1part21_part25part25part1_abbrevd11h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc11(12)
      complex(ki) :: Qspvak5e2
      complex(ki) :: Qspvae2k1
      complex(ki) :: Qspvak5k2
      complex(ki) :: Qspvak2k1
      complex(ki) :: QspQ
      Qspvak5e2 = dotproduct(Q,spvak5e2)
      Qspvae2k1 = dotproduct(Q,spvae2k1)
      Qspvak5k2 = dotproduct(Q,spvak5k2)
      Qspvak2k1 = dotproduct(Q,spvak2k1)
      QspQ = dotproduct(Q,Q)
      acc11(1)=abb11(4)
      acc11(2)=abb11(5)
      acc11(3)=abb11(6)
      acc11(4)=abb11(7)
      acc11(5)=abb11(8)
      acc11(6)=abb11(9)
      acc11(7)=abb11(12)
      acc11(8)=acc11(1)*Qspvak5e2
      acc11(8)=acc11(7)+acc11(8)
      acc11(8)=Qspvae2k1*acc11(8)
      acc11(9)=acc11(6)*Qspvak5e2
      acc11(10)=Qspvak5k2*acc11(5)
      acc11(11)=Qspvak2k1*acc11(4)
      acc11(12)=QspQ*acc11(2)
      brack=acc11(3)+acc11(8)+acc11(9)+acc11(10)+acc11(11)+acc11(12)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p3_part1part21_part25part25part1_d11h0l1_ninja")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki_nin
      use p3_part1part21_part25part25part1_globalsl1, only: epspow
      use p3_part1part21_part25part25part1_kinematics
      use p3_part1part21_part25part25part1_abbrevd11h0
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
      qshift = -k2
      Q(1:4)  =cmplx(real(-Q_ext(0:3)  -qshift(:),  ki_nin), aimag(-Q_ext(0:3))&
      &, ki)
      d11 = 0.0_ki
      d11 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d11, ki), aimag(d11), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p3_part1part21_part25part25part1_d11h0l1
