module     p9_part21part21_part25part25part21_d34h0l1
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d34h0l1.f90
   ! generator: buildfortran.py
   use p9_part21part21_part25part25part21_config, only: ki
   use p9_part21part21_part25part25part21_util, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_ninja
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd34h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc34(29)
      complex(ki) :: Qspvae5e1
      complex(ki) :: Qspvae1e5
      complex(ki) :: Qspvae5k2
      complex(ki) :: Qspvak2e5
      complex(ki) :: Qspvae5k1
      complex(ki) :: Qspvak1e5
      complex(ki) :: Qspvae1k5
      complex(ki) :: Qspvak5e1
      complex(ki) :: Qspvak5k2
      complex(ki) :: Qspvak5k1
      complex(ki) :: Qspvak2k5
      complex(ki) :: Qspvak1k5
      complex(ki) :: QspQ
      Qspvae5e1 = dotproduct(Q,spvae5e1)
      Qspvae1e5 = dotproduct(Q,spvae1e5)
      Qspvae5k2 = dotproduct(Q,spvae5k2)
      Qspvak2e5 = dotproduct(Q,spvak2e5)
      Qspvae5k1 = dotproduct(Q,spvae5k1)
      Qspvak1e5 = dotproduct(Q,spvak1e5)
      Qspvae1k5 = dotproduct(Q,spvae1k5)
      Qspvak5e1 = dotproduct(Q,spvak5e1)
      Qspvak5k2 = dotproduct(Q,spvak5k2)
      Qspvak5k1 = dotproduct(Q,spvak5k1)
      Qspvak2k5 = dotproduct(Q,spvak2k5)
      Qspvak1k5 = dotproduct(Q,spvak1k5)
      QspQ = dotproduct(Q,Q)
      acc34(1)=abb34(4)
      acc34(2)=abb34(5)
      acc34(3)=abb34(6)
      acc34(4)=abb34(7)
      acc34(5)=abb34(8)
      acc34(6)=abb34(9)
      acc34(7)=abb34(10)
      acc34(8)=abb34(11)
      acc34(9)=abb34(12)
      acc34(10)=abb34(13)
      acc34(11)=abb34(15)
      acc34(12)=abb34(16)
      acc34(13)=abb34(17)
      acc34(14)=abb34(18)
      acc34(15)=abb34(19)
      acc34(16)=abb34(20)
      acc34(17)=acc34(5)*Qspvae5e1
      acc34(17)=acc34(17)+acc34(4)
      acc34(17)=Qspvae1e5*acc34(17)
      acc34(18)=-Qspvae5k2*acc34(14)
      acc34(18)=acc34(13)+acc34(18)
      acc34(18)=Qspvak2e5*acc34(18)
      acc34(19)=Qspvae5k1*acc34(14)
      acc34(19)=acc34(6)+acc34(19)
      acc34(19)=Qspvak1e5*acc34(19)
      acc34(20)=acc34(3)*Qspvae5e1
      acc34(21)=acc34(10)*Qspvae5k1
      acc34(22)=acc34(11)*Qspvae5k2
      acc34(23)=Qspvae1k5*acc34(7)
      acc34(24)=Qspvak5e1*acc34(8)
      acc34(25)=Qspvak5k2*acc34(15)
      acc34(26)=Qspvak5k1*acc34(12)
      acc34(27)=Qspvak2k5*acc34(16)
      acc34(28)=Qspvak1k5*acc34(9)
      acc34(29)=QspQ*acc34(2)
      brack=acc34(1)+acc34(17)+acc34(18)+acc34(19)+acc34(20)+acc34(21)+acc34(22&
      &)+acc34(23)+acc34(24)+acc34(25)+acc34(26)+acc34(27)+acc34(28)+acc34(29)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p9_part21part21_part25part25part21_d34h0l1_ninja")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki_nin
      use p9_part21part21_part25part25part21_globalsl1, only: epspow
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_abbrevd34h0
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d34
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      Q(1:4)  =cmplx(real(-Q_ext(0:3),  ki_nin), aimag(-Q_ext(0:3)), ki)
      d34 = 0.0_ki
      d34 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d34, ki), aimag(d34), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p9_part21part21_part25part25part21_d34h0l1
