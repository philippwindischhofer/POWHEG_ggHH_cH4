module     p2_part21part21_part25part25part21_d42h0l1
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d42h0l1.f90
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
      use p2_part21part21_part25part25part21_abbrevd42h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc42(20)
      complex(ki) :: Qspvae1k5
      complex(ki) :: Qspvak5e1
      complex(ki) :: Qspvae1k2
      complex(ki) :: Qspvak2e1
      complex(ki) :: Qspvak5k1
      complex(ki) :: Qspvak2k1
      complex(ki) :: Qspvak1k5
      complex(ki) :: Qspvak1k2
      complex(ki) :: QspQ
      Qspvae1k5 = dotproduct(Q,spvae1k5)
      Qspvak5e1 = dotproduct(Q,spvak5e1)
      Qspvae1k2 = dotproduct(Q,spvae1k2)
      Qspvak2e1 = dotproduct(Q,spvak2e1)
      Qspvak5k1 = dotproduct(Q,spvak5k1)
      Qspvak2k1 = dotproduct(Q,spvak2k1)
      Qspvak1k5 = dotproduct(Q,spvak1k5)
      Qspvak1k2 = dotproduct(Q,spvak1k2)
      QspQ = dotproduct(Q,Q)
      acc42(1)=abb42(4)
      acc42(2)=abb42(5)
      acc42(3)=abb42(6)
      acc42(4)=abb42(7)
      acc42(5)=abb42(8)
      acc42(6)=abb42(9)
      acc42(7)=abb42(11)
      acc42(8)=abb42(12)
      acc42(9)=abb42(13)
      acc42(10)=abb42(14)
      acc42(11)=abb42(15)
      acc42(12)=-Qspvae1k5*acc42(9)
      acc42(12)=acc42(5)+acc42(12)
      acc42(12)=Qspvak5e1*acc42(12)
      acc42(13)=-Qspvae1k2*acc42(9)
      acc42(13)=acc42(10)+acc42(13)
      acc42(13)=Qspvak2e1*acc42(13)
      acc42(14)=acc42(3)*Qspvae1k5
      acc42(15)=acc42(7)*Qspvae1k2
      acc42(16)=Qspvak5k1*acc42(11)
      acc42(17)=Qspvak2k1*acc42(8)
      acc42(18)=Qspvak1k5*acc42(4)
      acc42(19)=Qspvak1k2*acc42(2)
      acc42(20)=-QspQ*acc42(6)
      brack=acc42(1)+acc42(12)+acc42(13)+acc42(14)+acc42(15)+acc42(16)+acc42(17&
      &)+acc42(18)+acc42(19)+acc42(20)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p2_part21part21_part25part25part21_d42h0l1_ninja")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki_nin
      use p2_part21part21_part25part25part21_globalsl1, only: epspow
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_abbrevd42h0
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d42
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k2
      Q(1:4)  =cmplx(real(+Q_ext(0:3)  -qshift(:),  ki_nin), aimag(+Q_ext(0:3))&
      &, ki)
      d42 = 0.0_ki
      d42 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d42, ki), aimag(d42), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p2_part21part21_part25part25part21_d42h0l1
