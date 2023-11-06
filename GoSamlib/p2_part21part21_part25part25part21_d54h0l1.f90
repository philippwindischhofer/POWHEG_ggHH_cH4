module     p2_part21part21_part25part25part21_d54h0l1
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/p2_part21part21_part25part &
   ! &25part21/helicity0d54h0l1.f90
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
      use p2_part21part21_part25part25part21_abbrevd54h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc54(44)
      complex(ki) :: Qspvak2e1
      complex(ki) :: Qspvak1e5
      complex(ki) :: Qspvak5e1
      complex(ki) :: Qspvae5e1
      complex(ki) :: QspQ
      complex(ki) :: Qspvae2e5
      complex(ki) :: Qspvae1e2
      complex(ki) :: Qspvae1k2
      complex(ki) :: Qspvae5k1
      complex(ki) :: Qspvae1k5
      complex(ki) :: Qspvae1e5
      complex(ki) :: Qspvae5e2
      complex(ki) :: Qspvae2e1
      Qspvak2e1 = dotproduct(Q,spvak2e1)
      Qspvak1e5 = dotproduct(Q,spvak1e5)
      Qspvak5e1 = dotproduct(Q,spvak5e1)
      Qspvae5e1 = dotproduct(Q,spvae5e1)
      QspQ = dotproduct(Q,Q)
      Qspvae2e5 = dotproduct(Q,spvae2e5)
      Qspvae1e2 = dotproduct(Q,spvae1e2)
      Qspvae1k2 = dotproduct(Q,spvae1k2)
      Qspvae5k1 = dotproduct(Q,spvae5k1)
      Qspvae1k5 = dotproduct(Q,spvae1k5)
      Qspvae1e5 = dotproduct(Q,spvae1e5)
      Qspvae5e2 = dotproduct(Q,spvae5e2)
      Qspvae2e1 = dotproduct(Q,spvae2e1)
      acc54(1)=abb54(2)
      acc54(2)=abb54(3)
      acc54(3)=abb54(4)
      acc54(4)=abb54(5)
      acc54(5)=abb54(6)
      acc54(6)=abb54(7)
      acc54(7)=abb54(8)
      acc54(8)=abb54(9)
      acc54(9)=abb54(10)
      acc54(10)=abb54(11)
      acc54(11)=abb54(12)
      acc54(12)=abb54(13)
      acc54(13)=abb54(14)
      acc54(14)=abb54(15)
      acc54(15)=abb54(16)
      acc54(16)=abb54(17)
      acc54(17)=abb54(18)
      acc54(18)=abb54(19)
      acc54(19)=abb54(20)
      acc54(20)=abb54(21)
      acc54(21)=abb54(22)
      acc54(22)=abb54(23)
      acc54(23)=abb54(24)
      acc54(24)=abb54(25)
      acc54(25)=abb54(26)
      acc54(26)=abb54(28)
      acc54(27)=abb54(29)
      acc54(28)=abb54(30)
      acc54(29)=abb54(31)
      acc54(30)=abb54(32)
      acc54(31)=abb54(33)
      acc54(32)=abb54(34)
      acc54(33)=abb54(35)
      acc54(34)=Qspvak2e1*acc54(33)
      acc54(35)=Qspvak1e5*acc54(4)
      acc54(36)=Qspvak5e1*acc54(31)
      acc54(37)=-Qspvae5e1*acc54(17)
      acc54(38)=QspQ*acc54(22)
      acc54(39)=-Qspvae5e1*acc54(6)
      acc54(39)=acc54(5)+acc54(39)
      acc54(39)=Qspvae2e5*acc54(39)
      acc54(34)=acc54(39)+acc54(38)+acc54(37)+acc54(36)+acc54(35)+acc54(1)+acc5&
      &4(34)
      acc54(34)=Qspvae1e2*acc54(34)
      acc54(35)=Qspvae1k2*acc54(32)
      acc54(36)=Qspvae5k1*acc54(25)
      acc54(37)=Qspvae1k5*acc54(30)
      acc54(38)=-Qspvae1e5*acc54(27)
      acc54(39)=QspQ*acc54(28)
      acc54(40)=-Qspvae1e5*acc54(6)
      acc54(40)=acc54(3)+acc54(40)
      acc54(40)=Qspvae5e2*acc54(40)
      acc54(35)=acc54(40)+acc54(39)+acc54(38)+acc54(37)+acc54(36)+acc54(10)+acc&
      &54(35)
      acc54(35)=Qspvae2e1*acc54(35)
      acc54(36)=Qspvae1k2*acc54(14)
      acc54(37)=-Qspvae5k1*acc54(13)
      acc54(38)=Qspvae1k5*acc54(20)
      acc54(39)=QspQ*acc54(8)
      acc54(36)=acc54(39)+acc54(38)+acc54(37)+acc54(16)+acc54(36)
      acc54(36)=Qspvae2e5*acc54(36)
      acc54(37)=Qspvak2e1*acc54(23)
      acc54(38)=-Qspvak1e5*acc54(15)
      acc54(39)=Qspvak5e1*acc54(18)
      acc54(40)=QspQ*acc54(24)
      acc54(37)=acc54(40)+acc54(39)+acc54(38)+acc54(19)+acc54(37)
      acc54(37)=Qspvae5e2*acc54(37)
      acc54(38)=Qspvak1e5*acc54(11)
      acc54(39)=Qspvae5k1*acc54(9)
      acc54(40)=Qspvak5e1*acc54(21)
      acc54(41)=Qspvae1k5*acc54(29)
      acc54(42)=Qspvae1e5*acc54(26)
      acc54(43)=Qspvae5e1*acc54(7)
      acc54(44)=QspQ*acc54(12)
      brack=acc54(2)+acc54(34)+acc54(35)+acc54(36)+acc54(37)+acc54(38)+acc54(39&
      &)+acc54(40)+acc54(41)+acc54(42)+acc54(43)+acc54(44)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p2_part21part21_part25part25part21_d54h0l1_ninja")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki_nin
      use p2_part21part21_part25part25part21_globalsl1, only: epspow
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_abbrevd54h0
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d54
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = -k2
      Q(1:4)  =cmplx(real(-Q_ext(0:3)  -qshift(:),  ki_nin), aimag(-Q_ext(0:3))&
      &, ki)
      d54 = 0.0_ki
      d54 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d54, ki), aimag(d54), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p2_part21part21_part25part25part21_d54h0l1
