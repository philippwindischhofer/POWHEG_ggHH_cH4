module     p9_part21part21_part25part25part21_d147h0l1
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d147h0l1.f90
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
      use p9_part21part21_part25part25part21_abbrevd147h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc147(44)
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
      acc147(1)=abb147(3)
      acc147(2)=abb147(4)
      acc147(3)=abb147(5)
      acc147(4)=abb147(6)
      acc147(5)=abb147(7)
      acc147(6)=abb147(8)
      acc147(7)=abb147(9)
      acc147(8)=abb147(10)
      acc147(9)=abb147(11)
      acc147(10)=abb147(12)
      acc147(11)=abb147(13)
      acc147(12)=abb147(14)
      acc147(13)=abb147(15)
      acc147(14)=abb147(16)
      acc147(15)=abb147(17)
      acc147(16)=abb147(18)
      acc147(17)=abb147(19)
      acc147(18)=abb147(20)
      acc147(19)=abb147(21)
      acc147(20)=abb147(23)
      acc147(21)=abb147(24)
      acc147(22)=abb147(25)
      acc147(23)=abb147(26)
      acc147(24)=abb147(27)
      acc147(25)=abb147(28)
      acc147(26)=abb147(29)
      acc147(27)=abb147(30)
      acc147(28)=abb147(31)
      acc147(29)=abb147(32)
      acc147(30)=abb147(33)
      acc147(31)=abb147(34)
      acc147(32)=abb147(35)
      acc147(33)=abb147(36)
      acc147(34)=Qspvak2e1*acc147(33)
      acc147(35)=Qspvak1e5*acc147(14)
      acc147(36)=Qspvak5e1*acc147(20)
      acc147(37)=Qspvae5e1*acc147(19)
      acc147(38)=QspQ*acc147(29)
      acc147(39)=Qspvae5e1*acc147(5)
      acc147(39)=acc147(1)+acc147(39)
      acc147(39)=Qspvae2e5*acc147(39)
      acc147(34)=acc147(39)+acc147(38)+acc147(37)+acc147(36)+acc147(35)+acc147(&
      &9)+acc147(34)
      acc147(34)=Qspvae1e2*acc147(34)
      acc147(35)=Qspvae1k2*acc147(32)
      acc147(36)=Qspvae5k1*acc147(23)
      acc147(37)=Qspvae1k5*acc147(31)
      acc147(38)=Qspvae1e5*acc147(15)
      acc147(39)=QspQ*acc147(26)
      acc147(40)=Qspvae1e5*acc147(5)
      acc147(40)=acc147(4)+acc147(40)
      acc147(40)=Qspvae5e2*acc147(40)
      acc147(35)=acc147(40)+acc147(39)+acc147(38)+acc147(37)+acc147(36)+acc147(&
      &13)+acc147(35)
      acc147(35)=Qspvae2e1*acc147(35)
      acc147(36)=Qspvae1k2*acc147(24)
      acc147(37)=Qspvae5k1*acc147(6)
      acc147(38)=-Qspvae1k5*acc147(12)
      acc147(39)=QspQ*acc147(25)
      acc147(36)=acc147(39)+acc147(38)+acc147(37)+acc147(10)+acc147(36)
      acc147(36)=Qspvae2e5*acc147(36)
      acc147(37)=Qspvak2e1*acc147(21)
      acc147(38)=Qspvak1e5*acc147(3)
      acc147(39)=-Qspvak5e1*acc147(16)
      acc147(40)=QspQ*acc147(22)
      acc147(37)=acc147(40)+acc147(39)+acc147(38)+acc147(7)+acc147(37)
      acc147(37)=Qspvae5e2*acc147(37)
      acc147(38)=Qspvak1e5*acc147(8)
      acc147(39)=Qspvae5k1*acc147(11)
      acc147(40)=Qspvak5e1*acc147(17)
      acc147(41)=Qspvae1k5*acc147(30)
      acc147(42)=Qspvae1e5*acc147(28)
      acc147(43)=Qspvae5e1*acc147(27)
      acc147(44)=QspQ*acc147(18)
      brack=acc147(2)+acc147(34)+acc147(35)+acc147(36)+acc147(37)+acc147(38)+ac&
      &c147(39)+acc147(40)+acc147(41)+acc147(42)+acc147(43)+acc147(44)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ subroutine numerator_ninja:
   subroutine numerator_ninja(ncut, Q_ext, mu2_ext, numerator) &
   & bind(c, name="p9_part21part21_part25part25part21_d147h0l1_ninja")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki_nin
      use p9_part21part21_part25part25part21_globalsl1, only: epspow
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_abbrevd147h0
      implicit none
      integer(c_int), intent(in) :: ncut
      complex(ki_nin), dimension(0:3), intent(in) :: Q_ext
      complex(ki_nin), intent(in) :: mu2_ext
      complex(ki_nin), intent(out) :: numerator
      complex(ki) :: d147
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = -k2
      Q(1:4)  =cmplx(real(-Q_ext(0:3)  -qshift(:),  ki_nin), aimag(-Q_ext(0:3))&
      &, ki)
      d147 = 0.0_ki
      d147 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d147, ki), aimag(d147), ki_nin)
   end subroutine numerator_ninja
   !------#] subroutine numerator_ninja:
!---#] numerator interfaces:
end module p9_part21part21_part25part25part21_d147h0l1
