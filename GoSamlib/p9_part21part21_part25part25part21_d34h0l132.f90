module     p9_part21part21_part25part25part21_d34h0l132
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d34h0l132.f90
   ! generator: buildfortran_tn3.py
   use p9_part21part21_part25part25part21_config, only: ki
   use p9_part21part21_part25part25part21_util, only: cond_t, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, parameter :: ninjaidxt2x0mu0 = 0
   integer, parameter :: ninjaidxt1x0mu0 = 1
   integer, parameter :: ninjaidxt1x1mu0 = 2
   integer, parameter :: ninjaidxt0x0mu0 = 3
   integer, parameter :: ninjaidxt0x0mu2 = 4
   integer, parameter :: ninjaidxt0x1mu0 = 5
   integer, parameter :: ninjaidxt0x2mu0 = 6
   public :: numerator_t2
contains
!---#[ subroutine brack_21:
   pure subroutine brack_21(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd34h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(10) :: acd34
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd34(1)=dotproduct(ninjaE3,spvae1e5)
      acd34(2)=dotproduct(ninjaE3,spvae5e1)
      acd34(3)=abb34(8)
      acd34(4)=dotproduct(ninjaE3,spvak1e5)
      acd34(5)=dotproduct(ninjaE3,spvae5k1)
      acd34(6)=abb34(18)
      acd34(7)=dotproduct(ninjaE3,spvak2e5)
      acd34(8)=dotproduct(ninjaE3,spvae5k2)
      acd34(9)=acd34(5)*acd34(4)
      acd34(10)=-acd34(8)*acd34(7)
      acd34(9)=acd34(10)+acd34(9)
      acd34(9)=acd34(6)*acd34(9)
      acd34(10)=acd34(3)*acd34(2)*acd34(1)
      acd34(9)=acd34(10)+acd34(9)
      brack(ninjaidxt2x0mu0)=0.0_ki
      brack(ninjaidxt1x0mu0)=acd34(9)
      brack(ninjaidxt1x1mu0)=0.0_ki
   end subroutine brack_21
!---#] subroutine brack_21:
!---#[ subroutine brack_22:
   pure subroutine brack_22(ninjaA0, ninjaA1, ninjaE3, ninjaE4, ninjaP, brack)
      use p9_part21part21_part25part25part21_model
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_color
      use p9_part21part21_part25part25part21_abbrevd34h0
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA0, ninjaA1, ninjaE3, ninja&
      &E4
      complex(ki), dimension(0:2), intent(in) :: ninjaP
      complex(ki), dimension(56) :: acd34
      complex(ki), dimension (0:*), intent(inout) :: brack
      complex(ki) :: ninjaP0, ninjaP1, ninjaP2
      ninjaP0 = ninjaP(0)
      ninjaP1 = ninjaP(1)
      ninjaP2 = ninjaP(2)
      acd34(1)=dotproduct(ninjaA1,spvae5e1)
      acd34(2)=dotproduct(ninjaE3,spvae1e5)
      acd34(3)=abb34(8)
      acd34(4)=dotproduct(ninjaA1,spvae1e5)
      acd34(5)=dotproduct(ninjaE3,spvae5e1)
      acd34(6)=dotproduct(ninjaA1,spvak1e5)
      acd34(7)=dotproduct(ninjaE3,spvae5k1)
      acd34(8)=abb34(18)
      acd34(9)=dotproduct(ninjaA1,spvae5k1)
      acd34(10)=dotproduct(ninjaE3,spvak1e5)
      acd34(11)=dotproduct(ninjaA1,spvae5k2)
      acd34(12)=dotproduct(ninjaE3,spvak2e5)
      acd34(13)=dotproduct(ninjaA1,spvak2e5)
      acd34(14)=dotproduct(ninjaE3,spvae5k2)
      acd34(15)=dotproduct(ninjaA0,ninjaE3)
      acd34(16)=abb34(5)
      acd34(17)=dotproduct(ninjaA0,spvae5e1)
      acd34(18)=dotproduct(ninjaA0,spvae1e5)
      acd34(19)=dotproduct(ninjaA0,spvak1e5)
      acd34(20)=dotproduct(ninjaA0,spvae5k1)
      acd34(21)=dotproduct(ninjaA0,spvae5k2)
      acd34(22)=dotproduct(ninjaA0,spvak2e5)
      acd34(23)=abb34(6)
      acd34(24)=abb34(7)
      acd34(25)=abb34(9)
      acd34(26)=dotproduct(ninjaE3,spvae1k5)
      acd34(27)=abb34(10)
      acd34(28)=dotproduct(ninjaE3,spvak5e1)
      acd34(29)=abb34(11)
      acd34(30)=dotproduct(ninjaE3,spvak1k5)
      acd34(31)=abb34(12)
      acd34(32)=abb34(13)
      acd34(33)=abb34(15)
      acd34(34)=dotproduct(ninjaE3,spvak5k1)
      acd34(35)=abb34(16)
      acd34(36)=abb34(17)
      acd34(37)=dotproduct(ninjaE3,spvak5k2)
      acd34(38)=abb34(19)
      acd34(39)=dotproduct(ninjaE3,spvak2k5)
      acd34(40)=abb34(20)
      acd34(41)=-acd34(14)*acd34(13)
      acd34(42)=-acd34(12)*acd34(11)
      acd34(43)=acd34(10)*acd34(9)
      acd34(44)=acd34(7)*acd34(6)
      acd34(41)=acd34(44)+acd34(43)+acd34(41)+acd34(42)
      acd34(41)=acd34(8)*acd34(41)
      acd34(42)=acd34(5)*acd34(4)
      acd34(43)=acd34(2)*acd34(1)
      acd34(42)=acd34(42)+acd34(43)
      acd34(42)=acd34(3)*acd34(42)
      acd34(41)=acd34(42)+acd34(41)
      acd34(42)=-acd34(14)*acd34(22)
      acd34(43)=-acd34(12)*acd34(21)
      acd34(44)=acd34(10)*acd34(20)
      acd34(45)=acd34(7)*acd34(19)
      acd34(42)=acd34(45)+acd34(44)+acd34(42)+acd34(43)
      acd34(42)=acd34(8)*acd34(42)
      acd34(43)=acd34(5)*acd34(18)
      acd34(44)=acd34(2)*acd34(17)
      acd34(43)=acd34(43)+acd34(44)
      acd34(43)=acd34(3)*acd34(43)
      acd34(44)=acd34(39)*acd34(40)
      acd34(45)=acd34(37)*acd34(38)
      acd34(46)=acd34(34)*acd34(35)
      acd34(47)=acd34(30)*acd34(31)
      acd34(48)=acd34(28)*acd34(29)
      acd34(49)=acd34(26)*acd34(27)
      acd34(50)=acd34(15)*acd34(16)
      acd34(51)=acd34(14)*acd34(33)
      acd34(52)=acd34(12)*acd34(36)
      acd34(53)=acd34(10)*acd34(25)
      acd34(54)=acd34(7)*acd34(32)
      acd34(55)=acd34(5)*acd34(23)
      acd34(56)=acd34(2)*acd34(24)
      acd34(42)=acd34(42)+acd34(43)+acd34(56)+acd34(55)+acd34(54)+acd34(53)+acd&
      &34(52)+acd34(51)+2.0_ki*acd34(50)+acd34(49)+acd34(48)+acd34(47)+acd34(46&
      &)+acd34(44)+acd34(45)
      brack(ninjaidxt0x0mu0)=acd34(42)
      brack(ninjaidxt0x0mu2)=0.0_ki
      brack(ninjaidxt0x1mu0)=acd34(41)
      brack(ninjaidxt0x2mu0)=0.0_ki
   end subroutine brack_22
!---#] subroutine brack_22:
!---#[ subroutine numerator_t2:
   subroutine numerator_t2(ncut, a0, a1, b, c, param, deg, coeffs) &
   & bind(c, name="p9_part21part21_part25part25part21_d34h0_ninja_t2")
      use iso_c_binding, only: c_int
      use ninjago_module, only: ki => ki_nin
      use p9_part21part21_part25part25part21_globalsl1, only: epspow
      use p9_part21part21_part25part25part21_kinematics
      use p9_part21part21_part25part25part21_abbrevd34h0
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a0, a1, b, c
      complex(ki), dimension(0:2), intent(in) :: param
      complex(ki), dimension(4) :: vecA0, vecA1, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
	     vecA0(1:4) = - a0(0:3)
	     vecA1(1:4) = - a1(0:3)
	     vecB(1:4) = - b(0:3)
	     vecC(1:4) = - c(0:3)
      if (deg.lt.0) return
      t1 = 0
      call cond_t(epspow.eq.t1,brack_21,vecA0,vecA1,vecB,vecC,param,coeffs)
      if (deg.le.(1+(0))) return
      call cond_t(epspow.eq.t1,brack_22,vecA0,vecA1,vecB,vecC,param,coeffs)
   end subroutine numerator_t2
!---#] subroutine numerator_t2:
end module     p9_part21part21_part25part25part21_d34h0l132
