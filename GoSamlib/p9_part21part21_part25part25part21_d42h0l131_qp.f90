module     p9_part21part21_part25part25part21_d42h0l131_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d42h0l131_qp.f90
   ! generator: buildfortran_tn3.py
   use p9_part21part21_part25part25part21_config, only: ki => ki_qp
   use p9_part21part21_part25part25part21_util_qp, only: cond_t, d => metric_te&
   &nsor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, parameter :: ninjaidxt3mu0 = 0
   integer, parameter :: ninjaidxt2mu0 = 1
   integer, parameter :: ninjaidxt1mu0 = 2
   integer, parameter :: ninjaidxt1mu2 = 3
   integer, parameter :: ninjaidxt0mu0 = 4
   integer, parameter :: ninjaidxt0mu2 = 5
   public :: numerator_t3
contains
!---#[ subroutine brack_31:
   pure subroutine brack_31(ninjaA, ninjaE3, ninjaE4, ninjaP, brack)
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_color_qp
      use p9_part21part21_part25part25part21_abbrevd42h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(7) :: acd42
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd42(1)=dotproduct(ninjaE3,spvak2e1)
      acd42(2)=dotproduct(ninjaE3,spvae1k2)
      acd42(3)=abb42(13)
      acd42(4)=dotproduct(ninjaE3,spvak5e1)
      acd42(5)=dotproduct(ninjaE3,spvae1k5)
      acd42(6)=-acd42(2)*acd42(1)
      acd42(7)=-acd42(5)*acd42(4)
      acd42(6)=acd42(6)+acd42(7)
      acd42(6)=acd42(3)*acd42(6)
      brack(ninjaidxt3mu0)=0.0_ki
      brack(ninjaidxt2mu0)=acd42(6)
   end subroutine brack_31
!---#] subroutine brack_31:
!---#[ subroutine brack_32:
   pure subroutine brack_32(ninjaA, ninjaE3, ninjaE4, ninjaP, brack)
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_color_qp
      use p9_part21part21_part25part25part21_abbrevd42h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: ninjaA, ninjaE3, ninjaE4
      complex(ki), intent(in) :: ninjaP
      complex(ki), dimension(45) :: acd42
      complex(ki), dimension (0:*), intent(inout) :: brack
      acd42(1)=dotproduct(ninjaE3,spvae1k5)
      acd42(2)=dotproduct(ninjaE4,spvak5e1)
      acd42(3)=abb42(13)
      acd42(4)=dotproduct(ninjaE3,spvak5e1)
      acd42(5)=dotproduct(ninjaE4,spvae1k5)
      acd42(6)=dotproduct(ninjaE3,spvae1k2)
      acd42(7)=dotproduct(ninjaE4,spvak2e1)
      acd42(8)=dotproduct(ninjaE3,spvak2e1)
      acd42(9)=dotproduct(ninjaE4,spvae1k2)
      acd42(10)=abb42(9)
      acd42(11)=dotproduct(ninjaA,ninjaE3)
      acd42(12)=dotproduct(ninjaA,spvae1k5)
      acd42(13)=dotproduct(ninjaA,spvak5e1)
      acd42(14)=dotproduct(ninjaA,spvae1k2)
      acd42(15)=dotproduct(ninjaA,spvak2e1)
      acd42(16)=dotproduct(ninjaE3,spvak1k2)
      acd42(17)=abb42(5)
      acd42(18)=abb42(6)
      acd42(19)=dotproduct(ninjaE3,spvak1k5)
      acd42(20)=abb42(7)
      acd42(21)=abb42(8)
      acd42(22)=abb42(11)
      acd42(23)=dotproduct(ninjaE3,spvak2k1)
      acd42(24)=abb42(12)
      acd42(25)=abb42(14)
      acd42(26)=dotproduct(ninjaE3,spvak5k1)
      acd42(27)=abb42(15)
      acd42(28)=dotproduct(ninjaA,ninjaA)
      acd42(29)=dotproduct(ninjaA,spvak1k2)
      acd42(30)=dotproduct(ninjaA,spvak1k5)
      acd42(31)=dotproduct(ninjaA,spvak2k1)
      acd42(32)=dotproduct(ninjaA,spvak5k1)
      acd42(33)=abb42(4)
      acd42(34)=acd42(8)*acd42(9)
      acd42(35)=acd42(6)*acd42(7)
      acd42(36)=acd42(4)*acd42(5)
      acd42(37)=acd42(1)*acd42(2)
      acd42(34)=acd42(34)+acd42(35)+acd42(36)+acd42(37)
      acd42(35)=-acd42(3)*acd42(34)
      acd42(35)=-acd42(10)+acd42(35)
      acd42(36)=-acd42(8)*acd42(14)
      acd42(37)=-acd42(6)*acd42(15)
      acd42(38)=-acd42(4)*acd42(12)
      acd42(39)=-acd42(1)*acd42(13)
      acd42(36)=acd42(39)+acd42(38)+acd42(36)+acd42(37)
      acd42(36)=acd42(3)*acd42(36)
      acd42(37)=acd42(27)*acd42(26)
      acd42(38)=acd42(24)*acd42(23)
      acd42(39)=acd42(20)*acd42(19)
      acd42(40)=acd42(17)*acd42(16)
      acd42(41)=acd42(10)*acd42(11)
      acd42(42)=acd42(8)*acd42(25)
      acd42(43)=acd42(6)*acd42(22)
      acd42(44)=acd42(4)*acd42(21)
      acd42(45)=acd42(1)*acd42(18)
      acd42(36)=acd42(36)+acd42(45)+acd42(44)+acd42(43)+acd42(42)-2.0_ki*acd42(&
      &41)+acd42(40)+acd42(39)+acd42(37)+acd42(38)
      acd42(34)=-ninjaP*acd42(34)
      acd42(37)=-acd42(14)*acd42(15)
      acd42(38)=-acd42(12)*acd42(13)
      acd42(34)=acd42(34)+acd42(37)+acd42(38)
      acd42(34)=acd42(3)*acd42(34)
      acd42(37)=-ninjaP-acd42(28)
      acd42(37)=acd42(10)*acd42(37)
      acd42(38)=acd42(27)*acd42(32)
      acd42(39)=acd42(24)*acd42(31)
      acd42(40)=acd42(20)*acd42(30)
      acd42(41)=acd42(17)*acd42(29)
      acd42(42)=acd42(15)*acd42(25)
      acd42(43)=acd42(14)*acd42(22)
      acd42(44)=acd42(13)*acd42(21)
      acd42(45)=acd42(12)*acd42(18)
      acd42(34)=acd42(34)+acd42(45)+acd42(44)+acd42(43)+acd42(42)+acd42(41)+acd&
      &42(40)+acd42(39)+acd42(33)+acd42(38)+acd42(37)
      brack(ninjaidxt1mu0)=acd42(36)
      brack(ninjaidxt1mu2)=0.0_ki
      brack(ninjaidxt0mu0)=acd42(34)
      brack(ninjaidxt0mu2)=acd42(35)
   end subroutine brack_32
!---#] subroutine brack_32:
!---#[ subroutine numerator_t3:
   subroutine numerator_t3(ncut, a, b, c, param, deg, coeffs) &
   & bind(c, name="p9_part21part21_part25part25part21_d42h0_qp_ninja_t3")
      use iso_c_binding, only: c_int
      use quadninjago_module, only: ki => ki_nin
      use p9_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_abbrevd42h0_qp
      implicit none
      integer(c_int), intent(in) :: ncut, deg
      complex(ki), dimension(0:3), intent(in) :: a, b, c
      complex(ki), intent(in) :: param
      complex(ki), dimension(4) :: vecA, vecB, vecC
      complex(ki), dimension(0:*), intent(out) :: coeffs
      integer :: t1
      complex(ki), dimension(4) :: qshift
      qshift = k2
	     vecA(1:4) = + a(0:3) - qshift(1:4)
	     vecB(1:4) = + b(0:3)
	     vecC(1:4) = + c(0:3)
      if (deg.lt.0) return
      t1 = 0
      call cond_t(epspow.eq.t1,brack_31,vecA,vecB,vecC,param,coeffs)
      if (deg.le.(1+(0))) return
      call cond_t(epspow.eq.t1,brack_32,vecA,vecB,vecC,param,coeffs)
   end subroutine numerator_t3
!---#] subroutine numerator_t3:
end module     p9_part21part21_part25part25part21_d42h0l131_qp
