module     p9_part21part21_part25part25part21_d42h0l1d_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part &
   ! &21_part25part25part21/helicity0d42h0l1d_qp.f90
   ! generator: buildfortran_d.py
   use p9_part21part21_part25part25part21_config, only: ki => ki_qp
   use p9_part21part21_part25part25part21_util_qp, only: cond, d => metric_tens&
   &or
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   integer, private :: iv3
   real(ki), dimension(4), private :: qshift
   public :: derivative
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_color_qp
      use p9_part21part21_part25part25part21_abbrevd42h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(29) :: acd42
      complex(ki) :: brack
      acd42(1)=dotproduct(qshift,qshift)
      acd42(2)=abb42(9)
      acd42(3)=dotproduct(qshift,spvak1k2)
      acd42(4)=abb42(5)
      acd42(5)=dotproduct(qshift,spvak1k5)
      acd42(6)=abb42(7)
      acd42(7)=dotproduct(qshift,spvak2k1)
      acd42(8)=abb42(12)
      acd42(9)=dotproduct(qshift,spvak5k1)
      acd42(10)=abb42(15)
      acd42(11)=dotproduct(qshift,spvak2e1)
      acd42(12)=dotproduct(qshift,spvae1k2)
      acd42(13)=abb42(13)
      acd42(14)=abb42(14)
      acd42(15)=abb42(11)
      acd42(16)=dotproduct(qshift,spvak5e1)
      acd42(17)=dotproduct(qshift,spvae1k5)
      acd42(18)=abb42(8)
      acd42(19)=abb42(6)
      acd42(20)=abb42(4)
      acd42(21)=-acd42(13)*acd42(11)
      acd42(21)=-acd42(15)+acd42(21)
      acd42(21)=acd42(12)*acd42(21)
      acd42(22)=-acd42(16)*acd42(13)
      acd42(22)=-acd42(19)+acd42(22)
      acd42(22)=acd42(17)*acd42(22)
      acd42(23)=-acd42(2)*acd42(1)
      acd42(24)=-acd42(4)*acd42(3)
      acd42(25)=-acd42(6)*acd42(5)
      acd42(26)=-acd42(8)*acd42(7)
      acd42(27)=-acd42(10)*acd42(9)
      acd42(28)=-acd42(14)*acd42(11)
      acd42(29)=-acd42(18)*acd42(16)
      brack=acd42(20)+acd42(21)+acd42(22)+acd42(23)+acd42(24)+acd42(25)+acd42(2&
      &6)+acd42(27)+acd42(28)+acd42(29)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_color_qp
      use p9_part21part21_part25part25part21_abbrevd42h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(33) :: acd42
      complex(ki) :: brack
      acd42(1)=qshift(iv1)
      acd42(2)=abb42(9)
      acd42(3)=spvak1k2(iv1)
      acd42(4)=abb42(5)
      acd42(5)=spvak1k5(iv1)
      acd42(6)=abb42(7)
      acd42(7)=spvak2k1(iv1)
      acd42(8)=abb42(12)
      acd42(9)=spvak5k1(iv1)
      acd42(10)=abb42(15)
      acd42(11)=spvak2e1(iv1)
      acd42(12)=dotproduct(qshift,spvae1k2)
      acd42(13)=abb42(13)
      acd42(14)=abb42(14)
      acd42(15)=spvae1k2(iv1)
      acd42(16)=dotproduct(qshift,spvak2e1)
      acd42(17)=abb42(11)
      acd42(18)=spvak5e1(iv1)
      acd42(19)=dotproduct(qshift,spvae1k5)
      acd42(20)=abb42(8)
      acd42(21)=spvae1k5(iv1)
      acd42(22)=dotproduct(qshift,spvak5e1)
      acd42(23)=abb42(6)
      acd42(24)=acd42(12)*acd42(11)
      acd42(25)=acd42(16)*acd42(15)
      acd42(26)=acd42(19)*acd42(18)
      acd42(27)=acd42(22)*acd42(21)
      acd42(24)=acd42(27)+acd42(26)+acd42(25)+acd42(24)
      acd42(24)=acd42(13)*acd42(24)
      acd42(25)=acd42(2)*acd42(1)
      acd42(26)=acd42(4)*acd42(3)
      acd42(27)=acd42(6)*acd42(5)
      acd42(28)=acd42(8)*acd42(7)
      acd42(29)=acd42(10)*acd42(9)
      acd42(30)=acd42(14)*acd42(11)
      acd42(31)=acd42(17)*acd42(15)
      acd42(32)=acd42(20)*acd42(18)
      acd42(33)=acd42(23)*acd42(21)
      brack=acd42(24)+2.0_ki*acd42(25)+acd42(26)+acd42(27)+acd42(28)+acd42(29)+&
      &acd42(30)+acd42(31)+acd42(32)+acd42(33)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_color_qp
      use p9_part21part21_part25part25part21_abbrevd42h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(15) :: acd42
      complex(ki) :: brack
      acd42(1)=d(iv1,iv2)
      acd42(2)=abb42(9)
      acd42(3)=spvak2e1(iv1)
      acd42(4)=spvae1k2(iv2)
      acd42(5)=abb42(13)
      acd42(6)=spvak2e1(iv2)
      acd42(7)=spvae1k2(iv1)
      acd42(8)=spvak5e1(iv1)
      acd42(9)=spvae1k5(iv2)
      acd42(10)=spvak5e1(iv2)
      acd42(11)=spvae1k5(iv1)
      acd42(12)=-acd42(4)*acd42(3)
      acd42(13)=-acd42(7)*acd42(6)
      acd42(14)=-acd42(9)*acd42(8)
      acd42(15)=-acd42(11)*acd42(10)
      acd42(12)=acd42(15)+acd42(14)+acd42(13)+acd42(12)
      acd42(12)=acd42(5)*acd42(12)
      acd42(13)=acd42(2)*acd42(1)
      brack=acd42(12)-2.0_ki*acd42(13)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_color_qp
      use p9_part21part21_part25part25part21_abbrevd42h0_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd42
      complex(ki) :: brack
      brack=0.0_ki
   end function brack_4
!---#] function brack_4:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3) result(numerator)
      use p9_part21part21_part25part25part21_globalsl1_qp, only: epspow
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_abbrevd42h0_qp
      implicit none
      complex(ki), intent(in) :: mu2
      integer, intent(in), optional :: i1
      integer, intent(in), optional :: i2
      integer, intent(in), optional :: i3
      complex(ki) :: numerator
      complex(ki) :: loc
      integer :: t1
      integer :: deg
      complex(ki), dimension(4), parameter :: Q = (/ (0.0_ki,0.0_ki),(0.0_ki,0.&
      &0_ki),(0.0_ki,0.0_ki),(0.0_ki,0.0_ki)/)
      qshift = k2
      numerator = 0.0_ki
      deg = 0
      if(present(i1)) then
          iv1=i1
          deg=1
      else
          iv1=1
      end if
      if(present(i2)) then
          iv2=i2
          deg=2
      else
          iv2=1
      end if
      if(present(i3)) then
          iv3=i3
          deg=3
      else
          iv3=1
      end if
      t1 = 0
      if(deg.eq.0) then
         numerator = cond(epspow.eq.t1,brack_1,Q,mu2)
         return
      end if
      if(deg.eq.1) then
         numerator = cond(epspow.eq.t1,brack_2,Q,mu2)
         return
      end if
      if(deg.eq.2) then
         numerator = cond(epspow.eq.t1,brack_3,Q,mu2)
         return
      end if
      if(deg.eq.3) then
         numerator = cond(epspow.eq.t1,brack_4,Q,mu2)
         return
      end if
   end function derivative
!---#] function derivative:
end module     p9_part21part21_part25part25part21_d42h0l1d_qp
