module     p3_part1part21_part25part25part1_d11h2l1d_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p3_part1part2 &
   ! &1_part25part25part1/helicity2d11h2l1d_qp.f90
   ! generator: buildfortran_d.py
   use p3_part1part21_part25part25part1_config, only: ki => ki_qp
   use p3_part1part21_part25part25part1_util_qp, only: cond, d => metric_tensor
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
      use p3_part1part21_part25part25part1_model_qp
      use p3_part1part21_part25part25part1_kinematics_qp
      use p3_part1part21_part25part25part1_color_qp
      use p3_part1part21_part25part25part1_abbrevd11h2_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(17) :: acd11
      complex(ki) :: brack
      acd11(1)=dotproduct(qshift,qshift)
      acd11(2)=abb11(5)
      acd11(3)=dotproduct(qshift,spvak1k2)
      acd11(4)=abb11(7)
      acd11(5)=dotproduct(qshift,spvak2k5)
      acd11(6)=abb11(8)
      acd11(7)=dotproduct(qshift,spvak1e2)
      acd11(8)=dotproduct(qshift,spvae2k5)
      acd11(9)=abb11(4)
      acd11(10)=abb11(12)
      acd11(11)=abb11(9)
      acd11(12)=abb11(6)
      acd11(13)=acd11(9)*acd11(7)
      acd11(13)=-acd11(11)+acd11(13)
      acd11(13)=acd11(8)*acd11(13)
      acd11(14)=acd11(2)*acd11(1)
      acd11(15)=-acd11(4)*acd11(3)
      acd11(16)=-acd11(6)*acd11(5)
      acd11(17)=-acd11(10)*acd11(7)
      brack=acd11(12)+acd11(13)+acd11(14)+acd11(15)+acd11(16)+acd11(17)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p3_part1part21_part25part25part1_model_qp
      use p3_part1part21_part25part25part1_kinematics_qp
      use p3_part1part21_part25part25part1_color_qp
      use p3_part1part21_part25part25part1_abbrevd11h2_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(18) :: acd11
      complex(ki) :: brack
      acd11(1)=qshift(iv1)
      acd11(2)=abb11(5)
      acd11(3)=spvak1k2(iv1)
      acd11(4)=abb11(7)
      acd11(5)=spvak2k5(iv1)
      acd11(6)=abb11(8)
      acd11(7)=spvak1e2(iv1)
      acd11(8)=dotproduct(qshift,spvae2k5)
      acd11(9)=abb11(4)
      acd11(10)=abb11(12)
      acd11(11)=spvae2k5(iv1)
      acd11(12)=dotproduct(qshift,spvak1e2)
      acd11(13)=abb11(9)
      acd11(14)=acd11(8)*acd11(9)
      acd11(14)=-acd11(10)+acd11(14)
      acd11(14)=acd11(7)*acd11(14)
      acd11(15)=acd11(12)*acd11(9)
      acd11(15)=-acd11(13)+acd11(15)
      acd11(15)=acd11(11)*acd11(15)
      acd11(16)=acd11(2)*acd11(1)
      acd11(17)=-acd11(4)*acd11(3)
      acd11(18)=-acd11(6)*acd11(5)
      brack=acd11(14)+acd11(15)+2.0_ki*acd11(16)+acd11(17)+acd11(18)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p3_part1part21_part25part25part1_model_qp
      use p3_part1part21_part25part25part1_kinematics_qp
      use p3_part1part21_part25part25part1_color_qp
      use p3_part1part21_part25part25part1_abbrevd11h2_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(9) :: acd11
      complex(ki) :: brack
      acd11(1)=d(iv1,iv2)
      acd11(2)=abb11(5)
      acd11(3)=spvak1e2(iv1)
      acd11(4)=spvae2k5(iv2)
      acd11(5)=abb11(4)
      acd11(6)=spvak1e2(iv2)
      acd11(7)=spvae2k5(iv1)
      acd11(8)=acd11(4)*acd11(3)
      acd11(9)=acd11(7)*acd11(6)
      acd11(8)=acd11(9)+acd11(8)
      acd11(8)=acd11(5)*acd11(8)
      acd11(9)=acd11(2)*acd11(1)
      brack=acd11(8)+2.0_ki*acd11(9)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p3_part1part21_part25part25part1_model_qp
      use p3_part1part21_part25part25part1_kinematics_qp
      use p3_part1part21_part25part25part1_color_qp
      use p3_part1part21_part25part25part1_abbrevd11h2_qp
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd11
      complex(ki) :: brack
      brack=0.0_ki
   end function brack_4
!---#] function brack_4:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3) result(numerator)
      use p3_part1part21_part25part25part1_globalsl1_qp, only: epspow
      use p3_part1part21_part25part25part1_kinematics_qp
      use p3_part1part21_part25part25part1_abbrevd11h2_qp
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
      qshift = -k2
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
end module     p3_part1part21_part25part25part1_d11h2l1d_qp
