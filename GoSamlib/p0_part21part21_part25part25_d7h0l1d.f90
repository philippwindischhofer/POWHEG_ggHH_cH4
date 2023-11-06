module     p0_part21part21_part25part25_d7h0l1d
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p0_part21part &
   ! &21_part25part25/helicity0d7h0l1d.f90
   ! generator: buildfortran_d.py
   use p0_part21part21_part25part25_config, only: ki
   use p0_part21part21_part25part25_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   integer, private :: iv3
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d7
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p0_part21part21_part25part25_model
      use p0_part21part21_part25part25_kinematics
      use p0_part21part21_part25part25_color
      use p0_part21part21_part25part25_abbrevd7h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(23) :: acd7
      complex(ki) :: brack
      acd7(1)=dotproduct(qshift,qshift)
      acd7(2)=abb7(4)
      acd7(3)=dotproduct(qshift,spvak1e2)
      acd7(4)=abb7(5)
      acd7(5)=dotproduct(qshift,spvae2k1)
      acd7(6)=abb7(6)
      acd7(7)=dotproduct(qshift,spvak2e1)
      acd7(8)=abb7(8)
      acd7(9)=dotproduct(qshift,spvae1k2)
      acd7(10)=abb7(10)
      acd7(11)=dotproduct(qshift,spvae1e2)
      acd7(12)=dotproduct(qshift,spvae2e1)
      acd7(13)=abb7(12)
      acd7(14)=abb7(11)
      acd7(15)=abb7(9)
      acd7(16)=abb7(3)
      acd7(17)=acd7(13)*acd7(11)
      acd7(17)=-acd7(15)+acd7(17)
      acd7(17)=acd7(12)*acd7(17)
      acd7(18)=acd7(2)*acd7(1)
      acd7(19)=-acd7(4)*acd7(3)
      acd7(20)=-acd7(6)*acd7(5)
      acd7(21)=-acd7(8)*acd7(7)
      acd7(22)=-acd7(10)*acd7(9)
      acd7(23)=-acd7(14)*acd7(11)
      brack=acd7(16)+acd7(17)+acd7(18)+acd7(19)+acd7(20)+acd7(21)+acd7(22)+acd7&
      &(23)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p0_part21part21_part25part25_model
      use p0_part21part21_part25part25_kinematics
      use p0_part21part21_part25part25_color
      use p0_part21part21_part25part25_abbrevd7h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(24) :: acd7
      complex(ki) :: brack
      acd7(1)=qshift(iv1)
      acd7(2)=abb7(4)
      acd7(3)=spvak1e2(iv1)
      acd7(4)=abb7(5)
      acd7(5)=spvae2k1(iv1)
      acd7(6)=abb7(6)
      acd7(7)=spvak2e1(iv1)
      acd7(8)=abb7(8)
      acd7(9)=spvae1k2(iv1)
      acd7(10)=abb7(10)
      acd7(11)=spvae1e2(iv1)
      acd7(12)=dotproduct(qshift,spvae2e1)
      acd7(13)=abb7(12)
      acd7(14)=abb7(11)
      acd7(15)=spvae2e1(iv1)
      acd7(16)=dotproduct(qshift,spvae1e2)
      acd7(17)=abb7(9)
      acd7(18)=-acd7(12)*acd7(13)
      acd7(18)=acd7(14)+acd7(18)
      acd7(18)=acd7(11)*acd7(18)
      acd7(19)=-acd7(16)*acd7(13)
      acd7(19)=acd7(17)+acd7(19)
      acd7(19)=acd7(15)*acd7(19)
      acd7(20)=acd7(2)*acd7(1)
      acd7(21)=acd7(4)*acd7(3)
      acd7(22)=acd7(6)*acd7(5)
      acd7(23)=acd7(8)*acd7(7)
      acd7(24)=acd7(10)*acd7(9)
      brack=acd7(18)+acd7(19)-2.0_ki*acd7(20)+acd7(21)+acd7(22)+acd7(23)+acd7(2&
      &4)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p0_part21part21_part25part25_model
      use p0_part21part21_part25part25_kinematics
      use p0_part21part21_part25part25_color
      use p0_part21part21_part25part25_abbrevd7h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(9) :: acd7
      complex(ki) :: brack
      acd7(1)=d(iv1,iv2)
      acd7(2)=abb7(4)
      acd7(3)=spvae1e2(iv1)
      acd7(4)=spvae2e1(iv2)
      acd7(5)=abb7(12)
      acd7(6)=spvae1e2(iv2)
      acd7(7)=spvae2e1(iv1)
      acd7(8)=acd7(4)*acd7(3)
      acd7(9)=acd7(7)*acd7(6)
      acd7(8)=acd7(9)+acd7(8)
      acd7(8)=acd7(5)*acd7(8)
      acd7(9)=acd7(2)*acd7(1)
      brack=acd7(8)+2.0_ki*acd7(9)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p0_part21part21_part25part25_model
      use p0_part21part21_part25part25_kinematics
      use p0_part21part21_part25part25_color
      use p0_part21part21_part25part25_abbrevd7h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd7
      complex(ki) :: brack
      brack=0.0_ki
   end function brack_4
!---#] function brack_4:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3) result(numerator)
      use p0_part21part21_part25part25_globalsl1, only: epspow
      use p0_part21part21_part25part25_kinematics
      use p0_part21part21_part25part25_abbrevd7h0
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
!---#[ subroutine reconstruct_d7:
   subroutine     reconstruct_d7(coeffs)
      use p0_part21part21_part25part25_groups, only: tensrec_info_group1
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      complex(ki), parameter :: ctwo = (2.0_ki, 0.0_ki)
      type(tensrec_info_group1), intent(out) :: coeffs
      ! rank 3 case :
      !---[# reconstruct coeffs%coeffs_7:
      coeffs%coeffs_7%c0 = derivative(czip)
      coeffs%coeffs_7%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_7%c1(1,2) = derivative(czip,1,1)/2.0_ki
      coeffs%coeffs_7%c1(1,3) = derivative(czip,1,1,1)/6.0_ki
      coeffs%coeffs_7%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_7%c1(2,2) = derivative(czip,2,2)/2.0_ki
      coeffs%coeffs_7%c1(2,3) = -derivative(czip,2,2,2)/6.0_ki
      coeffs%coeffs_7%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_7%c1(3,2) = derivative(czip,3,3)/2.0_ki
      coeffs%coeffs_7%c1(3,3) = -derivative(czip,3,3,3)/6.0_ki
      coeffs%coeffs_7%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_7%c1(4,2) = derivative(czip,4,4)/2.0_ki
      coeffs%coeffs_7%c1(4,3) = -derivative(czip,4,4,4)/6.0_ki
      coeffs%coeffs_7%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_7%c2(1,2) = derivative(czip,1,2,2)/2.0_ki
      coeffs%coeffs_7%c2(1,3) = -derivative(czip,1,1,2)/2.0_ki
      coeffs%coeffs_7%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_7%c2(2,2) = derivative(czip,1,3,3)/2.0_ki
      coeffs%coeffs_7%c2(2,3) = -derivative(czip,1,1,3)/2.0_ki
      coeffs%coeffs_7%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_7%c2(3,2) = derivative(czip,1,4,4)/2.0_ki
      coeffs%coeffs_7%c2(3,3) = -derivative(czip,1,1,4)/2.0_ki
      coeffs%coeffs_7%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_7%c2(4,2) = -derivative(czip,2,3,3)/2.0_ki
      coeffs%coeffs_7%c2(4,3) = -derivative(czip,2,2,3)/2.0_ki
      coeffs%coeffs_7%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_7%c2(5,2) = -derivative(czip,2,4,4)/2.0_ki
      coeffs%coeffs_7%c2(5,3) = -derivative(czip,2,2,4)/2.0_ki
      coeffs%coeffs_7%c2(6,1) = derivative(czip,3,4)
      coeffs%coeffs_7%c2(6,2) = -derivative(czip,3,4,4)/2.0_ki
      coeffs%coeffs_7%c2(6,3) = -derivative(czip,3,3,4)/2.0_ki
      coeffs%coeffs_7%c3(1,1) = derivative(czip,1,2,3)
      coeffs%coeffs_7%c3(2,1) = derivative(czip,1,2,4)
      coeffs%coeffs_7%c3(3,1) = derivative(czip,1,3,4)
      coeffs%coeffs_7%c3(4,1) = -derivative(czip,2,3,4)
      !---#] reconstruct coeffs%coeffs_7:
   end subroutine reconstruct_d7
!---#] subroutine reconstruct_d7:
end module     p0_part21part21_part25part25_d7h0l1d
