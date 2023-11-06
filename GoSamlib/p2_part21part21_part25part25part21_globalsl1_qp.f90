module     p2_part21part21_part25part25part21_globalsl1_qp
   use p2_part21part21_part25part25part21_config, only: ki => ki_qp
   use p2_part21part21_part25part25part21_color_qp, only: numcs
   implicit none

   private

   
   ! col0 is the color index to be returned in the virtual diagrams
   integer, public :: col0
   integer, dimension(numcs), public :: perm
   logical, public :: use_perm

   integer, public :: epspow

   interface ccontract
      module procedure ccontract_cc
      module procedure ccontract_rc
   end interface

   public :: ccontract

contains
   !---#[ function ccontract_cc:
   pure function ccontract_cc(color_vector1, color_vector2) result(amp)
      use p2_part21part21_part25part25part21_color_qp, only: cmat => CC
      implicit none
      complex(ki), dimension(numcs), intent(in) :: color_vector1
      complex(ki), dimension(numcs), intent(in) :: color_vector2
      complex(ki) :: amp
      complex(ki), dimension(numcs) :: v1, v2

      if (use_perm) then
         v1 = matmul(cmat, color_vector1(perm))
      else
         v1 = matmul(cmat, color_vector1)
      end if
      v2 = conjg(color_vector2)
      amp = sum(v1(:) * v2(:))
   end function  ccontract_cc
   !---#] function ccontract_cc:
   !---#[ function ccontract_rc:
   pure function ccontract_rc(color_vector1, color_vector2) result(amp)
      use p2_part21part21_part25part25part21_color_qp, only: cmat => CC
      implicit none
      real(ki), dimension(numcs), intent(in) :: color_vector1
      complex(ki), dimension(numcs), intent(in) :: color_vector2
      complex(ki) :: amp
      complex(ki), dimension(numcs) :: v1, v2

      if (use_perm) then
         v1 = matmul(cmat, color_vector1(perm))
      else
         v1 = matmul(cmat, color_vector1)
      end if
      v2 = conjg(color_vector2)
      amp = sum(v1(:) * v2(:))
   end function  ccontract_rc
   !---#] function ccontract_rc:
end module p2_part21part21_part25part25part21_globalsl1_qp
