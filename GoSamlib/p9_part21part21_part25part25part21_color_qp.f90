module     p9_part21part21_part25part25part21_color_qp
   ! file: /draco/u/lscyboz/POWHEG-BOX-V2/ggHH_EWChL/GoSam_POWHEG/p9_part21part21_part25part25part21/common/color.f90 
   ! generator: buildcolor.py
   use p9_part21part21_part25part25part21_config, only: ki => ki_qp
   use p9_part21part21_part25part25part21_model_qp, only: NC, Nf
   implicit none
   private
   save

   public :: init_color

   real(ki), parameter, public :: TR = 0.5_ki

   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   real(ki), parameter :: pi = &
   & 3.1415926535897932384626433832795028841971693993751058209749445920_ki
   real(ki), parameter :: pi6 = pi*pi/6.0_ki

   integer, parameter, public :: numcs = 2
   complex(ki), dimension(numcs,numcs), public :: T1T1
   complex(ki), dimension(numcs,numcs), public :: T1T2
   complex(ki), dimension(numcs,numcs), public :: T1T5
   complex(ki), dimension(numcs,numcs), public :: T2T2
   complex(ki), dimension(numcs,numcs), public :: T2T5
   complex(ki), dimension(numcs,numcs), public :: T5T5
   complex(ki), dimension(numcs,numcs), public :: CC
   real(ki), public :: incolors

   real(ki), public :: CA, CF, KA, KF, gammaA, gammaF

   real(ki) :: NA
   real(ki), dimension(8) :: cabb

   ! Basis vectors
   real(ki), dimension(numcs), parameter, public :: c1 = &
      & (/         1.0_ki         , 0.0_ki/)
   real(ki), dimension(numcs), parameter, public :: c2 = &
      & (/         0.0_ki,          1.0_ki/)
contains
   subroutine     init_color()
      implicit none
      NA=NC*NC-1.0_ki
      incolors=1.0_ki*NA*NA
      cabb(1)=NC**(-1)
      cabb(2)=NC**2
      cabb(3)=cabb(2)-3.0_ki
      cabb(4)=cabb(3)*NC
      cabb(4)=cabb(4)+2.0_ki*cabb(1)
      cabb(5)=TR**3
      cabb(4)=cabb(4)*cabb(5)
      cabb(6)=2.0_ki*cabb(5)
      cabb(7)=NC-cabb(1)
      cabb(7)=cabb(6)*cabb(7)
      cabb(3)=cabb(3)*cabb(2)
      cabb(3)=cabb(3)+2.0_ki
      cabb(5)=cabb(3)*cabb(5)
      cabb(2)=cabb(2)-1.0_ki
      cabb(6)=cabb(6)*cabb(2)
      cabb(8)=TR**4
      cabb(3)=cabb(3)*cabb(8)
      cabb(8)=2.0_ki*cabb(8)
      cabb(2)=cabb(8)*cabb(2)
      T1T1(1,1)=cabb(5)
      T1T1(1,2)=-cabb(6)
      T1T1(2,1)=-cabb(6)
      T1T1(2,1)=-cabb(6)
      T1T1(1,2)=-cabb(6)
      T1T1(2,2)=cabb(5)
      T1T2(1,1)=-cabb(3)
      T1T2(1,2)=cabb(2)
      T1T2(2,1)=cabb(2)
      T1T2(2,1)=cabb(2)
      T1T2(1,2)=cabb(2)
      T1T2(2,2)=-cabb(3)
      T1T5(1,1)=-cabb(3)
      T1T5(1,2)=cabb(2)
      T1T5(2,1)=cabb(2)
      T1T5(2,1)=cabb(2)
      T1T5(1,2)=cabb(2)
      T1T5(2,2)=-cabb(3)
      T2T2(1,1)=cabb(5)
      T2T2(1,2)=-cabb(6)
      T2T2(2,1)=-cabb(6)
      T2T2(2,1)=-cabb(6)
      T2T2(1,2)=-cabb(6)
      T2T2(2,2)=cabb(5)
      T2T5(1,1)=-cabb(3)
      T2T5(1,2)=cabb(2)
      T2T5(2,1)=cabb(2)
      T2T5(2,1)=cabb(2)
      T2T5(1,2)=cabb(2)
      T2T5(2,2)=-cabb(3)
      T5T5(1,1)=cabb(5)
      T5T5(1,2)=-cabb(6)
      T5T5(2,1)=-cabb(6)
      T5T5(2,1)=-cabb(6)
      T5T5(1,2)=-cabb(6)
      T5T5(2,2)=cabb(5)
      CC(1,1)=cabb(4)
      CC(1,2)=-cabb(7)
      CC(2,1)=-cabb(7)
      CC(2,2)=cabb(4)

      CA = NC
      CF = TR * NA / NC
      ! KA = Kg in (C.11) [Catani,Seymour]
      KA = (67.0_ki/18.0_ki - pi6) * CA &
         & - 10.0_ki/9.0_ki * TR * Nf
      ! KF = Kq in (C.11) [Catani,Seymour]
      KF = (3.5_ki - pi6) * CF
      ! gammaA = \gamma_g in (C.11) [Catani,Seymour]
      gammaA = 11.0_ki/6.0_ki * CA - 2.0_ki/3.0_ki * TR * Nf
      ! gammaF = \gamma_q in (C.11) [Catani,Seymour]
      gammaF = 1.5_ki * CF
   end subroutine ! end of "init_color"
   subroutine     inspect_color(unit)
      implicit none
      integer, intent(in) :: unit
      integer :: i, j
      character :: ch1, ch2, ch3

      ch3 = ","
      write (unit,'(A13)') "gosam_color=["
      do i=1,numcs
         do j=1,numcs
            if (j==1) then
               ch1 = "["
            else
               ch1 = " "
            endif

            if (j == numcs) then
               ch2 = "]"
               if (i == numcs) then
                  ch3 = "]"
               end if
            else
               ch2 = ","
            end if

            if (j == numcs) then
               write (unit,'(3x,A1,A8,G23.16,A1,G23.16,A1,A1,A1)') &
               & ch1, "complex(", real(CC(i,j)), ",", aimag(CC(i,j)), ")", &
               & ch2, ch3
            else
               write (unit,'(3x,A1,A8,G23.16,A1,G23.16,A1,A1)') &
               & ch1, "complex(", real(CC(i,j)), ",", aimag(CC(i,j)), ")", &
               & ch2
            end if
         enddo
      enddo
   end subroutine inspect_color
end module p9_part21part21_part25part25part21_color_qp
