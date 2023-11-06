module     p0_part21part21_part25part25_abbrevd2h0
   use p0_part21part21_part25part25_config, only: ki
   use p0_part21part21_part25part25_kinematics, only: epstensor
   use p0_part21part21_part25part25_globalsh0
   implicit none
   private
   complex(ki), dimension(11), public :: abb2
   complex(ki), public :: R2d2
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p0_part21part21_part25part25_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p0_part21part21_part25part25_kinematics
      use p0_part21part21_part25part25_model
      use p0_part21part21_part25part25_color, only: TR
      use p0_part21part21_part25part25_globalsl1, only: epspow
      implicit none
      abb2(1)=sqrt(mdlMT**2)
      abb2(2)=mdlGC7**2*mdlGC46*c1
      abb2(3)=abb2(2)*abb2(1)*spbe2e1
      abb2(4)=spae1e2*abb2(3)
      abb2(5)=spae2k2*spae1k1*spbk2k1
      abb2(6)=spae1e2*abb2(1)**2
      abb2(5)=-2.0_ki*abb2(6)+abb2(5)
      abb2(5)=spbe2e1*abb2(5)
      abb2(6)=spae1e2*spbk2e2*spbk1e1*spak1k2
      abb2(5)=abb2(6)+abb2(5)
      abb2(2)=abb2(2)*abb2(1)
      abb2(5)=abb2(5)*abb2(2)
      abb2(6)=-2.0_ki*abb2(4)
      abb2(7)=-spae1k1*spbe2k1
      abb2(8)=-spbk2e2*spae1k2
      abb2(7)=abb2(7)+abb2(8)
      abb2(7)=abb2(7)*abb2(2)
      abb2(8)=-spbk1e1*spak1e2
      abb2(9)=-spae2k2*spbk2e1
      abb2(8)=abb2(8)+abb2(9)
      abb2(8)=abb2(8)*abb2(2)
      abb2(9)=-4.0_ki*abb2(2)
      abb2(10)=spae2k2*abb2(3)
      abb2(2)=spae1e2*abb2(2)
      abb2(11)=spbk2e2*abb2(2)
      abb2(3)=spae1k1*abb2(3)
      abb2(2)=spbk1e1*abb2(2)
      R2d2=-abb2(4)
      rat2 = rat2 + R2d2
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='2' value='", &
          & R2d2, "'/>"
      end if
   end subroutine
end module p0_part21part21_part25part25_abbrevd2h0
