module     p0_part21part21_part25part25_abbrevd7h0
   use p0_part21part21_part25part25_config, only: ki
   use p0_part21part21_part25part25_kinematics, only: epstensor
   use p0_part21part21_part25part25_globalsh0
   implicit none
   private
   complex(ki), dimension(12), public :: abb7
   complex(ki), public :: R2d7
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
      abb7(1)=1.0_ki/(-mdlMh**2+es12)
      abb7(2)=sqrt(mdlMT**2)
      abb7(3)=mdlGC7**2*mdlGC45*mdlGC43*c1*i_
      abb7(4)=abb7(1)*abb7(2)
      abb7(5)=abb7(4)*abb7(3)
      abb7(6)=abb7(5)*spbe2e1
      abb7(7)=spae1e2*abb7(6)
      abb7(8)=spae1e2*spbk2e2
      abb7(9)=spbk1e1*spak1k2*abb7(8)
      abb7(10)=spbe2e1*spae2k2
      abb7(11)=spae1k1*spbk2k1*abb7(10)
      abb7(12)=abb7(2)**2*spbe2e1*spae1e2
      abb7(9)=-2.0_ki*abb7(12)+abb7(9)+abb7(11)
      abb7(3)=abb7(4)*abb7(3)*abb7(9)
      abb7(4)=-2.0_ki*abb7(7)
      abb7(9)=-spae1k1*spbe2k1
      abb7(11)=-spbk2e2*spae1k2
      abb7(9)=abb7(9)+abb7(11)
      abb7(9)=abb7(9)*abb7(5)
      abb7(11)=-spbk1e1*spak1e2
      abb7(12)=-spae2k2*spbk2e1
      abb7(11)=abb7(11)+abb7(12)
      abb7(11)=abb7(11)*abb7(5)
      abb7(12)=-4.0_ki*abb7(5)
      abb7(10)=abb7(5)*abb7(10)
      abb7(8)=abb7(5)*abb7(8)
      abb7(6)=spae1k1*abb7(6)
      abb7(5)=spbk1e1*abb7(5)*spae1e2
      R2d7=-abb7(7)
      rat2 = rat2 + R2d7
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='7' value='", &
          & R2d7, "'/>"
      end if
   end subroutine
end module p0_part21part21_part25part25_abbrevd7h0
