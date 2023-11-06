module     p10_part1part21_part25part25part1_abbrevd11h2_qp
   use p10_part1part21_part25part25part1_config, only: ki => ki_qp
   use p10_part1part21_part25part25part1_kinematics_qp, only: epstensor
   use p10_part1part21_part25part25part1_globalsh2_qp
   implicit none
   private
   complex(ki), dimension(14), public :: abb11
   complex(ki), public :: R2d11
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p10_part1part21_part25part25part1_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p10_part1part21_part25part25part1_kinematics_qp
      use p10_part1part21_part25part25part1_model_qp
      use p10_part1part21_part25part25part1_color_qp, only: TR
      use p10_part1part21_part25part25part1_globalsl1_qp, only: epspow
      implicit none
      abb11(1)=sqrt(mdlMT**2)
      abb11(2)=es51**(-1)
      abb11(3)=1.0_ki/(-mdlMh**2+es34)
      abb11(4)=mdlGC46*i_
      abb11(5)=mdlGC43*abb11(3)*mdlGC45
      abb11(4)=abb11(4)-abb11(5)
      abb11(5)=abb11(4)*spbk5e2*spak1e2
      abb11(6)=mdlGC7**3*abb11(2)*c1
      abb11(7)=abb11(6)*TR
      abb11(8)=abb11(1)*abb11(7)
      abb11(9)=2.0_ki*abb11(8)
      abb11(10)=abb11(5)*abb11(9)
      abb11(11)=spak1e2*spbk2e2
      abb11(12)=spbk5k1*abb11(11)*spak1k2
      abb11(13)=spbk5e2*spae2k2
      abb11(14)=spak1k5*abb11(13)*spbk5k2
      abb11(12)=abb11(12)+abb11(14)
      abb11(7)=-abb11(7)*abb11(12)*abb11(4)
      abb11(6)=2.0_ki*abb11(6)
      abb11(6)=TR*abb11(5)*abb11(1)**2*abb11(6)
      abb11(6)=abb11(7)+abb11(6)
      abb11(6)=2.0_ki*abb11(1)*abb11(6)
      abb11(5)=4.0_ki*abb11(5)*abb11(8)
      abb11(7)=-abb11(4)*abb11(9)
      abb11(9)=-spbk2e2*spak1k2*abb11(7)
      abb11(12)=spae2k2*spbk5k2*abb11(7)
      abb11(4)=-8.0_ki*abb11(4)*abb11(8)
      abb11(8)=abb11(11)*abb11(7)
      abb11(7)=-abb11(13)*abb11(7)
      R2d11=abb11(10)
      rat2 = rat2 + R2d11
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='11' value='", &
          & R2d11, "'/>"
      end if
   end subroutine
end module p10_part1part21_part25part25part1_abbrevd11h2_qp
