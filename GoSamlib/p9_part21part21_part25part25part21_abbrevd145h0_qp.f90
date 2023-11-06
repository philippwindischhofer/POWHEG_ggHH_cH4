module     p9_part21part21_part25part25part21_abbrevd145h0_qp
   use p9_part21part21_part25part25part21_config, only: ki => ki_qp
   use p9_part21part21_part25part25part21_kinematics_qp, only: epstensor
   use p9_part21part21_part25part25part21_globalsh0_qp
   implicit none
   private
   complex(ki), dimension(37), public :: abb145
   complex(ki), public :: R2d145
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p9_part21part21_part25part25part21_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p9_part21part21_part25part25part21_kinematics_qp
      use p9_part21part21_part25part25part21_model_qp
      use p9_part21part21_part25part25part21_color_qp, only: TR
      use p9_part21part21_part25part25part21_globalsl1_qp, only: epspow
      implicit none
      abb145(1)=1.0_ki/(-mdlMh**2+es34)
      abb145(2)=sqrt(mdlMT**2)
      abb145(3)=spbe2k1*spak1e5
      abb145(4)=spbk5e2*spae5k5
      abb145(5)=spak2e5*spbk2e2
      abb145(5)=abb145(3)+abb145(5)-abb145(4)
      abb145(5)=spbe5e1*spae1e2*abb145(5)
      abb145(6)=spae2k5*spbk5e5
      abb145(7)=spak1e2*spbe5k1
      abb145(6)=abb145(6)-abb145(7)
      abb145(7)=-spbe2e1*abb145(6)
      abb145(8)=spae2k2*spbe2e1
      abb145(9)=spbe5k2*abb145(8)
      abb145(7)=abb145(9)+abb145(7)
      abb145(7)=spae1e5*abb145(7)
      abb145(9)=spae2e5*spbe2e1
      abb145(10)=spae1k1*spbe5k1
      abb145(11)=abb145(10)*abb145(9)
      abb145(12)=spbk1e1*spak1e5
      abb145(13)=abb145(12)*spbe5e2
      abb145(14)=spae1e2*abb145(13)
      abb145(5)=abb145(14)+abb145(11)+abb145(7)+abb145(5)
      abb145(7)=c1-c2
      abb145(7)=abb145(1)*abb145(7)*mdlGC7**3*mdlGC43*mdlGC45
      abb145(11)=abb145(2)**3*abb145(7)
      abb145(5)=abb145(11)*abb145(5)
      abb145(14)=spbk5k2*spae5k5
      abb145(15)=-spbk2k1*spak1e5
      abb145(14)=abb145(15)+abb145(14)
      abb145(7)=abb145(2)*abb145(7)
      abb145(15)=abb145(10)*abb145(7)
      abb145(14)=abb145(15)*abb145(8)*abb145(14)
      abb145(16)=spak2k5*spbk5e5
      abb145(17)=spak1k2*spbe5k1
      abb145(16)=abb145(16)-abb145(17)
      abb145(17)=abb145(7)*spae1e2
      abb145(18)=abb145(16)*abb145(17)*spbk2e2
      abb145(19)=abb145(12)*abb145(18)
      abb145(5)=abb145(5)+abb145(19)+abb145(14)
      abb145(5)=1.0_ki/2.0_ki*abb145(5)
      abb145(13)=abb145(17)*abb145(13)
      abb145(14)=abb145(15)*abb145(9)
      abb145(13)=abb145(13)+abb145(14)
      abb145(13)=1.0_ki/2.0_ki*abb145(13)
      abb145(14)=-spbe5e1*abb145(11)*spae1e2
      abb145(17)=-spbe2e1*abb145(11)*spae1e5
      abb145(19)=abb145(11)*spbe5e2
      abb145(20)=spae1e2*abb145(19)
      abb145(18)=abb145(20)+abb145(18)
      abb145(18)=1.0_ki/2.0_ki*abb145(18)
      abb145(20)=1.0_ki/2.0_ki*abb145(7)
      abb145(21)=abb145(20)*spae1e2
      abb145(22)=spbe5e2*abb145(21)
      abb145(23)=abb145(7)*spae5k5
      abb145(24)=abb145(23)*spbk5k2
      abb145(25)=abb145(7)*spak1e5
      abb145(26)=abb145(25)*spbk2k1
      abb145(24)=abb145(24)-abb145(26)
      abb145(26)=abb145(8)*abb145(24)
      abb145(9)=abb145(11)*abb145(9)
      abb145(9)=abb145(9)+abb145(26)
      abb145(9)=1.0_ki/2.0_ki*abb145(9)
      abb145(26)=abb145(20)*spbe2e1
      abb145(27)=spae2e5*abb145(26)
      abb145(3)=abb145(3)-abb145(4)
      abb145(3)=abb145(15)*abb145(3)
      abb145(4)=-spae1e5*abb145(16)
      abb145(16)=spak2e5*abb145(10)
      abb145(4)=abb145(16)+abb145(4)
      abb145(16)=abb145(7)*spbk2e2
      abb145(4)=abb145(16)*abb145(4)
      abb145(19)=-spae1e5*abb145(19)
      abb145(3)=abb145(19)+abb145(4)+abb145(3)
      abb145(3)=1.0_ki/2.0_ki*abb145(3)
      abb145(4)=abb145(20)*spae1e5
      abb145(19)=-spbe5e2*abb145(4)
      abb145(23)=-spbk5e2*abb145(23)
      abb145(16)=spak2e5*abb145(16)
      abb145(25)=spbe2k1*abb145(25)
      abb145(16)=abb145(25)+abb145(23)+abb145(16)
      abb145(16)=1.0_ki/2.0_ki*abb145(16)
      abb145(23)=spbe5k2*spae2k2
      abb145(6)=abb145(6)-abb145(23)
      abb145(23)=abb145(12)*abb145(7)
      abb145(25)=-abb145(23)*abb145(6)
      abb145(24)=-spae2k2*abb145(24)
      abb145(11)=-spae2e5*abb145(11)
      abb145(11)=abb145(11)+abb145(24)
      abb145(11)=spbe5e1*abb145(11)
      abb145(11)=abb145(25)+abb145(11)
      abb145(11)=1.0_ki/2.0_ki*abb145(11)
      abb145(24)=abb145(20)*spbe5e1
      abb145(25)=-spae2e5*abb145(24)
      abb145(6)=-1.0_ki/2.0_ki*abb145(7)*abb145(6)
      abb145(28)=abb145(26)*abb145(10)
      abb145(29)=-spae5k5*abb145(28)
      abb145(30)=abb145(26)*spae5k5
      abb145(31)=spae5k5*abb145(24)
      abb145(12)=abb145(21)*abb145(12)
      abb145(32)=-spbk5e5*abb145(12)
      abb145(33)=-spbk5e5*abb145(21)
      abb145(34)=spbk5e5*abb145(4)
      abb145(8)=abb145(20)*abb145(8)
      abb145(10)=abb145(10)*abb145(8)
      abb145(20)=-spae2k2*abb145(24)
      abb145(35)=spbk2e2*abb145(12)
      abb145(36)=spbk2e2*abb145(21)
      abb145(37)=-spbk2e2*abb145(4)
      abb145(28)=-spak1e5*abb145(28)
      abb145(26)=-spak1e5*abb145(26)
      abb145(24)=spak1e5*abb145(24)
      abb145(12)=-spbe5k1*abb145(12)
      abb145(21)=-spbe5k1*abb145(21)
      abb145(4)=spbe5k1*abb145(4)
      R2d145=0.0_ki
      rat2 = rat2 + R2d145
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='145' value='", &
          & R2d145, "'/>"
      end if
   end subroutine
end module p9_part21part21_part25part25part21_abbrevd145h0_qp
