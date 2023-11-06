module     p2_part21part21_part25part25part21_abbrevd56h0_qp
   use p2_part21part21_part25part25part21_config, only: ki => ki_qp
   use p2_part21part21_part25part25part21_kinematics_qp, only: epstensor
   use p2_part21part21_part25part25part21_globalsh0_qp
   implicit none
   private
   complex(ki), dimension(35), public :: abb56
   complex(ki), public :: R2d56
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p2_part21part21_part25part25part21_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p2_part21part21_part25part25part21_kinematics_qp
      use p2_part21part21_part25part25part21_model_qp
      use p2_part21part21_part25part25part21_color_qp, only: TR
      use p2_part21part21_part25part25part21_globalsl1_qp, only: epspow
      implicit none
      abb56(1)=sqrt(mdlMT**2)
      abb56(2)=spbk2e2*spae2e5
      abb56(3)=spae1k2*abb56(2)
      abb56(4)=spae1k5*spbk5e2
      abb56(5)=-spae2e5*abb56(4)
      abb56(3)=abb56(5)+abb56(3)
      abb56(3)=spbe5e1*abb56(3)
      abb56(5)=spbk2e1*spae2k2
      abb56(6)=spbk5e1*spae2k5
      abb56(7)=-abb56(6)+abb56(5)
      abb56(8)=spbe5e2*spae1e5
      abb56(7)=abb56(8)*abb56(7)
      abb56(9)=spae1e2*spbe5e1
      abb56(10)=spae5k5*spbk5e2
      abb56(11)=-abb56(10)*abb56(9)
      abb56(12)=spbk5e5*spae2k5
      abb56(13)=abb56(12)*spbe2e1*spae1e5
      abb56(14)=spae1k1*spbe5e1
      abb56(15)=spbe2k1*spae2e5*abb56(14)
      abb56(16)=spbk1e1*spae1e5
      abb56(17)=spak1e2*spbe5e2*abb56(16)
      abb56(3)=abb56(17)+abb56(15)-abb56(13)+abb56(11)+abb56(7)+abb56(3)
      abb56(7)=c1-c2
      abb56(7)=abb56(7)*mdlGC7**3*i_*mdlGC46
      abb56(11)=-abb56(1)**3*abb56(7)
      abb56(3)=abb56(11)*abb56(3)
      abb56(7)=abb56(1)*abb56(7)
      abb56(15)=abb56(7)*spae2k2
      abb56(17)=spbk2k1*abb56(15)*abb56(14)
      abb56(18)=-abb56(10)*abb56(17)
      abb56(19)=abb56(7)*spbk2e2
      abb56(20)=spak1k2*abb56(19)*abb56(16)
      abb56(21)=-abb56(12)*abb56(20)
      abb56(3)=abb56(21)+abb56(18)+abb56(3)
      abb56(3)=1.0_ki/2.0_ki*abb56(3)
      abb56(13)=abb56(7)*abb56(13)
      abb56(18)=abb56(7)*abb56(10)
      abb56(21)=abb56(9)*abb56(18)
      abb56(13)=abb56(13)+abb56(21)
      abb56(13)=1.0_ki/2.0_ki*abb56(13)
      abb56(21)=spbk5k1*spae2k5
      abb56(22)=abb56(21)*abb56(7)*abb56(14)
      abb56(23)=-abb56(11)*abb56(9)
      abb56(17)=abb56(22)+abb56(23)-abb56(17)
      abb56(17)=1.0_ki/2.0_ki*abb56(17)
      abb56(22)=1.0_ki/2.0_ki*abb56(7)
      abb56(9)=abb56(9)*abb56(22)
      abb56(23)=spak1k5*spbk5e2
      abb56(24)=abb56(23)*abb56(16)*abb56(7)
      abb56(25)=abb56(11)*spbe2e1
      abb56(26)=-spae1e5*abb56(25)
      abb56(20)=abb56(24)+abb56(26)-abb56(20)
      abb56(20)=1.0_ki/2.0_ki*abb56(20)
      abb56(24)=1.0_ki/2.0_ki*spae1e5
      abb56(26)=abb56(24)*abb56(7)
      abb56(27)=spbe2e1*abb56(26)
      abb56(28)=spbe2k1*spae1k1
      abb56(29)=-spae1k2*spbk2e2
      abb56(29)=-abb56(28)+abb56(29)
      abb56(29)=abb56(12)*abb56(29)
      abb56(30)=spbk2k1*spae2k2
      abb56(21)=-abb56(21)+abb56(30)
      abb56(21)=abb56(21)*spae1k1*spbe5e2
      abb56(21)=abb56(21)+abb56(29)
      abb56(21)=abb56(7)*abb56(21)
      abb56(29)=spae1e2*spbe5e2*abb56(11)
      abb56(21)=abb56(29)+abb56(21)
      abb56(21)=1.0_ki/2.0_ki*abb56(21)
      abb56(29)=abb56(22)*spbe5e2
      abb56(30)=-spae1e2*abb56(29)
      abb56(4)=abb56(4)-abb56(28)
      abb56(4)=abb56(7)*abb56(4)
      abb56(28)=-spae1k2*abb56(19)
      abb56(4)=abb56(28)+abb56(4)
      abb56(4)=1.0_ki/2.0_ki*abb56(4)
      abb56(2)=spak1k2*abb56(2)
      abb56(23)=-abb56(23)*spae2e5
      abb56(2)=abb56(23)+abb56(2)
      abb56(2)=spbk1e1*abb56(2)
      abb56(23)=spak1e2*spbk1e1
      abb56(5)=-abb56(23)-abb56(5)
      abb56(5)=abb56(10)*abb56(5)
      abb56(2)=abb56(5)+abb56(2)
      abb56(2)=abb56(7)*abb56(2)
      abb56(5)=spae2e5*abb56(25)
      abb56(2)=abb56(5)+abb56(2)
      abb56(2)=1.0_ki/2.0_ki*abb56(2)
      abb56(5)=abb56(22)*spae2e5
      abb56(25)=-spbe2e1*abb56(5)
      abb56(6)=abb56(6)-abb56(23)
      abb56(6)=abb56(7)*abb56(6)
      abb56(15)=-spbk2e1*abb56(15)
      abb56(6)=abb56(15)+abb56(6)
      abb56(6)=1.0_ki/2.0_ki*abb56(6)
      abb56(8)=-abb56(11)*abb56(8)
      abb56(11)=-spbe5e1*abb56(11)*spae2e5
      abb56(15)=abb56(7)*abb56(12)
      abb56(23)=abb56(22)*spbe5e1
      abb56(28)=-spae2k5*abb56(23)
      abb56(31)=spae2k5*abb56(29)
      abb56(26)=-spbk5e2*abb56(26)
      abb56(32)=spbk5e2*abb56(5)
      abb56(23)=spae2k2*abb56(23)
      abb56(33)=abb56(10)*abb56(23)
      abb56(34)=-spae2k2*abb56(29)
      abb56(19)=abb56(19)*abb56(24)
      abb56(24)=abb56(12)*abb56(19)
      abb56(35)=-spbk2e2*abb56(5)
      abb56(14)=abb56(22)*abb56(14)
      abb56(10)=abb56(10)*abb56(14)
      abb56(29)=-spae1k1*abb56(29)
      abb56(16)=abb56(22)*abb56(16)
      abb56(12)=abb56(12)*abb56(16)
      abb56(5)=-spbk1e1*abb56(5)
      R2d56=0.0_ki
      rat2 = rat2 + R2d56
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='56' value='", &
          & R2d56, "'/>"
      end if
   end subroutine
end module p2_part21part21_part25part25part21_abbrevd56h0_qp
