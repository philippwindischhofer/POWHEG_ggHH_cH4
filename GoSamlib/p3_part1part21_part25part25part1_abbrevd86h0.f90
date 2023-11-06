module     p3_part1part21_part25part25part1_abbrevd86h0
   use p3_part1part21_part25part25part1_config, only: ki
   use p3_part1part21_part25part25part1_kinematics, only: epstensor
   use p3_part1part21_part25part25part1_globalsh0
   implicit none
   private
   complex(ki), dimension(42), public :: abb86
   complex(ki), public :: R2d86
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p3_part1part21_part25part25part1_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p3_part1part21_part25part25part1_kinematics
      use p3_part1part21_part25part25part1_model
      use p3_part1part21_part25part25part1_color, only: TR
      use p3_part1part21_part25part25part1_globalsl1, only: epspow
      implicit none
      abb86(1)=es51**(-1)
      abb86(2)=sqrt(mdlMT**2)
      abb86(3)=spak2l3**(-1)
      abb86(4)=spbl3k2**(-1)
      abb86(5)=spak2l4**(-1)
      abb86(6)=spbl4k2**(-1)
      abb86(7)=mdlGC7**3*abb86(1)*c1
      abb86(8)=mdlGC45**2
      abb86(9)=abb86(7)*abb86(8)*TR
      abb86(10)=-8.0_ki/3.0_ki*spae2k5*spbe2k1*abb86(9)
      abb86(11)=abb86(4)*spbk2e2
      abb86(12)=abb86(11)*spbl3k1
      abb86(13)=abb86(5)*abb86(6)
      abb86(14)=abb86(13)*spak2l3
      abb86(14)=abb86(14)+abb86(4)
      abb86(15)=-spbk2k1*abb86(14)*spbl3e2
      abb86(15)=abb86(15)+abb86(12)
      abb86(15)=spae2k5*abb86(15)
      abb86(16)=spal4k5*spbl4k2
      abb86(17)=spae2k2*spbe2k1
      abb86(18)=abb86(17)*abb86(16)
      abb86(19)=spae2k5*spbl4k1
      abb86(20)=spbk2e2*spak2l4*abb86(19)
      abb86(18)=abb86(18)+abb86(20)
      abb86(18)=abb86(4)*abb86(18)
      abb86(20)=spal3k5*abb86(17)
      abb86(18)=abb86(20)+abb86(18)
      abb86(18)=abb86(3)*abb86(18)
      abb86(20)=abb86(13)*spbl3k2
      abb86(20)=abb86(20)+abb86(3)
      abb86(21)=spbe2k1*spae2l3
      abb86(22)=-spak2k5*abb86(21)*abb86(20)
      abb86(15)=abb86(22)+abb86(18)+abb86(15)
      abb86(18)=mdlMh**2
      abb86(15)=abb86(15)*abb86(18)
      abb86(22)=spal4k5*spbl4l3
      abb86(23)=abb86(21)*abb86(22)
      abb86(24)=spbl3e2*spal3k5
      abb86(25)=spal4k5*spbl4e2
      abb86(26)=abb86(25)+2.0_ki*abb86(24)
      abb86(27)=spbk2k1*spae2k2
      abb86(28)=-abb86(26)*abb86(27)
      abb86(29)=spae2l3*spbl3k1
      abb86(30)=spbl4k1*spae2l4
      abb86(31)=abb86(30)+2.0_ki*abb86(29)
      abb86(32)=spak2k5*spbk2e2
      abb86(33)=-abb86(31)*abb86(32)
      abb86(34)=abb86(2)**2
      abb86(35)=abb86(34)*spae2k5
      abb86(36)=abb86(35)*spbe2k1
      abb86(37)=spbl3e2*spal3l4*abb86(19)
      abb86(15)=abb86(15)-2.0_ki*abb86(36)+abb86(33)+abb86(28)+abb86(23)+abb86(&
      &37)
      abb86(23)=2.0_ki*TR
      abb86(7)=abb86(7)*abb86(23)*abb86(8)
      abb86(8)=abb86(7)*abb86(34)
      abb86(15)=abb86(15)*abb86(8)
      abb86(23)=abb86(24)+abb86(25)
      abb86(24)=-abb86(23)*abb86(27)
      abb86(25)=abb86(29)+abb86(30)
      abb86(28)=-abb86(25)*abb86(32)
      abb86(24)=-6.0_ki*abb86(36)+abb86(24)+abb86(28)
      abb86(24)=abb86(24)*abb86(7)
      abb86(28)=abb86(3)*abb86(4)
      abb86(29)=abb86(28)*spae2k2
      abb86(30)=abb86(13)*spae2k2
      abb86(33)=abb86(30)+2.0_ki*abb86(29)
      abb86(36)=abb86(18)*spbk2k1
      abb86(33)=abb86(33)*abb86(36)
      abb86(31)=abb86(33)+3.0_ki*abb86(27)+abb86(31)
      abb86(31)=abb86(31)*abb86(8)
      abb86(29)=abb86(30)+abb86(29)
      abb86(29)=abb86(29)*abb86(36)
      abb86(25)=abb86(29)+abb86(27)+abb86(25)
      abb86(25)=abb86(25)*abb86(7)
      abb86(29)=abb86(34)*abb86(9)
      abb86(33)=4.0_ki*abb86(29)
      abb86(37)=-spbe2k1*spal4k5*abb86(33)
      abb86(19)=-abb86(19)*abb86(33)
      abb86(33)=spbl4k1*spal3l4
      abb86(38)=-abb86(32)*abb86(33)
      abb86(39)=abb86(13)*spbk2e2
      abb86(40)=spak2l3*abb86(39)
      abb86(11)=abb86(40)+abb86(11)
      abb86(40)=abb86(18)*spak2k5
      abb86(41)=abb86(40)*spbk2k1
      abb86(11)=abb86(11)*abb86(41)
      abb86(42)=spal3k5*abb86(34)*spbe2k1
      abb86(11)=abb86(11)+abb86(38)-3.0_ki*abb86(42)
      abb86(11)=abb86(11)*abb86(7)
      abb86(14)=-abb86(14)*abb86(36)
      abb86(14)=abb86(33)+abb86(14)
      abb86(14)=abb86(14)*abb86(7)
      abb86(33)=-abb86(27)*abb86(22)
      abb86(36)=abb86(3)*spae2k2
      abb86(30)=spbl3k2*abb86(30)
      abb86(30)=abb86(30)+abb86(36)
      abb86(30)=abb86(30)*abb86(41)
      abb86(36)=spbl3k1*abb86(35)
      abb86(30)=abb86(30)+abb86(33)-3.0_ki*abb86(36)
      abb86(30)=abb86(30)*abb86(7)
      abb86(18)=abb86(7)*abb86(18)
      abb86(33)=spbl4k1*spak2l4
      abb86(36)=abb86(28)*spbk2e2
      abb86(38)=-abb86(36)*abb86(33)
      abb86(13)=2.0_ki*abb86(13)+3.0_ki*abb86(28)
      abb86(13)=abb86(13)*abb86(34)
      abb86(41)=-spbe2k1*abb86(13)
      abb86(12)=abb86(41)-abb86(12)+abb86(38)
      abb86(12)=abb86(18)*spak2k5*abb86(12)
      abb86(33)=abb86(28)*abb86(33)
      abb86(38)=abb86(4)*spbl3k1
      abb86(33)=abb86(38)+abb86(33)
      abb86(33)=abb86(33)*abb86(18)
      abb86(13)=-spae2k5*abb86(13)
      abb86(16)=abb86(16)*abb86(4)
      abb86(16)=abb86(16)+spal3k5
      abb86(16)=abb86(16)*abb86(3)
      abb86(38)=-spae2k2*abb86(16)
      abb86(13)=abb86(38)+abb86(13)
      abb86(13)=abb86(18)*spbk2k1*abb86(13)
      abb86(38)=-abb86(39)-2.0_ki*abb86(36)
      abb86(38)=abb86(38)*abb86(40)
      abb86(26)=abb86(38)-3.0_ki*abb86(32)-abb86(26)
      abb86(26)=abb86(26)*abb86(8)
      abb86(36)=-abb86(39)-abb86(36)
      abb86(36)=abb86(36)*abb86(40)
      abb86(23)=abb86(36)-abb86(32)-abb86(23)
      abb86(23)=abb86(23)*abb86(7)
      abb86(29)=12.0_ki*abb86(29)
      abb86(9)=4.0_ki*abb86(9)
      abb86(20)=abb86(20)*abb86(40)
      abb86(20)=-abb86(22)+abb86(20)
      abb86(20)=abb86(20)*abb86(7)
      abb86(16)=-abb86(18)*abb86(16)
      abb86(8)=-abb86(21)*abb86(8)
      abb86(21)=abb86(7)*abb86(27)
      abb86(18)=abb86(28)*abb86(18)
      abb86(17)=-abb86(34)*abb86(18)*abb86(17)
      abb86(22)=abb86(27)*abb86(18)
      abb86(27)=spbl3e2*abb86(7)*abb86(35)
      abb86(28)=-abb86(7)*abb86(32)
      abb86(34)=spbk2e2*abb86(35)*abb86(18)
      abb86(32)=-abb86(32)*abb86(18)
      R2d86=abb86(10)
      rat2 = rat2 + R2d86
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='86' value='", &
          & R2d86, "'/>"
      end if
   end subroutine
end module p3_part1part21_part25part25part1_abbrevd86h0
