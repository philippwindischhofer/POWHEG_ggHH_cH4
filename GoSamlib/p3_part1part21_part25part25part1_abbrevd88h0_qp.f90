module     p3_part1part21_part25part25part1_abbrevd88h0_qp
   use p3_part1part21_part25part25part1_config, only: ki => ki_qp
   use p3_part1part21_part25part25part1_kinematics_qp, only: epstensor
   use p3_part1part21_part25part25part1_globalsh0_qp
   implicit none
   private
   complex(ki), dimension(41), public :: abb88
   complex(ki), public :: R2d88
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p3_part1part21_part25part25part1_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p3_part1part21_part25part25part1_kinematics_qp
      use p3_part1part21_part25part25part1_model_qp
      use p3_part1part21_part25part25part1_color_qp, only: TR
      use p3_part1part21_part25part25part1_globalsl1_qp, only: epspow
      implicit none
      abb88(1)=es51**(-1)
      abb88(2)=sqrt(mdlMT**2)
      abb88(3)=spak2l3**(-1)
      abb88(4)=spbl3k2**(-1)
      abb88(5)=spak2l4**(-1)
      abb88(6)=spbl4k2**(-1)
      abb88(7)=mdlGC7**3*abb88(1)*c1
      abb88(8)=mdlGC45**2
      abb88(9)=abb88(7)*abb88(8)*TR
      abb88(10)=-8.0_ki/3.0_ki*spae2k5*spbe2k1*abb88(9)
      abb88(11)=abb88(6)*spbk2e2
      abb88(12)=abb88(11)*spbl4k1
      abb88(13)=abb88(3)*abb88(4)
      abb88(14)=abb88(13)*spak2l4
      abb88(14)=abb88(14)+abb88(6)
      abb88(15)=-spbk2k1*abb88(14)*spbl4e2
      abb88(15)=abb88(15)+abb88(12)
      abb88(15)=spae2k5*abb88(15)
      abb88(16)=spal3k5*spbl3k2
      abb88(17)=spae2k2*spbe2k1
      abb88(18)=abb88(17)*abb88(16)
      abb88(19)=spae2k5*spbl3k1
      abb88(20)=spbk2e2*spak2l3*abb88(19)
      abb88(18)=abb88(18)+abb88(20)
      abb88(18)=abb88(6)*abb88(18)
      abb88(20)=spal4k5*abb88(17)
      abb88(18)=abb88(20)+abb88(18)
      abb88(18)=abb88(5)*abb88(18)
      abb88(20)=abb88(13)*spbl4k2
      abb88(20)=abb88(20)+abb88(5)
      abb88(21)=spbe2k1*spae2l4
      abb88(22)=-spak2k5*abb88(21)*abb88(20)
      abb88(15)=abb88(22)+abb88(18)+abb88(15)
      abb88(18)=mdlMh**2
      abb88(15)=abb88(15)*abb88(18)
      abb88(22)=spal3k5*spbl4l3
      abb88(23)=-abb88(21)*abb88(22)
      abb88(24)=spbl4e2*spal4k5
      abb88(25)=spal3k5*spbl3e2
      abb88(26)=abb88(25)+2.0_ki*abb88(24)
      abb88(27)=spbk2k1*spae2k2
      abb88(28)=-abb88(26)*abb88(27)
      abb88(29)=spae2l4*spbl4k1
      abb88(30)=spbl3k1*spae2l3
      abb88(31)=abb88(30)+2.0_ki*abb88(29)
      abb88(32)=spak2k5*spbk2e2
      abb88(33)=-abb88(31)*abb88(32)
      abb88(34)=abb88(2)**2
      abb88(35)=abb88(34)*spae2k5
      abb88(36)=abb88(35)*spbe2k1
      abb88(37)=-spbl4e2*spal3l4*abb88(19)
      abb88(15)=abb88(15)-2.0_ki*abb88(36)+abb88(33)+abb88(28)+abb88(23)+abb88(&
      &37)
      abb88(23)=2.0_ki*TR
      abb88(7)=abb88(7)*abb88(23)*abb88(8)
      abb88(8)=abb88(7)*abb88(34)
      abb88(15)=abb88(15)*abb88(8)
      abb88(23)=abb88(24)+abb88(25)
      abb88(24)=-abb88(23)*abb88(27)
      abb88(25)=abb88(29)+abb88(30)
      abb88(28)=-abb88(25)*abb88(32)
      abb88(24)=-6.0_ki*abb88(36)+abb88(24)+abb88(28)
      abb88(24)=abb88(24)*abb88(7)
      abb88(28)=abb88(5)*abb88(6)
      abb88(29)=abb88(28)*spae2k2
      abb88(30)=abb88(13)*spae2k2
      abb88(33)=abb88(30)+2.0_ki*abb88(29)
      abb88(36)=abb88(18)*spbk2k1
      abb88(33)=abb88(33)*abb88(36)
      abb88(31)=abb88(33)+3.0_ki*abb88(27)+abb88(31)
      abb88(31)=abb88(31)*abb88(8)
      abb88(29)=abb88(30)+abb88(29)
      abb88(29)=abb88(29)*abb88(36)
      abb88(25)=abb88(29)+abb88(27)+abb88(25)
      abb88(25)=abb88(25)*abb88(7)
      abb88(29)=spbl3k1*spal3l4
      abb88(33)=abb88(32)*abb88(29)
      abb88(37)=abb88(13)*spbk2e2
      abb88(38)=spak2l4*abb88(37)
      abb88(11)=abb88(38)+abb88(11)
      abb88(38)=abb88(18)*spak2k5
      abb88(39)=abb88(38)*spbk2k1
      abb88(11)=abb88(11)*abb88(39)
      abb88(40)=spal4k5*abb88(34)*spbe2k1
      abb88(11)=abb88(11)+abb88(33)-3.0_ki*abb88(40)
      abb88(11)=abb88(11)*abb88(7)
      abb88(14)=-abb88(14)*abb88(36)
      abb88(14)=-abb88(29)+abb88(14)
      abb88(14)=abb88(14)*abb88(7)
      abb88(29)=abb88(27)*abb88(22)
      abb88(33)=abb88(5)*spae2k2
      abb88(30)=spbl4k2*abb88(30)
      abb88(30)=abb88(30)+abb88(33)
      abb88(30)=abb88(30)*abb88(39)
      abb88(33)=spbl4k1*abb88(35)
      abb88(29)=abb88(30)+abb88(29)-3.0_ki*abb88(33)
      abb88(29)=abb88(29)*abb88(7)
      abb88(30)=abb88(34)*abb88(9)
      abb88(33)=4.0_ki*abb88(30)
      abb88(36)=-spbe2k1*spal3k5*abb88(33)
      abb88(19)=-abb88(19)*abb88(33)
      abb88(18)=abb88(7)*abb88(18)
      abb88(33)=spbl3k1*spak2l3
      abb88(39)=abb88(28)*spbk2e2
      abb88(40)=-abb88(39)*abb88(33)
      abb88(13)=2.0_ki*abb88(13)+3.0_ki*abb88(28)
      abb88(13)=abb88(13)*abb88(34)
      abb88(41)=-spbe2k1*abb88(13)
      abb88(12)=abb88(41)-abb88(12)+abb88(40)
      abb88(12)=abb88(18)*spak2k5*abb88(12)
      abb88(33)=abb88(28)*abb88(33)
      abb88(40)=abb88(6)*spbl4k1
      abb88(33)=abb88(40)+abb88(33)
      abb88(33)=abb88(33)*abb88(18)
      abb88(13)=-spae2k5*abb88(13)
      abb88(16)=abb88(16)*abb88(6)
      abb88(16)=abb88(16)+spal4k5
      abb88(16)=abb88(16)*abb88(5)
      abb88(40)=-spae2k2*abb88(16)
      abb88(13)=abb88(40)+abb88(13)
      abb88(13)=abb88(18)*spbk2k1*abb88(13)
      abb88(40)=-abb88(37)-2.0_ki*abb88(39)
      abb88(40)=abb88(40)*abb88(38)
      abb88(26)=abb88(40)-3.0_ki*abb88(32)-abb88(26)
      abb88(26)=abb88(26)*abb88(8)
      abb88(37)=-abb88(37)-abb88(39)
      abb88(37)=abb88(37)*abb88(38)
      abb88(23)=abb88(37)-abb88(32)-abb88(23)
      abb88(23)=abb88(23)*abb88(7)
      abb88(30)=12.0_ki*abb88(30)
      abb88(9)=4.0_ki*abb88(9)
      abb88(20)=abb88(20)*abb88(38)
      abb88(20)=abb88(22)+abb88(20)
      abb88(20)=abb88(20)*abb88(7)
      abb88(16)=-abb88(18)*abb88(16)
      abb88(8)=-abb88(21)*abb88(8)
      abb88(21)=abb88(7)*abb88(27)
      abb88(18)=abb88(28)*abb88(18)
      abb88(17)=-abb88(34)*abb88(18)*abb88(17)
      abb88(22)=abb88(27)*abb88(18)
      abb88(27)=spbl4e2*abb88(7)*abb88(35)
      abb88(28)=-abb88(7)*abb88(32)
      abb88(34)=spbk2e2*abb88(35)*abb88(18)
      abb88(32)=-abb88(32)*abb88(18)
      R2d88=abb88(10)
      rat2 = rat2 + R2d88
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='88' value='", &
          & R2d88, "'/>"
      end if
   end subroutine
end module p3_part1part21_part25part25part1_abbrevd88h0_qp
