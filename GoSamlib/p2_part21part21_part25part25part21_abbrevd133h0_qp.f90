module     p2_part21part21_part25part25part21_abbrevd133h0_qp
   use p2_part21part21_part25part25part21_config, only: ki => ki_qp
   use p2_part21part21_part25part25part21_kinematics_qp, only: epstensor
   use p2_part21part21_part25part25part21_globalsh0_qp
   implicit none
   private
   complex(ki), dimension(55), public :: abb133
   complex(ki), public :: R2d133
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
      abb133(1)=1.0_ki/(es34-es51-es12)
      abb133(2)=sqrt(mdlMT**2)
      abb133(3)=spak2l3**(-1)
      abb133(4)=spbl3k2**(-1)
      abb133(5)=spak2l4**(-1)
      abb133(6)=spbl4k2**(-1)
      abb133(7)=c2-c1
      abb133(8)=mdlGC45*mdlGC7
      abb133(8)=abb133(8)**2*abb133(1)*mdlGC6*spbe5e2*spae2e5
      abb133(9)=abb133(8)*i_
      abb133(10)=abb133(7)*abb133(9)
      abb133(11)=spbk5e1*spae1k5
      abb133(12)=abb133(11)*abb133(10)
      abb133(13)=-5.0_ki/6.0_ki*abb133(12)
      abb133(14)=spae1k1*spbk2k1
      abb133(15)=spbk5e1*spak2k5
      abb133(16)=abb133(14)*abb133(15)
      abb133(17)=spae1k5*spbk5k2
      abb133(18)=spbk1e1*spak1k2
      abb133(19)=abb133(17)*abb133(18)
      abb133(16)=abb133(16)+abb133(19)
      abb133(19)=spbk5l3*spae1l3
      abb133(20)=spak2k5*spbk2e1
      abb133(21)=abb133(19)*abb133(20)
      abb133(22)=spbk5k2*spae1k2
      abb133(23)=spal3k5*spbl3e1
      abb133(24)=abb133(22)*abb133(23)
      abb133(21)=-abb133(16)+abb133(21)+abb133(24)
      abb133(24)=abb133(5)*abb133(6)
      abb133(25)=1.0_ki/2.0_ki*abb133(24)
      abb133(21)=abb133(21)*abb133(25)
      abb133(25)=spbk5l4*spae1l4
      abb133(26)=abb133(25)*abb133(20)
      abb133(27)=spal4k5*spbl4e1
      abb133(28)=abb133(27)*abb133(22)
      abb133(26)=abb133(26)+abb133(28)
      abb133(28)=spae1k1*spbk5k1
      abb133(29)=abb133(28)*abb133(20)
      abb133(30)=spbk1e1*spae1k2
      abb133(31)=spak1k5*spbk5k2*abb133(30)
      abb133(26)=-abb133(29)-abb133(31)+1.0_ki/2.0_ki*abb133(26)
      abb133(29)=abb133(3)*abb133(4)
      abb133(26)=abb133(26)*abb133(29)
      abb133(21)=abb133(26)+abb133(21)
      abb133(26)=mdlMh**2
      abb133(21)=abb133(21)*abb133(26)
      abb133(31)=abb133(25)*abb133(23)
      abb133(32)=abb133(27)*abb133(19)
      abb133(33)=spak2l3*spbl4k2
      abb133(34)=spae1l4*abb133(33)*spbl3e1
      abb133(35)=spbl3k2*spak2l4
      abb133(36)=spbl4e1*abb133(35)*spae1l3
      abb133(31)=abb133(36)+abb133(34)+abb133(31)+abb133(32)
      abb133(32)=spae1k2*spbl4k2
      abb133(34)=spae1k5*spbk5l4
      abb133(32)=abb133(32)-abb133(34)
      abb133(34)=-spak1l4*abb133(32)
      abb133(36)=abb133(19)*spak1k5
      abb133(37)=spbl3k2*spae1l3
      abb133(38)=abb133(37)*spak1k2
      abb133(36)=abb133(36)-abb133(38)+1.0_ki/2.0_ki*abb133(34)
      abb133(36)=abb133(36)*spbk1e1
      abb133(38)=spbk2e1*spak2l4
      abb133(39)=spbk5e1*spal4k5
      abb133(38)=abb133(38)-abb133(39)
      abb133(39)=-spbl4k1*abb133(38)
      abb133(40)=abb133(23)*spbk5k1
      abb133(41)=spak2l3*spbl3e1
      abb133(42)=abb133(41)*spbk2k1
      abb133(40)=abb133(40)-abb133(42)+1.0_ki/2.0_ki*abb133(39)
      abb133(40)=abb133(40)*spae1k1
      abb133(42)=abb133(2)**2
      abb133(11)=abb133(11)*abb133(42)
      abb133(21)=abb133(21)-abb133(36)+1.0_ki/2.0_ki*abb133(31)-abb133(40)-abb1&
      &33(11)
      abb133(31)=-abb133(42)*abb133(7)
      abb133(9)=abb133(31)*abb133(9)
      abb133(21)=abb133(21)*abb133(9)
      abb133(36)=abb133(24)*abb133(26)
      abb133(16)=abb133(16)*abb133(36)
      abb133(39)=abb133(39)*spae1k1
      abb133(34)=abb133(34)*spbk1e1
      abb133(16)=abb133(16)+abb133(39)+abb133(34)
      abb133(11)=2.0_ki*abb133(11)+1.0_ki/2.0_ki*abb133(16)
      abb133(11)=abb133(11)*abb133(10)
      abb133(16)=abb133(20)*abb133(24)
      abb133(20)=abb133(29)*spak2k5
      abb133(34)=abb133(20)*spbk2e1
      abb133(16)=abb133(34)-abb133(16)
      abb133(16)=abb133(16)*abb133(26)
      abb133(34)=spbk1e1*spak1k5
      abb133(16)=-abb133(27)+abb133(23)+abb133(16)+2.0_ki*abb133(34)
      abb133(16)=-abb133(16)*abb133(9)
      abb133(22)=abb133(22)*abb133(24)
      abb133(23)=abb133(29)*spbk5k2
      abb133(27)=abb133(23)*spae1k2
      abb133(22)=abb133(27)-abb133(22)
      abb133(22)=abb133(22)*abb133(26)
      abb133(19)=abb133(22)+2.0_ki*abb133(28)-abb133(25)+abb133(19)
      abb133(19)=-abb133(19)*abb133(9)
      abb133(9)=-4.0_ki*abb133(9)
      abb133(22)=1.0_ki/2.0_ki*i_
      abb133(8)=abb133(8)*abb133(22)
      abb133(22)=abb133(31)*abb133(8)
      abb133(25)=abb133(38)*abb133(22)
      abb133(27)=abb133(7)*abb133(8)
      abb133(28)=-abb133(38)*abb133(27)
      abb133(31)=abb133(32)*abb133(22)
      abb133(32)=-abb133(32)*abb133(27)
      abb133(34)=spal3k5*abb133(36)*spbk5k2
      abb133(38)=abb133(34)*abb133(18)
      abb133(39)=spal3k5*spbk5l4
      abb133(33)=abb133(39)+abb133(33)
      abb133(39)=spbk1e1*spak1l4
      abb133(40)=abb133(33)*abb133(39)
      abb133(43)=spbk2e1*spak2l3
      abb133(44)=spbk5e1*spal3k5
      abb133(43)=abb133(43)-abb133(44)
      abb133(44)=abb133(43)*abb133(42)
      abb133(38)=-abb133(44)+abb133(38)+abb133(40)
      abb133(38)=-abb133(38)*abb133(27)
      abb133(40)=abb133(43)*abb133(27)
      abb133(33)=-abb133(33)*abb133(27)
      abb133(36)=spbk5l3*abb133(36)*spak2k5
      abb133(43)=abb133(36)*abb133(14)
      abb133(44)=spbk5l3*spal4k5
      abb133(35)=abb133(44)+abb133(35)
      abb133(44)=spae1k1*spbl4k1
      abb133(45)=abb133(35)*abb133(44)
      abb133(46)=spae1k2*spbl3k2
      abb133(47)=spae1k5*spbk5l3
      abb133(46)=abb133(46)-abb133(47)
      abb133(47)=abb133(46)*abb133(42)
      abb133(43)=-abb133(47)+abb133(43)+abb133(45)
      abb133(43)=-abb133(43)*abb133(27)
      abb133(45)=abb133(46)*abb133(27)
      abb133(35)=-abb133(35)*abb133(27)
      abb133(15)=abb133(15)*abb133(24)
      abb133(46)=abb133(20)*spbk5e1
      abb133(15)=abb133(46)-abb133(15)
      abb133(46)=1.0_ki/2.0_ki*abb133(26)
      abb133(47)=abb133(15)*abb133(46)
      abb133(48)=spak2l4*spbl4e1
      abb133(41)=-abb133(48)-2.0_ki*abb133(18)+abb133(47)+abb133(41)
      abb133(41)=abb133(41)*abb133(42)
      abb133(47)=-es12+es34-es51
      abb133(47)=abb133(47)*abb133(24)
      abb133(48)=abb133(29)*abb133(26)
      abb133(49)=abb133(47)*abb133(48)
      abb133(18)=abb133(49)*abb133(18)
      abb133(20)=abb133(20)*spbk5l4*abb133(39)
      abb133(18)=abb133(20)-abb133(18)
      abb133(18)=abb133(18)*abb133(46)
      abb133(18)=abb133(41)+abb133(18)
      abb133(18)=-abb133(18)*abb133(10)
      abb133(20)=abb133(27)*abb133(26)
      abb133(15)=-abb133(15)*abb133(20)
      abb133(26)=-abb133(48)*abb133(7)
      abb133(39)=-abb133(8)*spak2k5*abb133(26)
      abb133(41)=-spbk5l4*abb133(39)
      abb133(36)=-abb133(36)*abb133(27)
      abb133(17)=abb133(17)*abb133(24)
      abb133(24)=abb133(23)*spae1k5
      abb133(17)=abb133(24)-abb133(17)
      abb133(24)=abb133(17)*abb133(46)
      abb133(48)=spbl4k2*spae1l4
      abb133(24)=-abb133(48)-2.0_ki*abb133(14)+abb133(24)+abb133(37)
      abb133(24)=abb133(24)*abb133(42)
      abb133(14)=abb133(49)*abb133(14)
      abb133(23)=abb133(23)*spal4k5*abb133(44)
      abb133(14)=abb133(14)-abb133(23)
      abb133(14)=abb133(14)*abb133(46)
      abb133(14)=abb133(24)-abb133(14)
      abb133(14)=-abb133(14)*abb133(10)
      abb133(17)=-abb133(17)*abb133(20)
      abb133(20)=-abb133(8)*spbk5k2*abb133(26)
      abb133(23)=-spal4k5*abb133(20)
      abb133(24)=-abb133(34)*abb133(27)
      abb133(26)=abb133(29)*abb133(47)*mdlMh**4
      abb133(26)=abb133(26)+4.0_ki*abb133(42)
      abb133(10)=abb133(26)*abb133(10)
      abb133(26)=spbk2e1*spak1k2
      abb133(29)=spbk5e1*spak1k5
      abb133(26)=abb133(26)+abb133(29)
      abb133(29)=abb133(26)*abb133(22)
      abb133(26)=-abb133(26)*abb133(27)
      abb133(34)=spak1k2*spbl3k2
      abb133(37)=spbk5l3*spak1k5
      abb133(34)=abb133(34)-abb133(37)
      abb133(34)=-abb133(34)*abb133(27)
      abb133(37)=spak1k5*abb133(20)
      abb133(42)=spae1k2*spbk2k1
      abb133(44)=spae1k5*spbk5k1
      abb133(42)=abb133(42)+abb133(44)
      abb133(44)=abb133(42)*abb133(22)
      abb133(42)=-abb133(42)*abb133(27)
      abb133(46)=spbk2k1*spak2l3
      abb133(47)=spal3k5*spbk5k1
      abb133(46)=abb133(46)-abb133(47)
      abb133(46)=-abb133(46)*abb133(27)
      abb133(47)=spbk5k1*abb133(39)
      abb133(48)=abb133(22)*spae1k1
      abb133(49)=spbk5e1*abb133(48)
      abb133(7)=-abb133(8)*spae1k1*abb133(7)
      abb133(8)=spbk5e1*abb133(7)
      abb133(50)=-spbk5l3*abb133(7)
      abb133(20)=spae1k1*abb133(20)
      abb133(48)=spbk2e1*abb133(48)
      abb133(51)=spbk2e1*abb133(7)
      abb133(7)=spbl3k2*abb133(7)
      abb133(52)=spbk1e1*spae1k5
      abb133(53)=abb133(52)*abb133(22)
      abb133(52)=-abb133(52)*abb133(27)
      abb133(54)=abb133(27)*spbk1e1
      abb133(55)=spal3k5*abb133(54)
      abb133(39)=spbk1e1*abb133(39)
      abb133(22)=abb133(30)*abb133(22)
      abb133(27)=-abb133(30)*abb133(27)
      abb133(30)=-spak2l3*abb133(54)
      R2d133=abb133(13)
      rat2 = rat2 + R2d133
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='133' value='", &
          & R2d133, "'/>"
      end if
   end subroutine
end module p2_part21part21_part25part25part21_abbrevd133h0_qp
