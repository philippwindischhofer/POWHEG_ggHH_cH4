module     p2_part21part21_part25part25part21_abbrevd117h0
   use p2_part21part21_part25part25part21_config, only: ki
   use p2_part21part21_part25part25part21_kinematics, only: epstensor
   use p2_part21part21_part25part25part21_globalsh0
   implicit none
   private
   complex(ki), dimension(84), public :: abb117
   complex(ki), public :: R2d117
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p2_part21part21_part25part25part21_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p2_part21part21_part25part25part21_kinematics
      use p2_part21part21_part25part25part21_model
      use p2_part21part21_part25part25part21_color, only: TR
      use p2_part21part21_part25part25part21_globalsl1, only: epspow
      implicit none
      abb117(1)=es12**(-1)
      abb117(2)=sqrt(mdlMT**2)
      abb117(3)=spak2l3**(-1)
      abb117(4)=spbl3k2**(-1)
      abb117(5)=spak2l4**(-1)
      abb117(6)=spbl4k2**(-1)
      abb117(7)=c1-c2
      abb117(8)=mdlGC7*mdlGC45
      abb117(8)=i_*mdlGC6*abb117(8)**2
      abb117(9)=abb117(1)*abb117(8)
      abb117(10)=-abb117(9)*abb117(7)
      abb117(11)=spak1e2*spbe2k1
      abb117(12)=abb117(10)*abb117(11)
      abb117(13)=abb117(12)*spae1e5
      abb117(14)=spbe5e1*abb117(13)
      abb117(15)=abb117(7)*spae1e2*spbe2e1
      abb117(16)=abb117(9)*abb117(15)
      abb117(17)=1.0_ki/2.0_ki*abb117(16)
      abb117(18)=abb117(17)*spak1e5
      abb117(19)=spbe5k1*abb117(18)
      abb117(14)=abb117(14)+abb117(19)
      abb117(14)=5.0_ki/3.0_ki*abb117(14)
      abb117(9)=abb117(9)*abb117(2)**2
      abb117(19)=abb117(9)*abb117(15)
      abb117(20)=abb117(19)*spak2e5
      abb117(21)=spal4k5*spbk5e5
      abb117(22)=abb117(21)*abb117(20)
      abb117(23)=abb117(19)*spbe5l3
      abb117(24)=abb117(23)*spak2l3
      abb117(25)=-spal4e5*abb117(24)
      abb117(22)=abb117(22)+abb117(25)
      abb117(22)=spbl4k2*abb117(22)
      abb117(25)=abb117(19)*spbe5k2
      abb117(26)=spbk5l4*spae5k5
      abb117(27)=abb117(26)*abb117(25)
      abb117(28)=abb117(19)*spal3e5
      abb117(29)=abb117(28)*spbl3k2
      abb117(30)=-spbe5l4*abb117(29)
      abb117(27)=abb117(27)+abb117(30)
      abb117(27)=spak2l4*abb117(27)
      abb117(22)=abb117(22)+abb117(27)
      abb117(27)=abb117(6)*abb117(5)
      abb117(30)=mdlMh**2
      abb117(31)=abb117(27)*abb117(30)
      abb117(32)=abb117(11)*abb117(31)
      abb117(7)=-abb117(9)*abb117(7)
      abb117(9)=abb117(32)*abb117(7)
      abb117(33)=spae1k2*spbe5e1
      abb117(34)=abb117(33)*spae5k5
      abb117(35)=-abb117(9)*abb117(34)
      abb117(36)=abb117(19)*abb117(31)
      abb117(37)=spak1k2*spbe5k1
      abb117(38)=1.0_ki/2.0_ki*spae5k5
      abb117(39)=-abb117(38)*abb117(36)*abb117(37)
      abb117(40)=spae5k5*abb117(24)
      abb117(35)=abb117(40)+abb117(35)+abb117(39)
      abb117(35)=spbk5k2*abb117(35)
      abb117(39)=spbk2e1*spae1e5
      abb117(40)=abb117(39)*spbk5e5
      abb117(41)=-abb117(9)*abb117(40)
      abb117(42)=1.0_ki/2.0_ki*spbk5e5
      abb117(43)=-spbk2k1*spak1e5*abb117(42)*abb117(36)
      abb117(44)=spbk5e5*abb117(29)
      abb117(41)=abb117(44)+abb117(41)+abb117(43)
      abb117(41)=spak2k5*abb117(41)
      abb117(43)=spbk2e1*spak2e5
      abb117(44)=spae1l3*abb117(9)*abb117(43)
      abb117(45)=abb117(31)*spbk2k1
      abb117(46)=1.0_ki/2.0_ki*abb117(19)
      abb117(47)=abb117(46)*spak2e5
      abb117(48)=spak1l3*abb117(47)*abb117(45)
      abb117(44)=abb117(44)+abb117(48)
      abb117(44)=spbe5l3*abb117(44)
      abb117(48)=spae1k2*spbe5k2
      abb117(49)=spbl3e1*abb117(9)*abb117(48)
      abb117(50)=abb117(31)*spak1k2
      abb117(51)=abb117(46)*spbe5k2
      abb117(52)=spbl3k1*abb117(51)*abb117(50)
      abb117(49)=abb117(49)+abb117(52)
      abb117(49)=spal3e5*abb117(49)
      abb117(52)=abb117(7)*abb117(11)
      abb117(53)=abb117(4)*abb117(3)
      abb117(30)=abb117(53)*abb117(30)
      abb117(54)=abb117(30)*abb117(52)
      abb117(55)=abb117(43)*abb117(54)
      abb117(56)=spal3e5*spbl3e1
      abb117(57)=abb117(56)*abb117(52)
      abb117(55)=abb117(55)+abb117(57)
      abb117(57)=spae1l4*abb117(55)
      abb117(58)=abb117(19)*abb117(30)
      abb117(59)=spbk2k1*spak2e5
      abb117(60)=abb117(59)*abb117(58)
      abb117(28)=abb117(28)*spbl3k1
      abb117(60)=abb117(60)+abb117(28)
      abb117(61)=spak1l4*abb117(60)
      abb117(57)=abb117(57)+1.0_ki/2.0_ki*abb117(61)
      abb117(57)=spbe5l4*abb117(57)
      abb117(61)=abb117(48)*abb117(54)
      abb117(62)=spbe5l3*spae1l3
      abb117(63)=abb117(62)*abb117(52)
      abb117(61)=abb117(61)+abb117(63)
      abb117(63)=spbl4e1*abb117(61)
      abb117(64)=spak1k2*spbe5k2
      abb117(65)=abb117(64)*abb117(58)
      abb117(23)=abb117(23)*spak1l3
      abb117(65)=abb117(65)+abb117(23)
      abb117(66)=spbl4k1*abb117(65)
      abb117(63)=abb117(63)+1.0_ki/2.0_ki*abb117(66)
      abb117(63)=spal4e5*abb117(63)
      abb117(66)=2.0_ki*spbe5e1
      abb117(11)=abb117(11)*spae1e5
      abb117(67)=-abb117(11)*abb117(66)*abb117(10)
      abb117(68)=abb117(16)*spak1e5*spbe5k1
      abb117(67)=abb117(67)-abb117(68)
      abb117(67)=abb117(67)*abb117(2)**4
      abb117(69)=2.0_ki*spbk5e1
      abb117(69)=abb117(69)*spae5k5
      abb117(61)=-abb117(61)*abb117(69)
      abb117(70)=2.0_ki*spae1k5
      abb117(70)=abb117(70)*spbk5e5
      abb117(55)=-abb117(55)*abb117(70)
      abb117(65)=-spbk5k1*spae5k5*abb117(65)
      abb117(60)=-spak1k5*spbk5e5*abb117(60)
      abb117(71)=spbe5e1*abb117(52)
      abb117(26)=abb117(26)*spae1l4
      abb117(72)=-abb117(71)*abb117(26)
      abb117(7)=abb117(11)*abb117(7)
      abb117(11)=abb117(21)*spbl4e1
      abb117(21)=-abb117(7)*abb117(11)
      abb117(73)=abb117(19)*spak1e5
      abb117(74)=abb117(73)*abb117(42)
      abb117(75)=-spbl4k1*spal4k5*abb117(74)
      abb117(76)=abb117(19)*spbe5k1
      abb117(77)=abb117(76)*abb117(38)
      abb117(78)=-spak1l4*spbk5l4*abb117(77)
      abb117(21)=abb117(60)+abb117(65)+abb117(55)+abb117(61)+abb117(63)+abb117(&
      &57)+abb117(78)+abb117(75)+abb117(21)+abb117(72)+abb117(49)+abb117(44)+ab&
      &b117(41)+abb117(35)+abb117(67)+1.0_ki/2.0_ki*abb117(22)
      abb117(10)=abb117(32)*abb117(10)
      abb117(22)=-abb117(10)*abb117(34)
      abb117(32)=abb117(17)*spae5k5
      abb117(34)=abb117(32)*spbe5k1
      abb117(35)=-abb117(50)*abb117(34)
      abb117(22)=abb117(22)+abb117(35)
      abb117(22)=spbk5k2*abb117(22)
      abb117(35)=-abb117(10)*abb117(40)
      abb117(40)=abb117(17)*spbk5e5
      abb117(41)=abb117(40)*spak1e5
      abb117(44)=-abb117(45)*abb117(41)
      abb117(35)=abb117(35)+abb117(44)
      abb117(35)=spak2k5*abb117(35)
      abb117(44)=spak1l4*spbe5k1
      abb117(49)=spak2l4*spbe5k2
      abb117(49)=abb117(49)-abb117(44)
      abb117(55)=abb117(32)*spbk5l4
      abb117(49)=abb117(55)*abb117(49)
      abb117(57)=-abb117(7)*abb117(66)
      abb117(60)=-spbe5k1*abb117(73)
      abb117(57)=abb117(57)+abb117(60)
      abb117(60)=abb117(12)*spbe5e1
      abb117(61)=-abb117(60)*abb117(26)
      abb117(63)=-abb117(13)*abb117(11)
      abb117(65)=abb117(40)*spal4k5
      abb117(67)=abb117(65)*spbl4k1
      abb117(72)=-spak1e5*abb117(67)
      abb117(65)=abb117(65)*spbl4k2
      abb117(75)=spak2e5*abb117(65)
      abb117(22)=abb117(75)+abb117(72)+abb117(63)+abb117(61)+abb117(35)+2.0_ki*&
      &abb117(57)+abb117(22)+abb117(49)
      abb117(35)=-abb117(13)*abb117(66)
      abb117(35)=abb117(35)-abb117(68)
      abb117(9)=abb117(54)-abb117(9)
      abb117(48)=-abb117(9)*abb117(48)
      abb117(49)=spbe5l4*spae1l4
      abb117(49)=-abb117(70)+abb117(49)-abb117(62)
      abb117(49)=abb117(52)*abb117(49)
      abb117(48)=abb117(48)+abb117(49)
      abb117(48)=2.0_ki*abb117(48)
      abb117(43)=-abb117(9)*abb117(43)
      abb117(49)=spal4e5*spbl4e1
      abb117(49)=-abb117(69)+abb117(49)-abb117(56)
      abb117(49)=abb117(52)*abb117(49)
      abb117(43)=abb117(43)+abb117(49)
      abb117(43)=2.0_ki*abb117(43)
      abb117(49)=-8.0_ki*abb117(52)
      abb117(52)=spak2k5*abb117(51)
      abb117(54)=-spae1k5*abb117(71)
      abb117(56)=spak1k5*spbe5k1
      abb117(57)=-abb117(46)*abb117(56)
      abb117(52)=abb117(57)+abb117(52)+abb117(54)
      abb117(54)=abb117(17)*spbe5k2
      abb117(57)=spak2k5*abb117(54)
      abb117(61)=-spae1k5*abb117(60)
      abb117(56)=-abb117(17)*abb117(56)
      abb117(56)=abb117(56)+abb117(57)+abb117(61)
      abb117(57)=spbk5k2*abb117(47)
      abb117(61)=-spbk5e1*abb117(7)
      abb117(62)=1.0_ki/2.0_ki*spak1e5
      abb117(63)=abb117(62)*abb117(19)
      abb117(66)=-spbk5k1*abb117(63)
      abb117(57)=abb117(66)+abb117(57)+abb117(61)
      abb117(61)=abb117(17)*spak2e5
      abb117(66)=spbk5k2*abb117(61)
      abb117(68)=-spbk5e1*abb117(13)
      abb117(69)=-spbk5k1*abb117(18)
      abb117(66)=abb117(69)+abb117(66)+abb117(68)
      abb117(68)=spae5k5*abb117(71)
      abb117(69)=spae5k5*abb117(60)
      abb117(70)=spbk5e5*abb117(7)
      abb117(72)=spbk5e5*abb117(13)
      abb117(75)=-spae1l4*abb117(71)
      abb117(46)=-abb117(46)*abb117(44)
      abb117(78)=spak2l4*abb117(51)
      abb117(46)=abb117(78)+abb117(75)+abb117(46)
      abb117(75)=-spae1l4*abb117(60)
      abb117(44)=-abb117(17)*abb117(44)
      abb117(54)=spak2l4*abb117(54)
      abb117(44)=abb117(54)+abb117(75)+abb117(44)
      abb117(54)=-spbl4e1*abb117(7)
      abb117(63)=-spbl4k1*abb117(63)
      abb117(75)=spbl4k2*abb117(47)
      abb117(54)=abb117(75)+abb117(54)+abb117(63)
      abb117(63)=-spbl4e1*abb117(13)
      abb117(18)=-spbl4k1*abb117(18)
      abb117(61)=spbl4k2*abb117(61)
      abb117(18)=abb117(61)+abb117(63)+abb117(18)
      abb117(61)=abb117(10)*spbk2e1
      abb117(63)=spak2k5*spbk5e5
      abb117(75)=abb117(61)*abb117(63)
      abb117(71)=abb117(71)+abb117(75)
      abb117(71)=spae1l3*abb117(71)
      abb117(31)=abb117(31)*abb117(16)
      abb117(63)=spbk2k1*abb117(31)*abb117(63)
      abb117(63)=abb117(76)+abb117(63)
      abb117(63)=abb117(67)+1.0_ki/2.0_ki*abb117(63)
      abb117(63)=spak1l3*abb117(63)
      abb117(51)=-abb117(65)-abb117(51)
      abb117(51)=spak2l3*abb117(51)
      abb117(65)=abb117(12)*spae1l3
      abb117(75)=abb117(65)*abb117(11)
      abb117(51)=abb117(75)+abb117(51)+abb117(63)+abb117(71)
      abb117(60)=spae1l3*abb117(60)
      abb117(63)=abb117(17)*spak1l3
      abb117(71)=spbe5k1*abb117(63)
      abb117(75)=abb117(17)*spak2l3
      abb117(76)=-spbe5k2*abb117(75)
      abb117(60)=abb117(76)+abb117(60)+abb117(71)
      abb117(71)=-spbk5k2*abb117(75)
      abb117(76)=spbk5e1*abb117(65)
      abb117(78)=spbk5k1*abb117(63)
      abb117(71)=abb117(78)+abb117(71)+abb117(76)
      abb117(76)=-spbk5e5*abb117(65)
      abb117(65)=spbl4e1*abb117(65)
      abb117(78)=spbl4k1*abb117(63)
      abb117(75)=-spbl4k2*abb117(75)
      abb117(65)=abb117(75)+abb117(65)+abb117(78)
      abb117(75)=abb117(10)*spae1k2
      abb117(78)=spbk5k2*spae5k5
      abb117(79)=abb117(75)*abb117(78)
      abb117(7)=abb117(7)+abb117(79)
      abb117(7)=spbl3e1*abb117(7)
      abb117(79)=spak1l4*spbl3k1
      abb117(80)=-spak2l4*spbl3k2
      abb117(79)=abb117(80)+abb117(79)
      abb117(79)=abb117(55)*abb117(79)
      abb117(78)=spak1k2*abb117(31)*abb117(78)
      abb117(73)=abb117(73)+abb117(78)
      abb117(73)=spbl3k1*abb117(73)
      abb117(78)=abb117(12)*spbl3e1
      abb117(80)=abb117(26)*abb117(78)
      abb117(47)=-spbl3k2*abb117(47)
      abb117(7)=abb117(47)+abb117(80)+1.0_ki/2.0_ki*abb117(73)+abb117(79)+abb11&
      &7(7)
      abb117(13)=spbl3e1*abb117(13)
      abb117(47)=abb117(17)*spbl3k1
      abb117(73)=spak1e5*abb117(47)
      abb117(79)=abb117(17)*spbl3k2
      abb117(80)=-spak2e5*abb117(79)
      abb117(13)=abb117(80)+abb117(13)+abb117(73)
      abb117(73)=-spak2k5*abb117(79)
      abb117(80)=spae1k5*abb117(78)
      abb117(81)=spak1k5*abb117(47)
      abb117(73)=abb117(81)+abb117(73)+abb117(80)
      abb117(80)=-spae5k5*abb117(78)
      abb117(78)=spae1l4*abb117(78)
      abb117(81)=spak1l4*abb117(47)
      abb117(79)=-spak2l4*abb117(79)
      abb117(78)=abb117(79)+abb117(78)+abb117(81)
      abb117(79)=abb117(9)*abb117(33)
      abb117(36)=abb117(58)-abb117(36)
      abb117(37)=1.0_ki/2.0_ki*abb117(37)
      abb117(58)=abb117(36)*abb117(37)
      abb117(8)=abb117(53)*abb117(27)*mdlMh**4*abb117(15)*abb117(8)
      abb117(15)=abb117(42)*abb117(8)
      abb117(27)=2.0_ki*abb117(19)
      abb117(53)=abb117(27)*spbk5e5
      abb117(15)=abb117(53)+abb117(15)
      abb117(15)=spak2k5*abb117(15)
      abb117(12)=abb117(12)*abb117(30)
      abb117(81)=abb117(12)*spae1k2
      abb117(11)=abb117(11)*abb117(81)
      abb117(82)=abb117(30)*spak1k2
      abb117(67)=abb117(67)*abb117(82)
      abb117(83)=abb117(19)*spbe5l4
      abb117(84)=-spak2l4*abb117(83)
      abb117(11)=abb117(84)+abb117(24)+abb117(67)+abb117(11)+abb117(15)+abb117(&
      &79)+abb117(58)
      abb117(10)=abb117(10)-abb117(12)
      abb117(15)=-abb117(10)*abb117(33)
      abb117(16)=abb117(30)*abb117(16)
      abb117(16)=abb117(16)-abb117(31)
      abb117(24)=abb117(16)*abb117(37)
      abb117(15)=abb117(15)+abb117(24)
      abb117(24)=spbk5e1*abb117(81)
      abb117(17)=abb117(17)*abb117(30)
      abb117(31)=abb117(17)*spak1k2
      abb117(33)=spbk5k1*abb117(31)
      abb117(24)=abb117(24)+abb117(33)
      abb117(33)=-spbk5e5*abb117(81)
      abb117(37)=spbl4e1*abb117(81)
      abb117(31)=spbl4k1*abb117(31)
      abb117(31)=abb117(37)+abb117(31)
      abb117(37)=spbl3e1*abb117(75)
      abb117(47)=abb117(47)*abb117(50)
      abb117(37)=abb117(37)+abb117(47)
      abb117(9)=abb117(9)*abb117(39)
      abb117(47)=abb117(62)*spbk2k1
      abb117(50)=abb117(36)*abb117(47)
      abb117(58)=abb117(38)*abb117(8)
      abb117(27)=abb117(27)*spae5k5
      abb117(58)=abb117(27)+abb117(58)
      abb117(58)=spbk5k2*abb117(58)
      abb117(12)=abb117(12)*spbk2e1
      abb117(26)=abb117(26)*abb117(12)
      abb117(30)=abb117(30)*spbk2k1
      abb117(55)=spak1l4*abb117(55)*abb117(30)
      abb117(62)=abb117(19)*spal4e5
      abb117(67)=-spbl4k2*abb117(62)
      abb117(9)=abb117(67)+abb117(29)+abb117(55)+abb117(26)+abb117(58)+abb117(9&
      &)+abb117(50)
      abb117(10)=-abb117(10)*abb117(39)
      abb117(16)=abb117(16)*abb117(47)
      abb117(10)=abb117(10)+abb117(16)
      abb117(16)=spae1k5*abb117(12)
      abb117(17)=abb117(17)*spbk2k1
      abb117(26)=spak1k5*abb117(17)
      abb117(16)=abb117(16)+abb117(26)
      abb117(26)=-spae5k5*abb117(12)
      abb117(12)=spae1l4*abb117(12)
      abb117(17)=spak1l4*abb117(17)
      abb117(12)=abb117(12)+abb117(17)
      abb117(17)=spae1l3*abb117(61)
      abb117(29)=abb117(63)*abb117(45)
      abb117(17)=abb117(17)+abb117(29)
      abb117(19)=4.0_ki*abb117(19)
      abb117(8)=abb117(19)+abb117(8)
      abb117(29)=-abb117(36)*abb117(64)
      abb117(39)=spak1l4*abb117(83)
      abb117(45)=-spak1k5*abb117(53)
      abb117(23)=abb117(45)+abb117(39)+abb117(29)-abb117(23)
      abb117(29)=-abb117(36)*abb117(59)
      abb117(36)=spbl4k1*abb117(62)
      abb117(27)=-spbk5k1*abb117(27)
      abb117(27)=abb117(27)+abb117(36)+abb117(29)-abb117(28)
      abb117(20)=-abb117(42)*abb117(20)
      abb117(28)=-spak2e5*abb117(40)
      abb117(29)=spak2l3*abb117(40)
      abb117(36)=-spak1l3*abb117(40)
      abb117(39)=-abb117(40)*abb117(82)
      abb117(25)=-abb117(38)*abb117(25)
      abb117(38)=-spbe5k2*abb117(32)
      abb117(40)=spbl3k2*abb117(32)
      abb117(42)=-spbl3k1*abb117(32)
      abb117(30)=-abb117(32)*abb117(30)
      R2d117=abb117(14)
      rat2 = rat2 + R2d117
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='117' value='", &
          & R2d117, "'/>"
      end if
   end subroutine
end module p2_part21part21_part25part25part21_abbrevd117h0
