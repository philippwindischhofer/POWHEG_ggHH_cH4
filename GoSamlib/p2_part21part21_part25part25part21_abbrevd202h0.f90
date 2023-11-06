module     p2_part21part21_part25part25part21_abbrevd202h0
   use p2_part21part21_part25part25part21_config, only: ki
   use p2_part21part21_part25part25part21_kinematics, only: epstensor
   use p2_part21part21_part25part25part21_globalsh0
   implicit none
   private
   complex(ki), dimension(102), public :: abb202
   complex(ki), public :: R2d202
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
      abb202(1)=sqrt(mdlMT**2)
      abb202(2)=spak2l3**(-1)
      abb202(3)=spbl3k2**(-1)
      abb202(4)=spak2l4**(-1)
      abb202(5)=spbl4k2**(-1)
      abb202(6)=c2-c1
      abb202(6)=abb202(6)*mdlGC7**3
      abb202(7)=abb202(1)*mdlGC45
      abb202(7)=abb202(6)*abb202(7)**2
      abb202(8)=-mdlMh**2*abb202(7)
      abb202(9)=abb202(5)*abb202(4)
      abb202(10)=1.0_ki/2.0_ki*abb202(9)
      abb202(11)=abb202(8)*abb202(10)
      abb202(12)=abb202(3)*abb202(2)
      abb202(13)=abb202(8)*abb202(12)
      abb202(14)=abb202(11)+abb202(13)
      abb202(15)=spbk1e1*spak1e5
      abb202(16)=abb202(14)*abb202(15)
      abb202(17)=spae2k5*spbk5e5
      abb202(18)=abb202(17)*spbk2e2
      abb202(19)=abb202(16)*abb202(18)
      abb202(20)=spak1e2*spbk2e2
      abb202(21)=abb202(15)*spbe5k1
      abb202(22)=-abb202(21)*abb202(14)*abb202(20)
      abb202(23)=abb202(12)+abb202(10)
      abb202(24)=mdlMh*mdlGC45
      abb202(24)=-abb202(6)*abb202(24)**2
      abb202(25)=abb202(1)**4
      abb202(23)=-abb202(23)*abb202(25)*abb202(24)
      abb202(26)=spae2e5*spbk2e2
      abb202(27)=abb202(26)*spbe5e1
      abb202(28)=-abb202(23)*abb202(27)
      abb202(19)=abb202(28)+abb202(19)+abb202(22)
      abb202(19)=spae1k2*abb202(19)
      abb202(22)=spbe2k1*spak1e5
      abb202(28)=spbk5e2*spae5k5
      abb202(29)=abb202(22)-abb202(28)
      abb202(30)=abb202(8)*spbe5e1
      abb202(31)=abb202(29)*abb202(30)
      abb202(32)=abb202(15)*abb202(8)*spbe5e2
      abb202(31)=abb202(32)+abb202(31)
      abb202(32)=spae2k2*abb202(2)
      abb202(33)=-abb202(32)*abb202(31)
      abb202(22)=abb202(22)*spae2k2
      abb202(28)=abb202(28)*spae2k2
      abb202(22)=abb202(22)-abb202(28)
      abb202(34)=-spbe5e1*abb202(22)
      abb202(35)=spbe5e2*spae2k2
      abb202(36)=-abb202(15)*abb202(35)
      abb202(34)=abb202(36)+abb202(34)
      abb202(36)=abb202(8)*abb202(9)
      abb202(34)=spbl3k2*abb202(36)*abb202(34)
      abb202(33)=abb202(34)+abb202(33)
      abb202(34)=1.0_ki/2.0_ki*spae1l3
      abb202(33)=abb202(34)*abb202(33)
      abb202(37)=spak1e2*spbe5k1
      abb202(17)=abb202(37)-abb202(17)
      abb202(38)=abb202(8)*spae1e5
      abb202(39)=abb202(17)*abb202(38)
      abb202(40)=spae1k1*spbe5k1
      abb202(41)=abb202(40)*abb202(8)*spae2e5
      abb202(39)=abb202(41)+abb202(39)
      abb202(41)=spbk2e2*abb202(3)
      abb202(42)=-abb202(41)*abb202(39)
      abb202(37)=abb202(37)*spbk2e2
      abb202(37)=abb202(37)-abb202(18)
      abb202(43)=-spae1e5*abb202(37)
      abb202(44)=-abb202(40)*abb202(26)
      abb202(43)=abb202(44)+abb202(43)
      abb202(43)=spak2l3*abb202(36)*abb202(43)
      abb202(42)=abb202(43)+abb202(42)
      abb202(43)=1.0_ki/2.0_ki*spbl3e1
      abb202(42)=abb202(43)*abb202(42)
      abb202(44)=abb202(40)*abb202(14)
      abb202(45)=spak1e5*spbe2k1*spae2k2
      abb202(28)=abb202(45)-abb202(28)
      abb202(45)=-abb202(44)*abb202(28)
      abb202(46)=abb202(35)*spae1e5
      abb202(23)=-abb202(23)*abb202(46)
      abb202(23)=abb202(23)+abb202(45)
      abb202(23)=spbk2e1*abb202(23)
      abb202(45)=1.0_ki/2.0_ki*spbk2e1
      abb202(47)=abb202(45)*abb202(3)
      abb202(39)=abb202(39)*abb202(47)
      abb202(48)=-abb202(7)*abb202(15)
      abb202(49)=-abb202(48)*abb202(17)
      abb202(6)=-abb202(6)*mdlGC45**2
      abb202(25)=-abb202(25)*abb202(6)
      abb202(50)=abb202(25)*spbe5e1
      abb202(51)=abb202(50)*spae2e5
      abb202(49)=abb202(49)-abb202(51)
      abb202(51)=spae1l3*abb202(49)
      abb202(39)=abb202(39)+abb202(51)
      abb202(39)=spbl3e2*abb202(39)
      abb202(51)=1.0_ki/2.0_ki*spae1k2
      abb202(52)=abb202(51)*abb202(2)
      abb202(31)=abb202(31)*abb202(52)
      abb202(53)=-abb202(7)*abb202(40)
      abb202(54)=-abb202(53)*abb202(29)
      abb202(55)=abb202(25)*spae1e5
      abb202(56)=abb202(55)*spbe5e2
      abb202(54)=abb202(54)-abb202(56)
      abb202(56)=spbl3e1*abb202(54)
      abb202(31)=abb202(31)+abb202(56)
      abb202(31)=spae2l3*abb202(31)
      abb202(56)=1.0_ki/2.0_ki*spae2l4
      abb202(57)=abb202(56)*spbl4e1
      abb202(58)=abb202(54)*abb202(57)
      abb202(59)=1.0_ki/2.0_ki*spbl4e2
      abb202(60)=abb202(59)*spae1l4
      abb202(61)=abb202(49)*abb202(60)
      abb202(62)=-abb202(7)*spbe5e1
      abb202(63)=-abb202(62)*abb202(29)
      abb202(64)=abb202(48)*spbe5e2
      abb202(63)=abb202(63)-abb202(64)
      abb202(64)=abb202(56)*spbl4l3
      abb202(65)=abb202(64)*spae1l3
      abb202(66)=-abb202(63)*abb202(65)
      abb202(67)=-abb202(7)*spae1e5
      abb202(68)=-abb202(67)*abb202(17)
      abb202(69)=-abb202(7)*spae2e5
      abb202(70)=abb202(69)*abb202(40)
      abb202(68)=abb202(68)-abb202(70)
      abb202(70)=abb202(59)*spal3l4
      abb202(71)=abb202(70)*spbl3e1
      abb202(72)=-abb202(68)*abb202(71)
      abb202(73)=abb202(13)*spbe5e1
      abb202(74)=abb202(29)*abb202(73)
      abb202(75)=abb202(13)*spbe5e2
      abb202(76)=abb202(75)*abb202(15)
      abb202(74)=abb202(74)+abb202(76)
      abb202(56)=abb202(56)*spbl4k2
      abb202(76)=abb202(56)*spae1k2
      abb202(77)=abb202(74)*abb202(76)
      abb202(78)=abb202(13)*spae1e5
      abb202(79)=abb202(17)*abb202(78)
      abb202(80)=abb202(13)*spae2e5
      abb202(81)=abb202(80)*abb202(40)
      abb202(79)=abb202(79)+abb202(81)
      abb202(59)=abb202(59)*spak2l4
      abb202(81)=abb202(59)*spbk2e1
      abb202(82)=abb202(79)*abb202(81)
      abb202(83)=-abb202(55)*abb202(17)
      abb202(84)=abb202(25)*spae2e5
      abb202(85)=-abb202(40)*abb202(84)
      abb202(83)=abb202(85)+abb202(83)
      abb202(85)=1.0_ki/2.0_ki*spbe2e1
      abb202(83)=abb202(83)*abb202(85)
      abb202(86)=-abb202(50)*abb202(29)
      abb202(25)=abb202(25)*spbe5e2
      abb202(87)=-abb202(15)*abb202(25)
      abb202(86)=abb202(87)+abb202(86)
      abb202(87)=1.0_ki/2.0_ki*spae1e2
      abb202(86)=abb202(86)*abb202(87)
      abb202(19)=abb202(86)+abb202(83)+abb202(82)+abb202(77)+abb202(72)+abb202(&
      &66)+abb202(61)+abb202(58)+abb202(31)+abb202(39)+abb202(23)+abb202(19)+ab&
      &b202(42)+abb202(33)
      abb202(9)=abb202(24)*abb202(9)
      abb202(23)=abb202(9)*abb202(15)
      abb202(12)=abb202(24)*abb202(12)
      abb202(31)=abb202(12)*abb202(15)
      abb202(33)=abb202(23)+abb202(31)
      abb202(39)=abb202(33)*abb202(18)
      abb202(42)=abb202(9)+abb202(12)
      abb202(21)=-abb202(21)*abb202(42)*abb202(20)
      abb202(58)=abb202(36)+abb202(13)
      abb202(61)=abb202(58)*abb202(27)
      abb202(21)=abb202(61)+abb202(39)+abb202(21)
      abb202(21)=spae1k2*abb202(21)
      abb202(39)=abb202(9)*abb202(40)
      abb202(61)=abb202(12)*abb202(40)
      abb202(66)=abb202(39)+abb202(61)
      abb202(72)=-abb202(66)*abb202(28)
      abb202(58)=abb202(58)*abb202(46)
      abb202(58)=abb202(58)+abb202(72)
      abb202(58)=spbk2e1*abb202(58)
      abb202(72)=abb202(6)*abb202(40)
      abb202(77)=-abb202(72)*abb202(29)
      abb202(82)=abb202(67)*spbe5e2
      abb202(77)=abb202(82)+abb202(77)
      abb202(82)=spbl4e1*spae2l4
      abb202(83)=spae2l3*spbl3e1
      abb202(82)=abb202(82)+abb202(83)
      abb202(86)=abb202(77)*abb202(82)
      abb202(88)=abb202(6)*abb202(15)
      abb202(89)=-abb202(88)*abb202(17)
      abb202(90)=abb202(62)*spae2e5
      abb202(89)=abb202(90)+abb202(89)
      abb202(90)=spae1l4*spbl4e2
      abb202(91)=spbl3e2*spae1l3
      abb202(90)=abb202(90)+abb202(91)
      abb202(92)=abb202(89)*abb202(90)
      abb202(93)=3.0_ki*spbe2e1
      abb202(94)=-abb202(68)*abb202(93)
      abb202(95)=3.0_ki*spae1e2
      abb202(96)=-abb202(63)*abb202(95)
      abb202(21)=abb202(96)+abb202(94)+abb202(58)+abb202(21)+abb202(92)+abb202(&
      &86)
      abb202(21)=1.0_ki/2.0_ki*abb202(21)
      abb202(58)=abb202(32)*abb202(34)
      abb202(86)=spbl3k2*spae2k2
      abb202(92)=spae1l3*abb202(10)*abb202(86)
      abb202(58)=abb202(92)+abb202(58)
      abb202(58)=abb202(30)*abb202(58)
      abb202(92)=spbk2e1*spae2k2
      abb202(44)=abb202(44)*abb202(92)
      abb202(52)=-abb202(30)*abb202(52)
      abb202(94)=spbl3e1*abb202(53)
      abb202(52)=abb202(52)+abb202(94)
      abb202(52)=spae2l3*abb202(52)
      abb202(94)=abb202(53)*abb202(57)
      abb202(65)=-abb202(62)*abb202(65)
      abb202(96)=-abb202(73)*abb202(76)
      abb202(50)=abb202(50)*abb202(87)
      abb202(44)=abb202(50)+abb202(96)+abb202(65)+abb202(94)+abb202(52)+abb202(&
      &44)+abb202(58)
      abb202(50)=abb202(72)*abb202(82)
      abb202(52)=abb202(66)*abb202(92)
      abb202(58)=-abb202(62)*abb202(95)
      abb202(50)=abb202(58)+abb202(52)+abb202(50)
      abb202(50)=1.0_ki/2.0_ki*abb202(50)
      abb202(52)=abb202(41)*abb202(43)
      abb202(58)=spak2l3*spbk2e2
      abb202(10)=spbl3e1*abb202(10)*abb202(58)
      abb202(10)=abb202(10)+abb202(52)
      abb202(10)=abb202(38)*abb202(10)
      abb202(52)=spae1k2*spbk2e2
      abb202(16)=abb202(16)*abb202(52)
      abb202(47)=-abb202(38)*abb202(47)
      abb202(65)=spae1l3*abb202(48)
      abb202(47)=abb202(47)+abb202(65)
      abb202(47)=spbl3e2*abb202(47)
      abb202(65)=abb202(48)*abb202(60)
      abb202(66)=-abb202(67)*abb202(71)
      abb202(71)=-abb202(78)*abb202(81)
      abb202(55)=abb202(55)*abb202(85)
      abb202(10)=abb202(55)+abb202(71)+abb202(66)+abb202(65)+abb202(47)+abb202(&
      &16)+abb202(10)
      abb202(16)=abb202(88)*abb202(90)
      abb202(33)=abb202(33)*abb202(52)
      abb202(47)=-abb202(67)*abb202(93)
      abb202(16)=abb202(47)+abb202(33)+abb202(16)
      abb202(16)=1.0_ki/2.0_ki*abb202(16)
      abb202(8)=1.0_ki/2.0_ki*abb202(8)
      abb202(33)=abb202(8)*spbe5e2
      abb202(47)=-abb202(32)*abb202(33)
      abb202(55)=-abb202(7)*spbe5e2
      abb202(64)=abb202(55)*abb202(64)
      abb202(65)=-spbl3k2*abb202(11)*abb202(35)
      abb202(47)=abb202(65)+abb202(64)+abb202(47)
      abb202(47)=spae1l3*abb202(47)
      abb202(20)=abb202(20)*spbe5k1
      abb202(18)=abb202(20)-abb202(18)
      abb202(20)=-abb202(14)*abb202(18)
      abb202(64)=spae2l3*abb202(2)
      abb202(33)=abb202(33)*abb202(64)
      abb202(20)=abb202(33)+abb202(20)
      abb202(20)=spae1k2*abb202(20)
      abb202(33)=abb202(7)*abb202(17)
      abb202(60)=abb202(60)+abb202(91)
      abb202(65)=abb202(33)*abb202(60)
      abb202(66)=abb202(75)*abb202(76)
      abb202(25)=-abb202(87)*abb202(25)
      abb202(20)=abb202(25)+abb202(66)+abb202(47)+abb202(65)+abb202(20)
      abb202(25)=-abb202(6)*abb202(17)
      abb202(47)=abb202(25)*abb202(90)
      abb202(65)=-spae1k2*abb202(42)*abb202(18)
      abb202(66)=abb202(55)*abb202(95)
      abb202(47)=abb202(66)+abb202(65)+abb202(47)
      abb202(47)=1.0_ki/2.0_ki*abb202(47)
      abb202(60)=-abb202(7)*abb202(60)
      abb202(65)=abb202(14)*abb202(52)
      abb202(60)=abb202(65)+abb202(60)
      abb202(52)=abb202(42)*abb202(52)
      abb202(65)=abb202(6)*abb202(90)
      abb202(52)=abb202(52)+abb202(65)
      abb202(52)=1.0_ki/2.0_ki*abb202(52)
      abb202(8)=abb202(8)*spae2e5
      abb202(65)=-abb202(41)*abb202(8)
      abb202(66)=abb202(69)*abb202(70)
      abb202(11)=-spak2l3*abb202(11)*abb202(26)
      abb202(11)=abb202(11)+abb202(66)+abb202(65)
      abb202(11)=spbl3e1*abb202(11)
      abb202(65)=-abb202(14)*abb202(28)
      abb202(66)=spbl3e2*abb202(3)
      abb202(8)=abb202(8)*abb202(66)
      abb202(8)=abb202(8)+abb202(65)
      abb202(8)=spbk2e1*abb202(8)
      abb202(65)=abb202(7)*abb202(29)
      abb202(57)=abb202(57)+abb202(83)
      abb202(70)=abb202(65)*abb202(57)
      abb202(71)=abb202(80)*abb202(81)
      abb202(76)=-abb202(85)*abb202(84)
      abb202(8)=abb202(76)+abb202(71)+abb202(11)+abb202(70)+abb202(8)
      abb202(11)=-abb202(6)*abb202(29)
      abb202(70)=abb202(11)*abb202(82)
      abb202(71)=-spbk2e1*abb202(42)*abb202(28)
      abb202(76)=abb202(69)*abb202(93)
      abb202(70)=abb202(76)+abb202(71)+abb202(70)
      abb202(70)=1.0_ki/2.0_ki*abb202(70)
      abb202(57)=-abb202(7)*abb202(57)
      abb202(14)=abb202(14)*abb202(92)
      abb202(14)=abb202(14)+abb202(57)
      abb202(42)=abb202(42)*abb202(92)
      abb202(57)=abb202(6)*abb202(82)
      abb202(42)=abb202(42)+abb202(57)
      abb202(42)=1.0_ki/2.0_ki*abb202(42)
      abb202(54)=-3.0_ki/2.0_ki*abb202(54)
      abb202(57)=1.0_ki/2.0_ki*abb202(77)
      abb202(53)=-3.0_ki/2.0_ki*abb202(53)
      abb202(71)=1.0_ki/2.0_ki*abb202(72)
      abb202(65)=-3.0_ki/2.0_ki*abb202(65)
      abb202(76)=1.0_ki/2.0_ki*abb202(11)
      abb202(7)=-3.0_ki/2.0_ki*abb202(7)
      abb202(81)=1.0_ki/2.0_ki*abb202(6)
      abb202(49)=-3.0_ki/2.0_ki*abb202(49)
      abb202(82)=1.0_ki/2.0_ki*abb202(89)
      abb202(48)=-3.0_ki/2.0_ki*abb202(48)
      abb202(83)=1.0_ki/2.0_ki*abb202(88)
      abb202(33)=-3.0_ki/2.0_ki*abb202(33)
      abb202(84)=1.0_ki/2.0_ki*abb202(25)
      abb202(85)=-spae2l4*abb202(63)
      abb202(87)=-spae2l4*abb202(62)
      abb202(90)=spae2l4*abb202(55)
      abb202(91)=-spbl4e2*abb202(68)
      abb202(92)=-spbl4e2*abb202(67)
      abb202(93)=spbl4e2*abb202(69)
      abb202(94)=abb202(63)*abb202(34)
      abb202(95)=abb202(62)*abb202(34)
      abb202(34)=-abb202(34)*abb202(55)
      abb202(96)=abb202(68)*abb202(43)
      abb202(97)=abb202(67)*abb202(43)
      abb202(43)=-abb202(69)*abb202(43)
      abb202(98)=abb202(24)*abb202(15)
      abb202(99)=abb202(17)*abb202(98)
      abb202(30)=abb202(30)*spae2e5
      abb202(30)=abb202(30)-abb202(99)
      abb202(99)=abb202(41)*abb202(30)
      abb202(100)=3.0_ki*spae2l3
      abb202(63)=-abb202(63)*abb202(100)
      abb202(101)=spal3l4*spbl4e2
      abb202(89)=-abb202(89)*abb202(101)
      abb202(102)=-abb202(23)*abb202(37)
      abb202(27)=abb202(36)*abb202(27)
      abb202(27)=abb202(27)+abb202(102)
      abb202(27)=spak2l3*abb202(27)
      abb202(27)=abb202(27)+abb202(89)+abb202(63)+abb202(99)
      abb202(27)=1.0_ki/2.0_ki*abb202(27)
      abb202(62)=-3.0_ki/2.0_ki*spae2l3*abb202(62)
      abb202(63)=abb202(98)*abb202(41)
      abb202(88)=-abb202(88)*abb202(101)
      abb202(23)=abb202(23)*abb202(58)
      abb202(23)=abb202(23)+abb202(63)+abb202(88)
      abb202(23)=1.0_ki/2.0_ki*abb202(23)
      abb202(41)=abb202(41)*abb202(24)
      abb202(63)=-abb202(41)*abb202(17)
      abb202(37)=-spak2l3*abb202(9)*abb202(37)
      abb202(55)=abb202(55)*abb202(100)
      abb202(25)=-abb202(25)*abb202(101)
      abb202(25)=abb202(37)+abb202(25)+abb202(55)+abb202(63)
      abb202(25)=1.0_ki/2.0_ki*abb202(25)
      abb202(37)=-abb202(6)*abb202(101)
      abb202(55)=abb202(9)*abb202(58)
      abb202(37)=abb202(55)+abb202(41)+abb202(37)
      abb202(37)=1.0_ki/2.0_ki*abb202(37)
      abb202(41)=abb202(24)*abb202(40)
      abb202(55)=abb202(29)*abb202(41)
      abb202(38)=abb202(38)*spbe5e2
      abb202(38)=abb202(38)-abb202(55)
      abb202(55)=abb202(32)*abb202(38)
      abb202(58)=3.0_ki*spbl3e2
      abb202(63)=-abb202(68)*abb202(58)
      abb202(68)=spbl4l3*spae2l4
      abb202(77)=-abb202(77)*abb202(68)
      abb202(88)=-abb202(39)*abb202(22)
      abb202(46)=abb202(36)*abb202(46)
      abb202(46)=abb202(46)+abb202(88)
      abb202(46)=spbl3k2*abb202(46)
      abb202(46)=abb202(46)+abb202(77)+abb202(63)+abb202(55)
      abb202(46)=1.0_ki/2.0_ki*abb202(46)
      abb202(55)=abb202(41)*abb202(32)
      abb202(63)=-abb202(72)*abb202(68)
      abb202(39)=abb202(39)*abb202(86)
      abb202(39)=abb202(39)+abb202(55)+abb202(63)
      abb202(39)=1.0_ki/2.0_ki*abb202(39)
      abb202(55)=-3.0_ki/2.0_ki*spbl3e2*abb202(67)
      abb202(32)=abb202(32)*abb202(24)
      abb202(63)=-abb202(32)*abb202(29)
      abb202(22)=-spbl3k2*abb202(9)*abb202(22)
      abb202(58)=abb202(69)*abb202(58)
      abb202(11)=-abb202(11)*abb202(68)
      abb202(11)=abb202(22)+abb202(11)+abb202(58)+abb202(63)
      abb202(11)=1.0_ki/2.0_ki*abb202(11)
      abb202(6)=-abb202(6)*abb202(68)
      abb202(9)=abb202(9)*abb202(86)
      abb202(6)=abb202(9)+abb202(32)+abb202(6)
      abb202(6)=1.0_ki/2.0_ki*abb202(6)
      abb202(9)=-abb202(74)*abb202(51)
      abb202(22)=abb202(73)*abb202(51)
      abb202(32)=-abb202(75)*abb202(51)
      abb202(51)=-abb202(79)*abb202(45)
      abb202(58)=abb202(78)*abb202(45)
      abb202(45)=-abb202(80)*abb202(45)
      abb202(13)=abb202(36)+3.0_ki/2.0_ki*abb202(13)
      abb202(36)=abb202(13)*spbe5e1
      abb202(28)=abb202(36)*abb202(28)
      abb202(35)=abb202(13)*abb202(35)
      abb202(15)=abb202(15)*abb202(35)
      abb202(63)=1.0_ki/2.0_ki*spbl3e2
      abb202(30)=-abb202(63)*abb202(3)*abb202(30)
      abb202(67)=abb202(17)*abb202(31)
      abb202(68)=abb202(73)*spae2e5
      abb202(67)=abb202(68)-abb202(67)
      abb202(68)=-abb202(67)*abb202(59)
      abb202(15)=abb202(68)+abb202(30)+abb202(15)+abb202(28)
      abb202(28)=-spae2k2*abb202(36)
      abb202(30)=-abb202(98)*abb202(66)
      abb202(36)=spak2l4*spbl4e2
      abb202(66)=-abb202(31)*abb202(36)
      abb202(30)=abb202(30)+abb202(66)
      abb202(30)=1.0_ki/2.0_ki*abb202(30)
      abb202(66)=abb202(24)*abb202(3)
      abb202(63)=abb202(63)*abb202(66)*abb202(17)
      abb202(17)=abb202(17)*abb202(12)
      abb202(59)=abb202(17)*abb202(59)
      abb202(35)=abb202(59)+abb202(35)+abb202(63)
      abb202(59)=-spbl3e2*abb202(66)
      abb202(36)=-abb202(12)*abb202(36)
      abb202(36)=abb202(59)+abb202(36)
      abb202(36)=1.0_ki/2.0_ki*abb202(36)
      abb202(59)=1.0_ki/2.0_ki*abb202(67)
      abb202(31)=1.0_ki/2.0_ki*abb202(31)
      abb202(17)=-1.0_ki/2.0_ki*abb202(17)
      abb202(63)=1.0_ki/2.0_ki*abb202(12)
      abb202(66)=abb202(13)*spae1e5
      abb202(18)=abb202(66)*abb202(18)
      abb202(13)=abb202(13)*abb202(26)
      abb202(26)=abb202(40)*abb202(13)
      abb202(40)=1.0_ki/2.0_ki*spae2l3
      abb202(38)=-abb202(40)*abb202(2)*abb202(38)
      abb202(67)=abb202(29)*abb202(61)
      abb202(68)=abb202(78)*spbe5e2
      abb202(67)=abb202(68)-abb202(67)
      abb202(68)=-abb202(67)*abb202(56)
      abb202(18)=abb202(68)+abb202(38)+abb202(26)+abb202(18)
      abb202(26)=-abb202(41)*abb202(64)
      abb202(38)=spbl4k2*spae2l4
      abb202(41)=-abb202(61)*abb202(38)
      abb202(26)=abb202(26)+abb202(41)
      abb202(26)=1.0_ki/2.0_ki*abb202(26)
      abb202(41)=-spbk2e2*abb202(66)
      abb202(24)=abb202(24)*abb202(2)
      abb202(40)=abb202(40)*abb202(24)*abb202(29)
      abb202(29)=abb202(29)*abb202(12)
      abb202(56)=abb202(29)*abb202(56)
      abb202(13)=abb202(56)+abb202(13)+abb202(40)
      abb202(24)=-spae2l3*abb202(24)
      abb202(12)=-abb202(12)*abb202(38)
      abb202(12)=abb202(24)+abb202(12)
      abb202(12)=1.0_ki/2.0_ki*abb202(12)
      abb202(24)=1.0_ki/2.0_ki*abb202(67)
      abb202(38)=1.0_ki/2.0_ki*abb202(61)
      abb202(29)=-1.0_ki/2.0_ki*abb202(29)
      R2d202=0.0_ki
      rat2 = rat2 + R2d202
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='202' value='", &
          & R2d202, "'/>"
      end if
   end subroutine
end module p2_part21part21_part25part25part21_abbrevd202h0
