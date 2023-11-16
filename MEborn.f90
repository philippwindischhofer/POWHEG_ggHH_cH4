! adapted from gitrepo/HH
!subroutine ME2born_gbij(p, mh2, mt2, muren2, chhh, modg1, modg2)
module Born_amplitudes
  use precision_golem
  implicit none
  private

  public :: ME2born_gbij, compPolVectors, physics_initialise_f90
!  public :: eps
!-----------------------------------------------------------------------  
!>> variables for F1quartic form factor:
  complex(ki), parameter:: czip=(0.0_ki, 0.0_ki)
  complex(ki), parameter:: cone=(1.0_ki, 0.0_ki)
  complex(ki), parameter:: imag=(0.0_ki, 1.0_ki)
  real(ki), parameter:: pi=3.1415926535897932384626433832795_ki
  integer, parameter:: maxbin=999
  real(ki), public:: sqrts(maxbin), f1real(maxbin), f1imag(maxbin)
  !>> physics parameters
  real(ki), public:: f90_alphaem,f90_wmass,f90_sthw2,f90_cH4

  contains

subroutine ME2born_gbij(p, mh2, mt2, muren2, chhh, ct, ctt, cg, cgg, mpol, imode)
  !o
  use precision_golem
  use parametre
  use matrice_s  
  use form_factor_type
  use form_factor_4p
  use form_factor_3p
  use cache
  use constante  
  use array
  use spinor, only: scalar  
  !
  !
  implicit none
  !
  !
  real(ki) :: s12,s23,s13,s1,s2,s3,s4,mw,mh2,mt2,muren2
  real(ki) :: m1sq,m2sq,m3sq,m4sq
  integer:: imode !>> this is to identify if Born was called from 'setborn' or 'setvirt'
!  complex(ki) :: m1sq,m2sq,m3sq,m4sq
  real(ki) :: coef,smin,smax,full,modg1,modg2,ME2
  real(ki), dimension(4) :: mass,mass_int_sq
  real(ki), dimension(4) :: p1,p2,p3,p4,p5,p6,k1,k2,k3,k4,p12,p13,p23,p34
  real(ki), dimension(4,4) :: p
  real(ki), dimension(6,4) :: vecs
  complex(ki), dimension(-1:1,-1:1) :: mpol
  integer :: i,j,k,ep,iep,nbpoints,nlegs,jmax,ievt,nsteps,l,lp
  complex(ki) :: cs1,cs2,cs3,cm0,cm1,cm2,cm3,cm4,gauge1,gauge2
  complex(ki) :: Gquartic(2)
  character (len=3) :: bb0,bb1
  type(form_factor) :: D123,D132,D213,C12,C13,C23,C34,D0temp
  type(form_factor) :: t0,gauge1tri,gauge1box,gauge2box,g1trilim,g1boxlim,g2boxlim
  real(ki) :: chhh,ct,ctt,cg,cgg
  complex*16, external :: D04
  real*8 :: mt
  logical*1 HTL
  logical, parameter:: wolfram_mathematica_check=.false.
  real(ki):: factor
  real(ki):: lambda, sqrttau, tauh, partreal, partimag
  complex(ki):: correction

  !
  ! delta=1.e-8_ki
  !
  mu2_scale_par=muren2

  s1=    0._ki
  s2=    0._ki
  s3=    mh2
  s4=    mh2
  m1sq=   mt2
  m2sq=   mt2
  m3sq=   mt2  
  m4sq=   mt2  
  mt=sqrt(mt2)

  k1 = p(:,1)
  k2 = p(:,2)
  k3 = -p(:,3)
  k4 = -p(:,4)

  p12=k1+k2
  p13=k1+k3
  p23=k2+k3

  ! 
  call initgolem95(4)

  s12 = 2.0_ki*scalar(k1,k2)
  s13 = scalar(p13,p13)
  s23 = scalar(p23,p23)

  D123 = D04(0d0,0d0,mh2,mh2,s12,s23,mt,mt,mt,mt)
  D213 = D04(0d0,0d0,mh2,mh2,s12,s13,mt,mt,mt,mt)
  D132 = D04(0d0,mh2,0d0,mh2,s23,s13,mt,mt,mt,mt)

  !if(.not.HTL) then

  ! print*, "at the beginning"
  ! print*, "s12=", s12
  ! print*, "s13=", s13
  ! print*, "s23=", s23
  ! print*, "mh2=", mh2
  ! print*, "mt2=", mt2

  
  ! "planar" box:
  s_mat(1,:) = (/-m1sq*2._ki,s2-m1sq-m2sq,s23-m1sq-m3sq,s1-m1sq-m4sq/)
  s_mat(2,:) = (/s2-m1sq-m2sq,-m2sq*2._ki,s3-m2sq-m3sq,s12-m2sq-m4sq/)
  s_mat(3,:) = (/s23-m1sq-m3sq,s3-m2sq-m3sq,-m3sq*2._ki,s4-m3sq-m4sq/)
  s_mat(4,:) = (/s1-m1sq-m4sq,s12-m2sq-m4sq,s4-m3sq-m4sq,-m4sq*2._ki/)
  !
  !
  
  call preparesmatrix()
  !
  !D0temp=A40(s_null)
  !write(6,*) 'D123=',D123,D0temp

  C34=A30((/1/))
  C23=A30((/2/))
  C12=A30((/3/))

  call exitgolem95()

! crossed graphs: 1 <-> 2

  !call initgolem95(4)

  !s12 = 2.0_ki*scalar(k1,k2)
  !s23 = scalar(p13,p13)
  !s13 = scalar(p23,p23)

  !! box D213:
  !s_mat(1,:) = (/-m1sq*2._ki,s2-m1sq-m2sq,s23-m1sq-m3sq,s1-m1sq-m4sq/)
  !s_mat(2,:) = (/s2-m1sq-m2sq,-m2sq*2._ki,s3-m2sq-m3sq,s12-m2sq-m4sq/)
  !s_mat(3,:) = (/s23-m1sq-m3sq,s3-m2sq-m3sq,-m3sq*2._ki,s4-m3sq-m4sq/)
  !s_mat(4,:) = (/s1-m1sq-m4sq,s12-m2sq-m4sq,s4-m3sq-m4sq,-m4sq*2._ki/)

  !!
  !call preparesmatrix()
  !!
  !!D0temp=A40(s_null)
  !!write(6,*) 'D213=',D213, D0temp


  !call exitgolem95()

  ! crossed graphs 2 <-> 3

  call initgolem95(4)

  s13 = scalar(p12,p12)
  s12 = scalar(p13,p13)
  s23 = scalar(p23,p23)

  ! box D132:
  s_mat(1,:) = (/-m1sq*2._ki,s3-m1sq-m2sq,s23-m1sq-m3sq,s1-m1sq-m4sq/)
  s_mat(2,:) = (/s3-m1sq-m2sq,-m2sq*2._ki,s2-m2sq-m3sq,s12-m2sq-m4sq/)
  s_mat(3,:) = (/s23-m1sq-m3sq,s2-m2sq-m3sq,-m3sq*2._ki,s4-m3sq-m4sq/)
  s_mat(4,:) = (/s1-m1sq-m4sq,s12-m2sq-m4sq,s4-m3sq-m4sq,-m4sq*2._ki/)

  !
  call preparesmatrix()
  !
  !D0temp=A40(s_null)
  !write(6,*) 'D132=',D132, D0temp
  C13=A30((/3/))

  call exitgolem95()

!!!!!! amplitudes !!!!!!!!!!

  s12 = scalar(p12,p12)
  s13 = scalar(p13,p13)
  s23 = scalar(p23,p23)

  if( imode.eq.0 )then
     !>> include kappa3 bit

     !gauge1tri=12._ki*mh2*m1sq/(s12-mh2)*(2._ki+(4._ki*m1sq-s12)*C12)*(ct*chhh &
     !     &  + 2._ki/3*(s12-mh2)/mh2*ctt) + 6._ki*s12*mh2/(s12-mh2)*cg*chhh

     !>> a correction due to 1-loop diagrams for three-Higgs-vertex [ numerical expresstion valid for: 0 < tauh < 1 ]
     lambda = f90_alphaem * pi * mh2 / (2d0 * f90_Wmass**2 * f90_sthw2 )
     tauh= 4d0*mh2/s12
     sqrttau= sqrt(abs( 1d0 - tauh ))
     !>> eq. 3.8, 3.9
     partreal= + f90_cH4 * lambda * ( -18d0 + 2d0*sqrt(3d0)*pi + 3*sqrttau*log((1d0+sqrttau)/(1d0-sqrttau)) ) / (16d0*pi**2)
     partimag= - f90_cH4 * lambda * 3*sqrttau * pi / (16d0*pi**2)
     !>> eq. 3.10, 3.11
     partreal= partreal + f90_cH4 * lambda / chhh * ( -7d0 ) / (16d0*pi**2)
     correction = cone + partreal*cone + partimag*imag

     gauge1tri=12._ki*mh2*m1sq/(s12-mh2)*(2._ki+(4._ki*m1sq-s12)*C12)*(ct*chhh*correction &
          &  + 2._ki/3*(s12-mh2)/mh2*ctt) + 6._ki*s12*mh2/(s12-mh2)*cg*chhh
  else
     !>> do NOT include kappa3 bit for counterterms
     gauge1tri=12._ki*mh2*m1sq/(s12-mh2)*(2._ki+(4._ki*m1sq-s12)*C12)*(ct*1d0 &
          &  + 2._ki/3*(s12-mh2)/mh2*ctt) + 6._ki*s12*mh2/(s12-mh2)*cg*1d0
  endif


  gauge1box= 4._ki*m1sq*ct**2*( m1sq*(8._ki*m1sq-s12-2._ki*mh2)*(D123+D213+D132) + &
       &  (s13*s23-mh2**2)/s12*(-mh2+4._ki*m1sq)*D132+2._ki+4._ki*m1sq*C12 &
       &  +2._ki/s12*(mh2-4._ki*m1sq)*((s13-mh2)*C13+(s23-mh2)*C23) ) + 4._ki*cgg*s12

  gauge2box= 2._ki*m1sq*ct**2*( 2._ki*(8._ki*m1sq+s12-2._ki*mh2)*(m1sq*(D123+D213+D132)-C34) &
       &  -2._ki*(s12*C12+(s13-mh2)*C13+(s23-mh2)*C23)  &
       &  + 1._ki/(s13*s23-mh2**2)*( &
       &  + s12*s23*(8._ki*s23*m1sq-s23**2-mh2**2)*D123 &
       &  + s12*s13*(8._ki*s13*m1sq-s13**2-mh2**2)*D213 &
       &  + (8._ki*m1sq+s12-2._ki*mh2)*( s12*(s12-2._ki*mh2)*C12 + s12*(s12-4._ki*mh2)*C34 &
       &  + 2._ki*s13*(mh2-s13)*C13 + 2._ki*s23*(mh2-s23)*C23 ) ) )

  ! g1trilim= 4._ki*mh2*s12/(s12-mh2)
  ! g1boxlim= -4._ki/3._ki*s12
  ! g2boxlim= -11._ki/45._ki*(s13*s23-mh2**2)/m1sq

!  else ! (HTL)

!  gauge1tri = 4._ki*mh2*s12/(s12-mh2)*chhh
!  gauge1box = -4._ki/3._ki*s12
!  gauge2box = 0._ki

!  endif

  gauge1=(gauge1tri%c+gauge1box%c)/2._ki
  gauge2=gauge2box%c/2._ki

  !!print*,'sqrt(s)=',sqrt(s12)
  if( imode.eq.0 )then 
     !>> include kappa4 bit
     call HHxquartic(sqrt(s12),Gquartic)
  else 
     !>> do not include kappa4 bit
     Gquartic=0.0_ki
  endif
  
  factor= (f90_alphaem/f90_Wmass**2/f90_sthw2)/pi * (mh2*mt2)
  factor= factor * (-2d0)  !!>> to match PWHG vs. FormCalc (for Born ggHH amplitude)

  !!>> the ratio between Luca+Wojtek and Uli (note Uli's numbers are for 1/3*fhat ) 
  factor= factor*(3d0/32d0)

  factor= factor*(f90_cH4)
  
  Gquartic= factor*Gquartic

  if( wolfram_mathematica_check )then
     !!>> this is to check compatibility with FormCalc
     print*,'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
     print*,'k1=',k1
     print*,'k2=',k2
     print*,'k3=',k3
     print*,'k4=',k4
     print*,''
     print*,'s=',s12
     print*,'t=',s13
     print*,'u=',s23
     print*,'>> check:  s+t+u=',s12+s13+s23
     print*,'>> check: 2*mh^2=',2d0*mh2
     print*,''
     print*,'gauge1=',gauge1/(-2d0),'   <--- kappa4=',0d0
     print*,'gauge2=',gauge2/(-2d0),'   <--- kappa4=',0d0
     gauge1= gauge1 + Gquartic(1)
     gauge2= gauge2 + Gquartic(2)
     print*,''
     print*,'quartic=',Gquartic(:)
     print*,''
     print*,'gauge1=',gauge1/(-2d0),'   <--- kappa4=',10d0
     print*,'gauge2=',gauge2/(-2d0),'   <--- kappa4=',10d0
     print*,''
     print*,''
     stop
  else
     gauge1= gauge1 + Gquartic(1)
     gauge2= gauge2 + Gquartic(2)
  endif
 
  mpol = 0.0_ki

  do l=-1,1,2
     do lp=-1,1,2
        if(l*lp.eq.1)  then
           mpol(l,lp)=-gauge1
        elseif(l*lp.eq.-1) then
           mpol(l,lp)=-gauge2
        else
           write(*,*) "problem with polarized born"
           call pwhg_exit(-1)
        endif
     enddo
  enddo

end subroutine ME2born_gbij


subroutine compPolVectors(p, epsten)
!  use p2_part21part21_part25part25part21_kinematics, only: epsi
  use p0_part21part21_part25part25_kinematics, only: epsi
  implicit none
  real(ki), intent(in) :: p(4,4)
  complex(ki), intent(out) :: epsten(2, 0:3, -1:1)
  complex(ki) :: eps1(0:3, -1:1), eps2(0:3, -1:1)

  eps1(:,-1) = epsi(p(:,1),p(:,2),1)
  eps1(:,1)  = epsi(p(:,1),p(:,2),-1)
  eps2(:,-1) = epsi(p(:,2),p(:,1),1)
  eps2(:,1)  = epsi(p(:,2),p(:,1),-1)

  epsten(1,:,:)=eps1
  epsten(2,:,:)=eps2

  return
end subroutine compPolVectors

subroutine HHxquartic(xsqrt,gauge)
  implicit none

  logical iniF1
  data iniF1/.true./
  save iniF1

  complex(ki):: gauge(2)
  real(ki):: tau,xsqrt

  integer k,nbins,ios,ibin
!  parameter (maxbin=999)
  character*200 xline
!  real(ki):: sqrts(maxbin),f1real(maxbin),f1imag(maxbin)
!  common/F1TABLE/nbins,sqrts,f1real,f1imag

  if( iniF1 )then
     open(unit=77, file='form-factor-1.dat',status='old')
     do k=1,maxbin
        read(unit=77,fmt='(a)',end=111) xline
        read(unit=xline,fmt=*,iostat=ios) sqrts(k), f1real(k), f1imag(k)
     enddo
111  nbins=k-1
     print*,'## Additional FORM FACTOR:'
     do k=1,nbins
        print*,'>>>',sqrts(k), f1real(k), f1imag(k)
     enddo
     iniF1=.false.
  endif

  gauge(1)= czip !>> zero, placeholder for F1
  gauge(2)= czip !>> zero, placeholder for F2

  if( xsqrt.lt.sqrts(1) ) return
  if( xsqrt.gt.sqrts(nbins) ) return

  do k=1,nbins
     ibin=k
     if( xsqrt.lt.sqrts(k) ) exit !!>> exit loop
  enddo

  tau= ( xsqrt - sqrts(ibin-1) ) / ( sqrts(ibin) - sqrts(ibin-1) )
  gauge(1)= gauge(1) + cone*( f1real(ibin)*tau + f1real(ibin-1)*(1d0-tau) )
  gauge(1)= gauge(1) + imag*( f1imag(ibin)*tau + f1imag(ibin-1)*(1d0-tau) )

  return
end subroutine HHxquartic

subroutine physics_initialise_f90(alphaem,wmass,sthw2,cH4)
  implicit none
  real(ki):: alphaem,wmass,sthw2,cH4
  
  f90_alphaem= alphaem
  f90_wmass=   wmass
  f90_sthw2=   sthw2
  f90_cH4=     cH4

  write(*,*)'======================================================='
  write(*,*)'>> Variables rewritten to f90:'
  write(*,*)'   --> f90_alphaem=',f90_alphaem
  write(*,*)'   --> f90_wmass=  ',f90_wmass
  write(*,*)'   --> f90_sthw2=  ',f90_sthw2
  write(*,*)'   --> f90_ch4=    ',f90_cH4
  write(*,*)'======================================================='
  return
end subroutine physics_initialise_f90

end module Born_amplitudes
