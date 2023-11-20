! adapted from gitrepo/HH
!subroutine ME2born_gbij(p, mh2, mt2, muren2, chhh, modg1, modg2)
module Born_amplitudes
  use precision_golem
  implicit none
  private

  public :: ME2born_gbij, compPolVectors
!  public :: eps
   
  contains

subroutine ME2born_gbij(p, mh2, mt2, muren2, chhh, ct, ctt, cg, cgg, mpol)
  !
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
  real(ki) :: s12,s23,s13,s1,s2,s3,s4,mw,mh2,mt2,muren2
  real(ki) :: m1sq,m2sq,m3sq,m4sq
!  complex(ki) :: m1sq,m2sq,m3sq,m4sq
  real(ki) :: coef,smin,smax,full,modg1,modg2,ME2
  real(ki), dimension(4) :: mass,mass_int_sq
  real(ki), dimension(4) :: p1,p2,p3,p4,p5,p6,k1,k2,k3,k4,p12,p13,p23,p34
  real(ki), dimension(4,4) :: p
  real(ki), dimension(6,4) :: vecs
  complex(ki), dimension(-1:1,-1:1) :: mpol
  integer :: i,j,k,ep,iep,nbpoints,nlegs,jmax,ievt,nsteps,l,lp
  complex(ki) :: cs1,cs2,cs3,cm0,cm1,cm2,cm3,cm4,gauge1,gauge2
  character (len=3) :: bb0,bb1
  type(form_factor) :: D123,D132,D213,C12,C13,C23,C34,D0temp
  type(form_factor) :: t0,gauge1tri,gauge1box,gauge2box,g1trilim,g1boxlim,g2boxlim
  real(ki) :: chhh,ct,ctt,cg,cgg
  complex*16, external :: D04
  real*8 :: mt
  logical*1 HTL

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

  gauge1tri=12._ki*mh2*m1sq/(s12-mh2)*(2._ki+(4._ki*m1sq-s12)*C12)*(ct*chhh &
       &  + 2._ki/3*(s12-mh2)/mh2*ctt) + 6._ki*s12*mh2/(s12-mh2)*cg*chhh

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

end module Born_amplitudes
