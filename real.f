      subroutine setreal(p,rflav,amp2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_kn.h'
      integer nlegs
      parameter (nlegs=nlegreal)
      real * 8 p(0:3,nlegs), p_amp(0:3,nlegs)
      integer rflav(nlegs)
      real * 8 amp2,amp21,amp22
      integer i,j, nu
      integer mtdep
      common/cmtdep/mtdep
      logical ini,rescue
      data ini/.true./
      save ini,rescue
      real *8 powheginput
      external powheginput

      if (ini) then
         rescue=powheginput("#rescue_reals").ne.0d0
         ini=.false.
      endif

      amp2 = 0d0

c******************************************************
c     Choose the type of computation:
c     mtdep = 0: heavy top limit (htl)
c     mtdep = 1: Born improved heavy top limit (full theory born)
c     mtdep = 2: full theory approx (with htl virtual)
c     mtdep = 3: full theory
c     mtdep = 4,5: checks
c******************************************************

C--  Real Cross section in the large top mass limit
      if((mtdep.eq.0).or.(mtdep.eq.1)) then
         call compreal_htl(p,rflav,amp2)

C--  Real Cross section with finite top mass
      else if((mtdep.eq.2).or.(mtdep.eq.3)) then
         call compreal_top(p,rflav,amp2)

C--  Testing purposes
      else if(mtdep.eq.4) then
         call compreal_htl(p,rflav,amp21)
         call compreal_top(p,rflav,amp22)
         write(*,*) "--> Real: should be one:", amp22/amp21
         amp2=amp22
c         pause
      else
         write(*,*) "Unknown value of 'largemtlim' for real, abort!"
         call exit(1)
      endif
      end

      subroutine compreal_htl(p,rflav,amp2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      integer nlegs
      parameter (nlegs=nlegreal)
      real * 8 p(0:3,nlegs), p_amp(0:3,nlegs)
      integer rflav(nlegs)
      real * 8 amp2
      integer i,j

      i=rflav(1)
      j=rflav(2)

      p_amp(:,1) = p(:,1)
      p_amp(:,2) = p(:,2)
      p_amp(:,3) = p(:,5)
      p_amp(:,4) = p(:,3)
      p_amp(:,5) = p(:,4)

      if((i.eq.j).and.(i.eq.0)) then
c     g g -> h g
         call ME2realgg_htl(p_amp,st_alpha/pi,amp2)
      elseif(((i*j.eq.0).and.(i.ne.j)).and.(i.ne.0)) then
c     q g -> h q
         call ME2realqg_htl(p_amp,st_alpha/pi,amp2)
      elseif(((i*j.eq.0).and.(i.ne.j)).and.(j.ne.0)) then
c     g q -> h q
         call ME2realgq_htl(p_amp,st_alpha/pi,amp2)
      elseif((i+j.eq.0).and.(i.ne.j)) then
c     q aq -> h g
         call ME2realqqbar_htl(p_amp,st_alpha/pi,amp2)
      else
         write(*,*) 'ERROR setreal: unammissible flavour!'
         write(*,*) rflav(1),' ',rflav(2),'->',rflav(3),' ',rflav(4)
         call exit(1)
      endif

      if(amp2.eq.0d0) then
         write(*,*) 'WARNING setreal: returning 0 amplitude!'
         write(*,*) rflav(1),' ',rflav(2),'->',rflav(3),' ',rflav(4)
      endif

c     Cancel as/(2pi) associated with amp2. It will be put back by real_ampsq
      amp2 = amp2/(st_alpha/(2d0*pi))

      end subroutine compreal_htl


      subroutine compreal_top(p,rflav,amp2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      integer nlegs
      parameter (nlegs=nlegreal)
c     input
      real * 8 p(0:3,nlegs)
      integer rflav(nlegs)
c     output
      real * 8 amp2
c     additional variables
      integer rflavloc(nlegs), rflav_gosam(1:nlegs,0:6)
      logical equalintlists
      external equalintlists
      integer dim_mom_array,processid,i
      parameter (dim_mom_array=50)
      real * 8 pgosam(dim_mom_array)
      real * 8 params(10),muren,res(4)
      real :: start1, finish1, start2, finish2
      integer icheck
      data(rflav_gosam(i,    0),i=1,nlegs)/
     $      0,
     $      0,
     $     25,
     $     25,
     $      0/
      data(rflav_gosam(i,    1),i=1,nlegs)/
     $      1,
     $      0,
     $     25,
     $     25,
     $      1/
      data(rflav_gosam(i,    2),i=1,nlegs)/
     $      0,
     $      1,
     $     25,
     $     25,
     $      1/
      data(rflav_gosam(i,    3),i=1,nlegs)/
     $     -1,
     $      0,
     $     25,
     $     25,
     $     -1/
      data(rflav_gosam(i,    4),i=1,nlegs)/
     $      0,
     $     -1,
     $     25,
     $     25,
     $     -1/
      data(rflav_gosam(i,    5),i=1,nlegs)/
     $      1,
     $     -1,
     $     25,
     $     25,
     $      0/
      data(rflav_gosam(i,    6),i=1,nlegs)/
     $     -1,
     $      1,
     $     25,
     $     25,
     $      0/


c     transfer the flavor list into a list with only down quark or antiquark
      rflavloc=rflav
      do i=1,nlegs
         if ((rflavloc(i).ne.0).and.(rflavloc(i).ne.25)) then
            if (abs(rflavloc(i)).ne.1) then
               rflavloc(i)=sign(1,rflavloc(i))
            endif
         endif
      enddo

      do i=0,6
         if (equalintlists(nlegs,rflavloc,rflav_gosam(1,i))) then
            processid=i+2
            goto 222
         endif
      enddo
 111  write(*,*) 'NO matching flavour string between POWHEG and GoSam'
      write(*,*) 'PROGRAM ABORT'
      call exit(1)

 222  call gosam_momenta_r(p,pgosam)
      muren=sqrt(st_muren2)
      params(1)=1d0
c      call cpu_time(start1)
      call OLP_EvalSubProcess(processid,pgosam,muren,params,res,icheck)
c      call cpu_time(finish1)
      amp2=res(3)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     GOSAM returns this result with NO gs factor ==>
C     virt_gosam ->  virt_gosam * (gs^2)^AlphasPower =
C     virt_gosam * (4*pi*st_alpha)^AlphasPower
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      amp2=amp2  * (4d0*pi*st_alpha)**3
c     The libraries for the  scalar integrals need an extra factor of 1/(16*pi^2)
c     Since this contribution is the square of a one-loop calculation, this factor has to be squared
     $          * 2d0/(16*pi**2)**2

c     add SM/BSM interference and tree contributions
      processid=processid+7
      call gosam_momenta_r(p,pgosam)
      muren=sqrt(st_muren2)
      params(1)=1d0
c      call cpu_time(start2)
      call OLP_EvalSubProcess(processid,pgosam,muren,params,res,icheck)
c      call cpu_time(finish2)
      amp2 = amp2 + res(4) * (4d0*pi*st_alpha)**3
     $                     * 1d0/(4d0*pi)**3
     $            + res(3) * (4d0*pi*st_alpha)**3
     $                     * 1d0/(4d0*pi)**3

c    Check processing time (SM x BSM automatic rescue if icheck = 3)
c
c      write(*,*) '#####################################################'
c      write(*,*) "Time elapsed SM: ", finish1-start1
c      write(*,*) "Time elapsed SM x BSM: ", finish2-start2
c     
c      endif
c      write(*,*) '(alpha_EW = ', ph_alphaem, ')'
c      write(*,*) '(alpha = ', st_alpha, ')'

c     cancel as/(2pi) associated with amp2. It will be put back by real_ampsq
      amp2 = amp2/(st_alpha/(2d0*pi))

      end subroutine compreal_top


      subroutine gosam_momenta_r(p,pgosam)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      real * 8 p(0:3,nlegreal)
c     In GoSam the array of the momenta has dimension 50.
c     It accounts for 10 particles at most
      integer dim_mom_array
      parameter (dim_mom_array=50)
      real * 8 pgosam(dim_mom_array)
      integer i

      if (nlegborn*5 .gt. dim_mom_array) then
         write(*,*) 'The dimension of the pgosam array in the '//
     $        'pwhg_gosam.f file NEEDS to be increased'
         write(*,*) 'PROGRAM ABORT'
         call exit(1)
      endif

      do i=1,nlegreal
         pgosam(1+5*(i-1))=p(0,i)
         pgosam(2+5*(i-1))=p(1,i)
         pgosam(3+5*(i-1))=p(2,i)
         pgosam(4+5*(i-1))=p(3,i)
c         write(*,*) i,p(0,i)**2-p(1,i)**2-p(2,i)**2-p(3,i)**2
         pgosam(5+5*(i-1))=kn_masses(i)
      enddo
      end

      subroutine ME2realgg_htl(p, aopi, ME2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'PhysPars.h'
      include 'pwhg_math.h'
      integer nlegs
      parameter (nlegs=nlegreal)
c     input:
      real * 8 p(0:3,nlegs)
      real * 8 aopi
c     output:
      real * 8 ME2
c     invariants, abbreviations:
      real * 8 square_Lorentz_vector
      external square_Lorentz_vector
      real * 8 s,t,u,sHH,s14,s24,s34,denH,MH2
      real * 8 NA,nf
      parameter(NA=8.0d0)

      nf = st_nlight
      MH2 = ph_Hmass2

      s = square_Lorentz_vector(p(:,1)+p(:,2))
      t = square_Lorentz_vector(p(:,1)-p(:,3))
      u = square_Lorentz_vector(p(:,2)-p(:,3))
      sHH = square_Lorentz_vector(p(:,4)+p(:,5))
      s14 = square_Lorentz_vector(p(:,1)-p(:,4))
      s24 = square_Lorentz_vector(p(:,2)-p(:,4))
      s34 = square_Lorentz_vector(p(:,3)+p(:,4))
      denH = 1d0/(sHH-MH2)*ph_mdlchhh

      ME2 =  + NA*CA*TF**2 * ( 4.0d0/9.0d0*s**(-1)*t**(-1)*u**3 - 8.0d0
     $     /3.0d0*s**(-1)*t**(-1)*u**3*MH2*denH + 4.0d0*s**(-1)*t**(-1)
     $     *u**3*MH2**2*denH**2 + 8.0d0/9.0d0*s**(-1)*u**2 - 16.0d0
     $     /3.0d0*s**(-1)*u**2*MH2*denH + 8.0d0*s**(-1)*u**2*MH2**2*denH
     $     **2 + 4.0d0/3.0d0*s**(-1)*t*u - 8.0d0*s**(-1)*t*u*MH2*denH +
     $     12.0d0*s**(-1)*t*u*MH2**2*denH**2 + 8.0d0/9.0d0*s**(-1)*t**2
     $     - 16.0d0/3.0d0*s**(-1)*t**2*MH2*denH + 8.0d0*s**(-1)*t**2*MH2
     $     **2*denH**2 + 4.0d0/9.0d0*s**(-1)*t**3*u**(-1) - 8.0d0/3.0d0
     $     *s**(-1)*t**3*u**(-1)*MH2*denH + 4.0d0*s**(-1)*t**3*u**(-1)
     $     *MH2**2*denH**2 + 8.0d0/9.0d0*t**(-1)*u**2 - 16.0d0/3.0d0*t
     $     **(-1)*u**2*MH2*denH + 8.0d0*t**(-1)*u**2*MH2**2*denH**2 +
     $     8.0d0/3.0d0*u - 16.0d0*u*MH2*denH +24.0d0*u*MH2**2*denH**2 +
     $     8.0d0/3.0d0*t - 16.0d0*t*MH2*denH + 24.0d0*t*MH2**2*denH**2 +
     $     8.0d0/9.0d0*t**2*u**(-1)- 16.0d0/3.0d0*t**2*u**(-1)*MH2*denH
     $     + 8.0d0*t**2*u**(-1)*MH2**2*denH**2 + 4.0d0/3.0d0*s*t**(-1)*u
     $     - 8.0d0*s*t**(-1)*u*MH2*denH + 12.0d0*s*t**(-1)*u*MH2**2*denH
     $     **2 + 8.0d0/3.0d0*s - 16.0d0*s*MH2*denH + 24.0d0*s*MH2**2
     $     *denH**2+ 4.0d0/3.0d0*s*t*u**(-1) - 8.0d0*s*t*u**(-1)*MH2
     $     *denH +12.0d0*s*t*u**(-1)*MH2**2*denH**2 + 8.0d0/9.0d0*s**2*t
     $     **(-1) - 16.0d0/3.0d0*s**2*t**(-1)*MH2*denH + 8.0d0*s**2*t**(
     $     -1)*MH2**2*denH**2 + 8.0d0/9.0d0*s**2*u**(-1) - 16.0d0/3.0d0
     $     *s**2*u**(-1)*MH2*denH + 8.0d0*s**2*u**(-1)*MH2**2*denH**2 +
     $     4.0d0/9.0d0*s**3*t**(-1)*u**(-1) - 8.0d0/3.0d0*s**3*t**(-1)*u
     $     **(-1)*MH2*denH + 4.0d0*s**3*t**(-1)*u**(-1)*MH2**2*denH**2 )

      ME2 = ME2*aopi**3*ph_GF**2
     $     *16.d0*pi**2         ! from gs^2 = 4*pi*alpha_s
     $     /2.0d0               ! Symmetry factor
     $     /256.0d0             ! spin/colour average

      end subroutine ME2realgg_htl


      subroutine ME2realqg_htl(p, aopi, ME2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'PhysPars.h'
      include 'pwhg_math.h'
      integer nlegs
      parameter (nlegs=nlegreal)
c     input:
      real * 8 p(0:3,nlegs)
      real * 8 aopi
c     output:
      real * 8 ME2
c     invariants, abbreviations:
      real * 8 square_Lorentz_vector
      external square_Lorentz_vector
      real * 8 s,t,u,sHH,s14,s24,s34,denH,MH2
      real * 8 NA,nf
      parameter(NA=8.0d0)

      MH2 = ph_Hmass2

      s = square_Lorentz_vector(p(:,1)+p(:,2))
      t = square_Lorentz_vector(p(:,1)-p(:,3))
      u = square_Lorentz_vector(p(:,2)-p(:,3))
      sHH = square_Lorentz_vector(p(:,4)+p(:,5))
      s14 = square_Lorentz_vector(p(:,1)-p(:,4))
      s24 = square_Lorentz_vector(p(:,2)-p(:,4))
      s34 = square_Lorentz_vector(p(:,3)+p(:,4))
      denH = 1d0/(sHH-MH2)*ph_mdlchhh

      ME2 =  + NA*TF**3 * (  - 2.0d0/9.0d0*t**(-1)*u**2 + 4.0d0/3.0d0*t
     $     **(-1)*u**2*MH2*denH - 2.0d0*t**(-1)*u**2*MH2**2*denH**2 -
     $     2.0d0/9.0d0*s**2*t**(-1) + 4.0d0/3.0d0*s**2*t**(-1)*MH2* denH
     $     - 2.0d0*s**2*t**(-1)*MH2**2*denH**2 )

      ME2 = ME2*aopi**3*ph_GF**2
     $     *16.d0*pi**2         ! from gs^2 = 4*pi*alpha_s
     $     /2.0d0               ! Symmetry factor
     $     /96.0d0              ! spin/colour average

      end subroutine ME2realqg_htl


      subroutine ME2realgq_htl(p, aopi, ME2)
      implicit none
      include 'nlegborn.h'
      integer nlegs
      parameter (nlegs=nlegreal)
c     input:
      real * 8 p(0:3,nlegs)
      real * 8 aopi
c     output:
      real * 8 ME2
c     abbreviations:
      real * 8 p_cross12(0:3,nlegs)

c     cross p1 <-> p2
      p_cross12(:,1) = p(:,2)
      p_cross12(:,2) = p(:,1)
      p_cross12(:,3) = p(:,3)
      p_cross12(:,4) = p(:,4)
      p_cross12(:,5) = p(:,5)

      call ME2realqg_htl(p_cross12, aopi, ME2)

      end subroutine ME2realgq_htl


      subroutine ME2realqqbar_htl(p, aopi, ME2)
      implicit none
      include 'nlegborn.h'
      integer nlegs
      parameter (nlegs=nlegreal)
c     input:
      real * 8 p(0:3,nlegs)
      real * 8 aopi
c     output:
      real * 8 ME2
c     abbreviations:
      real * 8 p_cross23(0:3,nlegs)

c     cross p2 <-> -p3
      p_cross23(:,1) = p(:,1)
      p_cross23(:,2) = -p(:,3)
      p_cross23(:,3) = -p(:,2)
      p_cross23(:,4) = p(:,4)
      p_cross23(:,5) = p(:,5)

      call ME2realqg_htl(p_cross23, aopi, ME2)

      ME2 = -ME2*8.d0/3.d0      ! Correct spin average, minus sign from crossing fermion line

      end subroutine ME2realqqbar_htl
