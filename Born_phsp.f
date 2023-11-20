      subroutine born_phsp(xborn)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      real * 8 xborn(ndiminteg-3)
      real * 8 xjac,wt,sqrts
      real * 8 expon,tau,taumin,lntaum,ymax,ycm
      real * 8 xx(2),p1(4),p2(4),p3(4),p4(4),p12(4)
      real * 8 beta,vec(3)
      integer mu,k
      logical debug
      parameter(debug=.false.)
      logical ini
      data ini/.true./
      save ini
      if(ini) then
c     set initial- and final-state masses for Born and real
         do k=1,nlegborn
            kn_masses(k)=0
         enddo
         kn_masses(nlegreal)=0
         
c     Assign the Higgs bosons mass
         kn_masses(3)=ph_Hmass
         kn_masses(4)=ph_Hmass
         write(*,*) '*******************************************'
         write(*,*) '*******************************************'
         write(*,*) '****        Higgs bosons generated     ****'
         write(*,*) '****            on mass shell          ****'
         write(*,*) '*******************************************'
         write(*,*) '*******************************************'
         ini=.false.
      endif

      sqrts = sqrt(kn_sbeams)
      xjac=1d0

      expon=1d0/3d0
      taumin = (2d0*ph_Hmass/sqrts)**2
      lntaum = log(taumin)      
      tau = exp(lntaum*(1d0-xborn(1))**expon)
      xjac = xjac*(-lntaum*tau)*expon*(1d0-xborn(1))**(expon-1)

      kn_sborn = kn_sbeams*tau

      ymax=-0.5d0*log(tau)
      ycm=(2d0*xborn(2)**3-3d0*xborn(2)**2+2d0*xborn(2))*2d0*ymax-ymax
      xjac = xjac*2d0*ymax*(6d0*xborn(2)**2-6d0*xborn(2)+2d0)
      xx(1)=sqrt(tau)*exp(ycm)
      xx(2)=tau/xx(1)

c---if x's out of normal range abort
      if   ((xx(1) .gt. 1d0)
     & .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. 1d-8)
     & .or. (xx(2) .lt. 1d-8)) then
         write(*,*) ' error in Born phase space!, x1,x2 our of range'
         write(*,*) xx(1),xx(2)
         kn_jacborn = 0
         return
      endif

      p1(4)=xx(1)*sqrts*0.5d0
      p1(1)=0d0
      p1(2)=0d0
      p1(3)=xx(1)*sqrts*0.5d0

      p2(4)=xx(2)*sqrts*0.5d0
      p2(1)=0d0
      p2(2)=0d0
      p2(3)=-xx(2)*sqrts*0.5d0

C     total incoming momentum 
      p12 = p1+p2 

      call phi3m(xborn(3),xborn(4),p12,p3,p4,ph_Hmass,ph_Hmass,wt)
      xjac=xjac*wt

      kn_pborn(0,1) = p1(4)
      kn_pborn(0,2) = p2(4)
      kn_pborn(0,3) = p3(4)
      kn_pborn(0,4) = p4(4)

      kn_pborn(1:3,1) = p1(1:3)
      kn_pborn(1:3,2) = p2(1:3)
      kn_pborn(1:3,3) = p3(1:3)
      kn_pborn(1:3,4) = p4(1:3)

      kn_jacborn = xjac

c     now boost everything BACK along z-axis 
      kn_xb1 = xx(1)
      kn_xb2 = xx(2)
      beta=(kn_xb1-kn_xb2)/(kn_xb1+kn_xb2)
      vec(1)=0
      vec(2)=0
      vec(3)=-1
      call mboost(nlegborn-2,vec,beta,kn_pborn(:,3:),
     1     kn_cmpborn(:,3:))
      do mu=0,3
         kn_cmpborn(mu,1)=sqrt(kn_xb1*kn_xb2)*kn_beams(mu,1)
         kn_cmpborn(mu,2)=sqrt(kn_xb1*kn_xb2)*kn_beams(mu,2)
      enddo
      kn_minmass=2d0*ph_Hmass
      end


      function mass(p)
      implicit none
      real * 8 p(0:3),mass
      mass = sqrt(abs(p(0)**2-p(1)**2-p(2)**2-p(3)**2))
      end


      subroutine born_suppression(fact)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      logical ini
      data ini/.true./
      real * 8 fact,pt
      real * 8 powheginput
      external powheginput
      if (ini) then
         pt = powheginput("#ptsupp")         
         if(pt.gt.0) then
            write(*,*) ' ******** WARNING: ptsupp is deprecated'
            write(*,*) ' ******** Replace it with bornsuppfact'
         else
            pt = powheginput("#bornsuppfact")
         endif
         if(pt.ge.0) then
            write(*,*) '**************************'
            write(*,*) 'No Born suppression factor'
            write(*,*) '**************************'
         endif
         ini=.false.
      endif
      fact=1d0
      end


      subroutine set_fac_ren_scales(muf,mur)
      implicit none
      include 'PhysPars.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      real * 8 muf,mur
      logical ini
      data ini/.true./
      real *8 muref
      real *8 dotp
      external dotp
      real * 8 p34(0:3)
      logical fixedscales
      save fixedscales
      real * 8 powheginput
      external powheginput
      if(ini) then
         if(powheginput('#fixedscale').ge.1) then
            fixedscales=.true.
         else
            fixedscales=.false.
         endif
         write(*,*) '****************************************'
         write(*,*) '****************************************'
         write(*,*) '** Factorization and renormalization  **'   
         if(fixedscales) then
            write(*,*) '** scales set to: 2 * m_H             **'
            write(*,*) '****************************************'
            write(*,*) '****************************************'
         else
            write(*,*) '** scales set dynamically to:  m_HH/2 **'
            write(*,*) '****************************************'
            write(*,*) '****************************************'
            fixedscales=.false.
         endif
         ini=.false.
      endif
      if (fixedscales) then
         muref=ph_Hmass+ph_Hmass
      else
         if ((flg_btildepart.eq.'b').or.(flg_btildepart.eq.'c')) then
            p34=kn_pborn(:,3)+kn_pborn(:,4)
         elseif ((flg_btildepart.eq.'r')) then
            p34=kn_preal(:,3)+kn_preal(:,4)
         end if
         muref=sqrt(dotp(p34,p34))/2d0
      endif
      muf=muref
      mur=muref
      end
