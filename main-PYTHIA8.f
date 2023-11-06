      program main_pythia8
      implicit none
      include 'LesHouches.h'
      include "pwhg_physpar.h"
      include "pwhg_st.h"
      include 'hepevt.h'
      include 'pwhg_rnd.h'
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      integer j,k,l,m,iret
      integer maxev
      integer hdecaymode
      common/mcmaxev/maxev
      logical changescalup
      character * 4 rnd_cseed
      integer py8tune,nohad,hepmc,noMPI,bstable
      common/cpy8pars/hdecaymode, py8tune,nohad,hepmc,noMPI,bstable,
     $     rnd_cseed
      real * 8 bmass
      real * 8 powheginput,scalupfac
      external powheginput
      integer iun
      common/c_unit/iun

      if(powheginput("#changescalup").eq.1) then
         changescalup=.true.
      else
         changescalup=.false.
      endif
      hdecaymode=powheginput('#hdecaymode')
      hepmc = powheginput("#py8hepmc")
      nohad = powheginput("#py8nohad")
      bstable = powheginput("#py8bstable")
      py8tune = powheginput("#py8tune")
      nompi = powheginput("#py8nompi")

      WHCPRG='PYTHIA'

      bmass = powheginput("#bmass")
      if (bmass.le.0) then
         ! same as in init_couplings
         bmass = 4.5d0
      endif

c particle masses for reshuffling
      st_nlight = 5
      physpar_ml(1) = 0.511d-3
      physpar_ml(2) = 0.1057d0
      physpar_ml(3) = 1.777d0
      physpar_mq(1) = 0.33d0     ! up
      physpar_mq(2) = 0.33d0     ! down
      physpar_mq(3) = 0.50d0     ! strange
      physpar_mq(4) = 1.50d0     ! charm
      physpar_mq(5) = bmass      ! bottom

      if(hepmc.le.0) call init_hist

      call getmaxev(maxev)

      rnd_cseed = rnd_cwhichseed

      call lhefreadhdr(iun)
      call lhefinitemasses

      call pythia_init
     
      nevhep=0

      do l=1,maxev

         call lhefreadev(iun)

c         print*, 'old ',scalup
         if(changescalup) call py_change_scalup
c         print*, 'new ',scalup

         do m=1,5
c Insist to shower this event;
            call pythia_next(iret)
            if(iret.ne.1) then
               write(*,*) ' return code ',iret
               if(m.eq.1) then
                  write(*,*) ' Pythia could not shower this event'
                  call printleshouches
               endif
               write(*,*) ' retry'
            else
               exit
            endif
         enddo

         if(iret.eq.1) then
            call pythia_to_hepevt(nmxhep,nhep,isthep,idhep,jmohep,
     1           jdahep,phep,vhep)
c--- Uncomment to see events:
c            if(nevhep.lt.6) then
c               do j=1,nhep
c                  write(*,100)j,isthep(j),idhep(j),jmohep(1,j),
c     1           jmohep(2,j),jdahep(1,j),jdahep(2,j), (phep(k,j),k=1,5)
c               enddo
c            endif
            if(hepmc.le.0) then
               call pyanal
            else
               nevhep=nevhep+1
            endif
            if(nevhep.gt.0.and.mod(nevhep,20000).eq.0) then
               write(*,*)'# of events processed=',nevhep
               write(*,*)'# of events generated=',nevhep
               call pyaend
            endif
            if(nevhep.gt.0.and.mod(nevhep,1000).eq.0) then
               write(*,*)'# of events processed=',nevhep
               write(*,*)'# of events generated=',nevhep
            endif
         endif
      enddo

      write(*,*) 'At the end NEVHEP is ',nevhep

      call pythia_stat
!:      write(*,*) 'At the end: #warnings= ',mstu(27),' #errors= ',mstu(23)
c---user's terminal calculations
      if(hepmc.le.0) call pyaend
 100  format(i4,2x,i5,2x,i5,2x,i4,1x,i4,2x,i4,1x,i4,2x,5(d10.4,1x))
      end


      subroutine pyanal
      implicit none
      include 'LesHouches.h'
      include 'hepevt.h'
c     check parameters
      logical verbose
      parameter (verbose=.false.)
      real * 8 powheginput
      external powheginput
      nevhep=nevhep+1
      if(abs(idwtup).eq.3) xwgtup=xwgtup*xsecup(1)
      call analysis(xwgtup)
      call pwhgaccumup
      end

      subroutine getmaxev(maxev)
      integer maxev
      integer iun
      common/c_unit/iun
C--- Opens input file and counts number of events, setting MAXEV;
      call opencountunit(maxev,iun)
      end

      subroutine pyaend
      character * 6 vetoname
      character * 20 pwgprefix
      character * 100 filename
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      include 'pwhg_rnd.h'
      logical analysis_jetveto
      common/canalysis_jetveto/analysis_jetveto

      if(analysis_jetveto) then
         vetoname=trim(adjustl('_veto'))
      else
         vetoname=trim(adjustl(''))
      endif
      if(rnd_cwhichseed.ne.'none') then
         filename=pwgprefix(1:lprefix)//'POWHEG+PYTHIA8-output-'
     1        //rnd_cwhichseed//vetoname
      else
         filename=pwgprefix(1:lprefix)//'POWHEG+PYTHIA8-output'
     $        //vetoname
      endif

      call pwhgsetout
      call pwhgtopout(filename)
      end

      subroutine printleshouches
c useful for debugging
      call lhefwritev(6)
      end

c...lhefeader(nlf)
c...writes event information to a les houches events file on unit nlf.
      subroutine lhefwritev(nlf)
      implicit none
      integer nlf
      include 'LesHouches.h'
      include 'pwhg_flg.h'
      integer i,j
      write(nlf,'(a)')'<event>'
      write(nlf,210) nup,idprup,xwgtup,scalup,aqedup,aqcdup
      do 200 i=1,nup
         write(nlf,220) idup(i),istup(i),mothup(1,i),
     & mothup(2,i),icolup(1,i),icolup(2,i),(pup(j,i),j=1,5),
     & vtimup(i),spinup(i)
 200  continue
      write(nlf,'(a)')'</event>'
 210  format(1p,2(1x,i8),4(1x,e12.5))
 220  format(1p,i8,5(1x,i5),5(1x,e16.9),1x,e12.5,1x,e10.3)
      end

      subroutine py_change_scalup
      implicit none
      include 'LesHouches.h'
c      include 'hepevt.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      real * 8 ptmin,ptmin2,pcm(0:3,maxnup),beta,vec(3),ptkj2
      integer k,mu,j
      logical ini
      save ini
      data ini/.true./
      real * 8 dotp,powheginput,maxscalup
      external dotp,powheginput
      save maxscalup
      integer npart
      npart=0
      do k=3,nup
c     only light partons
         if (idup(k).eq.21.or.abs(idup(k)).le.5) then
            npart=npart+1
            do mu=1,3
               pcm(mu,npart)=pup(mu,k)
            enddo
            pcm(0,npart)=pup(4,k)
         endif
      enddo
c     compute min pt of light partons with respect to the incoming beam
      ptmin2=1d30
      do k=1,npart
         ptmin2=min(ptmin2,pcm(1,k)**2+pcm(2,k)**2)
      enddo

c     compute pt's of the final state partons with respect to each other
      beta=-(pup(3,1)+pup(3,2))/(pup(4,1)+pup(4,2))
      vec(1)=0
      vec(2)=0
      vec(3)=1
c     go in the CM frame
      call mboost(npart,vec,beta,pcm,pcm)
      do k=1,npart-1
         do j=k+1,npart
            ptkj2 = 2*dotp(pcm(0,k),pcm(0,j))*
     $           pcm(0,k)*pcm(0,j)/(pcm(0,k)+pcm(0,j))**2
            ptmin2=min(ptmin2,ptkj2)
         enddo
      enddo
      ptmin=sqrt(ptmin2)
      if(scalup.gt.ptmin) then
         if (ini) then
            write(*,*) '*************************************'
            write(*,*) 'scalup set to the min pt in the event'
            maxscalup=powheginput("#maxscalup")
            if(maxscalup.gt.0) then
               write(*,*) 'and limited below ',maxscalup,' GeV'
            endif
            write(*,*) '*************************************'
            ini=.false.
         endif
         if(maxscalup.gt.0) then
            scalup = min(ptmin,maxscalup)
         else
            scalup=ptmin
         endif
      endif
      end
