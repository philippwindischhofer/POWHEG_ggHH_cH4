c     The next subroutines, open some histograms and prepare them
c     to receive data
c     You can substitute these  with your favourite ones
c     init   :  opens the histograms
c     topout :  closes them
c     pwhgfill  :  fills the histograms with data

      subroutine init_hist
      implicit none
      include 'LesHouches.h'
      include 'pwhg_bookhist-multi-new.h'
      include 'PhysPars.h'
      include 'pwhg_math.h'

      call inihists

      call bookupeqbins('sigtot',1d0,0d0,1d0)
      call bookupeqbins('mHH',25d0,0d0,2500d0)
      call bookupeqbins('ptH-max',10d0,0d0,1000d0)
      call bookupeqbins('ptH-min',10d0,0d0,1000d0)
      call bookupeqbins('ptHH',20d0,0d0,1000d0)
      call bookupeqbins('mHH-zoom',10d0,0d0,1000d0)
      return

      end

      subroutine analysis(dsig0)
      implicit none
      include 'hepevt.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      include 'pwhg_bookhist-multi-new.h'
      include 'pwhg_weights.h'
      include 'pwhg_rwl.h'
C     allow multiweights
c      real * 8 dsig0,dsig(1:weights_max)
      real * 8 dsig0,dsig(1:rwl_maxweights)
      logical ini
      data ini/.true./
      save ini
c     we need to tell to this analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      integer i,j
c     arrays to reconstruct jets
      integer maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=2048)
      real *8 ptrack(4,maxtrack)
      real *8 pjet(4,maxjet)
      real *8 ph1(4),ph2(4),phh(4)
      real *8 ktj(2)
      real *8 pt2h1,pt2h2
      real *8 R,ptmin,palg
      real *8 y,pt,eta,mhh
      real *8 dy,deta,dphi,dr
      integer jetvec(maxtrack),hlike(2),j1
      integer mu,jpart,jjet,found,njets,
     1     ihep,ntracks,numjets
      logical buildjets
      parameter (buildjets=.true.)
      integer nptmin
      parameter (nptmin=3)
      character * 4 cptmin(nptmin)
      real * 8 ptminarr(nptmin)
      common/infohist/ptminarr,cptmin
      save /infohist/

      if (ini) then
         write(*,*) '*****************************'
         if(whcprg.eq.'NLO'.or.whcprg.eq.'LHE') then
            write(*,*) '       NLO ANALYSIS'
         elseif(WHCPRG.eq.'LHE   ') then
            write(*,*) '       LHE ANALYSIS'
         elseif(WHCPRG.eq.'HERWIG') then
            write (*,*) '           HERWIG ANALYSIS            '
         elseif(WHCPRG.eq.'PYTHIA') then
            write (*,*) '           PYTHIA ANALYSIS            '
         endif
         write(*,*) '*****************************'

         if((weights_num.eq.0).and.(rwl_num_weights.eq.0)) then
            call setupmulti(1)
         else
            if (weights_num.gt.0) then
               call setupmulti(weights_num)
            else
               call setupmulti(rwl_num_weights)
            endif
         endif
         ini=.false.
      endif

      dsig=0
      if((weights_num.eq.0).and.(rwl_num_weights.eq.0)) then
         dsig(1)=dsig0
      else
         if (weights_num.gt.0) then
            dsig(1:weights_num)=weights_val(1:weights_num)
         else
            dsig(1:rwl_num_weights)=rwl_weights(1:rwl_num_weights)
         endif
      endif
      if(sum(abs(dsig)).eq.0) return


      found=0

c     Loop over final state particles to find Higgses
      do ihep=1,nhep
         if (((isthep(ihep).eq.1).or.(isthep(ihep).eq.2)
     $        .or.(isthep(ihep).eq.155).or.(isthep(ihep).eq.195))
     $        .and.(idhep(ihep).eq.25)) then
            found=found+1
            hlike(found)=ihep
         endif
      enddo

      if(found.lt.2) then
         write(*,*) 'ERROR: Not enough Higgs-like particles found'
         call exit(1)
      elseif(found.gt.2) then
         write(*,*) 'ERROR: more than 2 Higgs-like particles found'
         call exit(1)
      endif


c     HIGGSES:
      pt2h1 = phep(1,hlike(1))**2 + phep(2,hlike(1))**2
      pt2h2 = phep(1,hlike(2))**2 + phep(2,hlike(2))**2
      if(pt2h1.ge.pt2h2) then
         ph1 = phep(1:4,hlike(1))
         ph2 = phep(1:4,hlike(2))
      else
         ph1 = phep(1:4,hlike(2))
         ph2 = phep(1:4,hlike(1))
      endif

      !>> fill plots
      call filld('sigtot',0.5d0,dsig)
      pt=sqrt( ph1(1)**2 + ph1(2)**2 )
      call filld('ptH-max', pt, dsig)
      pt=sqrt( ph2(1)**2 + ph2(2)**2 )
      call filld('ptH-min', pt, dsig)

      !>> hh-system
      phh = ph1 + ph2
      mhh=sqrt(abs( phh(4)**2 - phh(1)**2 - phh(2)**2 - phh(3)**2 ))
      call filld('mHH',mhh,dsig)

      return
      end
