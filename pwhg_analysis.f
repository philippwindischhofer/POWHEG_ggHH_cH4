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
      integer i
      real *8 dy,dpt,dpta,dr,dra,dm
      real *8 powheginput
      external powheginput
      integer nptmin
      parameter (nptmin=3)
      character * 4 cptmin(nptmin)
      real * 8 ptminarr(nptmin)
      data cptmin/  '-000',  '-020',  '-080'/
      data ptminarr/   0d0,    20d0,    80d0/
      common/infohist/ptminarr,cptmin
      save /infohist/

      call inihists

      dy  = 0.2d0
      dpt = 40d0
      dpta= 20d0
      dm  = 40d0
      dr  = 0.2d0
      dra = 0.4d0

      do i=1,nptmin
      call bookupeqbins('sigtot'//cptmin(i),1d0,0.5d0,1.5d0)
      call bookupeqbins('Njet'//cptmin(i),1d0,-0.5d0,3.5d0)
      call bookupeqbins('Nevents'//cptmin(i),2.5d0,0,2d2)
      call bookupeqbins('HH-y'//cptmin(i),dy,-4d0,4d0)
      call bookupeqbins('HH-y-paper'//cptmin(i),0.4d0,-3.4d0,3.4d0)
      call bookupeqbins('HH-pt'//cptmin(i),dpt,0d0,700d0)
      call bookupeqbins('HH-pt-a'//cptmin(i),dpta,0d0,800d0)
      call bookupeqbins('HH-ptzoom'//cptmin(i),2d0,0d0,50d0)
      call bookupeqbins('HH-m'//cptmin(i),dm,250d0,1010d0)
      call bookupeqbins('HH-m-paper'//cptmin(i),25d0,250d0,1400d0)
      call bookupeqbins('H-y'//cptmin(i),dy,-4d0,4d0)
      call bookupeqbins('H-y-paper'//cptmin(i),0.4d0,-3.4d0,3.4d0)
      call bookupeqbins('H-pt'//cptmin(i),dpt,0d0,720d0)
      call bookupeqbins('H-pt-a'//cptmin(i),dpta,0d0,1000d0)
      call bookupeqbins('H-pt-paper'//cptmin(i),20d0,0d0,700d0)
      call bookupeqbins('H1-pt'//cptmin(i),dpt,0d0,1000d0)
      call bookupeqbins('H1-pt-a'//cptmin(i),dpta,0d0,1000d0)
      call bookupeqbins('H1-pt-paper'//cptmin(i),40d0,0d0,520d0)
      call bookupeqbins('H1-y'//cptmin(i),dy,-4d0,4d0)
      call bookupeqbins('H1-y-paper'//cptmin(i),0.4d0,-3.4d0,3.4d0)
      call bookupeqbins('H2-pt'//cptmin(i),dpt,0d0,1000d0)
      call bookupeqbins('H2-pt-a'//cptmin(i),dpta,0d0,1000d0)
      call bookupeqbins('H2-pt-paper'//cptmin(i),40d0,0d0,520d0)
      call bookupeqbins('H2-y'//cptmin(i),dy,-4d0,4d0)
      call bookupeqbins('H2-y-paper'//cptmin(i),0.4d0,-3.4d0,3.4d0)
      call bookupeqbins('H1-H2-dy'//cptmin(i),dy,-4d0,4d0)
      call bookupeqbins('H1-H2-deta'//cptmin(i),dy,-4d0,4d0)
      call bookupeqbins('H1-H2-dphi'//cptmin(i),pi/20,0d0,pi)
      call bookupeqbins('H1-H2-dr'//cptmin(i),dr,0d0,10d0)
      call bookupeqbins('H1-H2-dr-a'//cptmin(i),dra,0d0,10d0)
      call bookupeqbins('j1-pt'//cptmin(i),dpt,0d0,1000d0)
      call bookupeqbins('j1-pt-a'//cptmin(i),dpta,0d0,1000d0)
      call bookupeqbins('j1-y'//cptmin(i),dy,-4d0,4d0)
      call bookupeqbins('j2-pt'//cptmin(i),dpt,0d0,1000d0)
      call bookupeqbins('j2-pt-a'//cptmin(i),dpta,0d0,1000d0)
      call bookupeqbins('j2-y'//cptmin(i),dy,-4d0,4d0)
      call bookupeqbins('H1-j1-dy'//cptmin(i),dy,-4d0,4d0)
      call bookupeqbins('H1-j1-deta'//cptmin(i),dy,-4d0,4d0)
      call bookupeqbins('H1-j1-dphi'//cptmin(i),pi/20,0d0,pi)
      call bookupeqbins('H1-j1-dr'//cptmin(i),dr,0d0,10d0)
      call bookupeqbins('H1-j1-dr-a'//cptmin(i),dra,0d0,10d0)
      call bookupeqbins('HH-pt-resum1'//cptmin(i),5d0,0d0,400d0)
      call bookupeqbins('HH-pt-resum2'//cptmin(i),2d0,0d0,100d0)
      call bookupeqbins('HH-pt-mHH1'//cptmin(i),dpta,0d0,1000d0)
      call bookupeqbins('HH-pt-mHH2'//cptmin(i),dpta,0d0,1000d0)
      call bookupeqbins('HH-pt-mHH3'//cptmin(i),dpta,0d0,1000d0)
      enddo

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
      real *8 y,pt,eta,m
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
     #.or.(isthep(ihep).eq.155).or.(isthep(ihep).eq.195))
     #.and.(idhep(ihep).eq.25)) then
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

c     JETS:
c     set up arrays for jet finding
      do jpart=1,maxtrack
         do mu=1,4
            ptrack(mu,jpart)=0d0
         enddo
         jetvec(jpart)=0
      enddo
      do jjet=1,maxjet
         do mu=1,4
            pjet(mu,jjet)=0d0
         enddo
      enddo
      j1=0
      ntracks=0
      njets=0
c     Loop over final state particles to find jets
      do ihep=1,nhep
         if ((isthep(ihep).eq.1).and.
     1    (((abs(idhep(ihep)).le.10).or.(abs(idhep(ihep)).ge.40))
c     exclude leptons, gauge and higgs bosons
     2    .or.(abs(idhep(ihep)).eq.21)))
c     but  include gluons
     3           then
            if(ntracks.eq.maxtrack) then
               write(*,*)
     $              'hh analysis: too many particles, increase maxtrack'
               stop
            endif
c     copy momenta to construct jets
            ntracks=ntracks+1
            do mu=1,4
               ptrack(mu,ntracks)=phep(mu,ihep)
            enddo
         endif
      enddo

      if(buildjets.and.ntracks.gt.0) then
************************************************************************
*     fastkt algorithm
************************************************************************
c     palg=1 is standard kt, -1 is antikt
c.....run the clustering
      palg  = -1d0
      R     = 0.4d0
      ptmin = 0d0
      call fastjetppgenkt(ptrack,ntracks,R,palg,ptmin,pjet,njets,jetvec)
c.....now we have the jets

      endif

      do i=1,nptmin
         numjets=0

         if(njets.gt.0) then
            do j=1,min(2,njets)
               ktj(j) = sqrt(pjet(1,j)**2 + pjet(2,j)**2 )
               if (ktj(j).gt.ptminarr(i)) then
                  numjets=numjets+1
               endif
            enddo
         endif
	 
c HACK
c	 if(numjets.eq.0) return
c HACK end

         call filld('sigtot'//cptmin(i),1d0,dsig)
         call filld('Nevents'//cptmin(i),abs(dsig),1d0)

         if(numjets.eq.0) then
            call filld('Njet'//cptmin(i),0d0,dsig)
         elseif(numjets.eq.1) then
            call filld('Njet'//cptmin(i),1d0,dsig)
         elseif(numjets.eq.2) then
            call filld('Njet'//cptmin(i),2d0,dsig)
         elseif(numjets.eq.3) then
            call filld('Njet'//cptmin(i),3d0,dsig)
         else
            write(*,*) 'Strange number of jets:',numjets
         endif

c     Average Higgs:
         call getyetaptmass(phep(:,hlike(1)),y,eta,pt,m)
         call filld('H-y'//cptmin(i), y, dsig)
         call filld('H-pt'//cptmin(i), pt, dsig)
         call filld('H-pt-a'//cptmin(i), pt, dsig)
         call filld('H-y-paper'//cptmin(i), y, dsig)
         call filld('H-pt-paper'//cptmin(i), pt, dsig)


c     Leading pt Higgs:
         call getyetaptmass(ph1,y,eta,pt,m)
         call filld('H1-y'//cptmin(i), y, dsig)
         call filld('H1-pt'//cptmin(i), pt, dsig)
         call filld('H1-pt-a'//cptmin(i), pt, dsig)
         call filld('H1-y-paper'//cptmin(i), y, dsig)
         call filld('H1-pt-paper'//cptmin(i), pt, dsig)

c     Subleading pt Higgs:
         call getyetaptmass(ph2,y,eta,pt,m)
         call filld('H2-y'//cptmin(i), y, dsig)
         call filld('H2-pt'//cptmin(i), pt, dsig)
         call filld('H2-pt-a'//cptmin(i), pt, dsig)
         call filld('H2-y-paper'//cptmin(i), y, dsig)
         call filld('H2-pt-paper'//cptmin(i), pt, dsig)


c     Double Higgs system:
         phh = ph1 + ph2
         call getyetaptmass(phh,y,eta,pt,m)
         call filld('HH-y'//cptmin(i), y, dsig)
         call filld('HH-pt'//cptmin(i), pt, dsig)
         call filld('HH-pt-a'//cptmin(i), pt, dsig)
         call filld('HH-ptzoom'//cptmin(i), pt, dsig)
         call filld('HH-m'//cptmin(i), m, dsig)
         call filld('HH-y-paper'//cptmin(i), y, dsig)
         call filld('HH-m-paper'//cptmin(i), m, dsig)

         if (m.le.500d0) then
            call filld('HH-pt-mHH1'//cptmin(i), pt, dsig)
         else if ((m.gt.500d0).and.(m.lt.800d0)) then
            call filld('HH-pt-mHH2'//cptmin(i), pt, dsig)
         else
            call filld('HH-pt-mHH3'//cptmin(i), pt, dsig)
         endif

         if ((m.gt.300d0).and.(m.lt.500d0)) then
            call filld('HH-pt-resum1'//cptmin(i), pt, dsig)
            call filld('HH-pt-resum2'//cptmin(i), pt, dsig)
         endif

         call deltaplot(ph1,ph2,dsig,'H1-H2',cptmin(i))

         if(numjets.gt.0) then
            call getyetaptmass(pjet(:,1),y,eta,pt,m)
            call filld('j1-y'//cptmin(i), y, dsig)
            call filld('j1-pt'//cptmin(i), pt, dsig)
            call filld('j1-pt-a'//cptmin(i), pt, dsig)
            call deltaplot(ph1,pjet(:,1),dsig,'H1-j1',cptmin(i))
         endif

         if(numjets.gt.1) then
            call getyetaptmass(pjet(:,2),y,eta,pt,m)
            call filld('j2-y'//cptmin(i), y, dsig)
            call filld('j2-pt'//cptmin(i), pt, dsig)
            call filld('j2-pt-a'//cptmin(i), pt, dsig)
         endif

      enddo
      end

      subroutine deltaplot(p1,p2,dsig,prefix,postfix)
      implicit none
      real * 8 p1(4),p2(4),dsig(7)
      character *(*) prefix,postfix
      real * 8 dy,deta,delphi,dr
      call getdydetadphidr(p1,p2,dy,deta,delphi,dr)
      call filld(prefix//'-dy'//postfix,dy,dsig)
      call filld(prefix//'-deta'//postfix,deta,dsig)
      call filld(prefix//'-dphi'//postfix,delphi,dsig)
      call filld(prefix//'-dr'//postfix,dr,dsig)
      call filld(prefix//'-dr-a'//postfix,dr,dsig)
      end

      subroutine getyetaptmass(p,y,eta,pt,mass)
      implicit none
      real * 8 p(4),y,eta,pt,mass
      call pwhg_getrapidity(p,y)
      pt=sqrt(p(1)**2+p(2)**2)
      call pwhg_getpseudorapidity(p,eta)
      call pwhg_getinvmass(p,mass)
      end

      subroutine getdydetadphidr(p1,p2,dy,deta,dphi,dr)
      implicit none
      include 'pwhg_math.h'
      real * 8 p1(*),p2(*),dy,deta,dphi,dr
      real * 8 y1,eta1,pt1,mass1,phi1
      real * 8 y2,eta2,pt2,mass2,phi2
      call getyetaptmass(p1,y1,eta1,pt1,mass1)
      call getyetaptmass(p2,y2,eta2,pt2,mass2)
      dy=y1-y2
      deta=eta1-eta2
      phi1=atan2(p1(1),p1(2))
      phi2=atan2(p2(1),p2(2))
      dphi=abs(phi1-phi2)
      dphi=min(dphi,2d0*pi-dphi)
      dr=sqrt(deta**2+dphi**2)
      end

      subroutine pwhg_getinvmass(p,m)
      implicit none
      real * 8 p(4),m
      real * 8 arg
      arg = p(4)**2-p(1)**2-p(2)**2-p(3)**2
c     SIGN(A,B) returns the value of A with the sign of B
      m=sign(1.d0,arg)*sqrt(abs(arg))
      end

      subroutine pwhg_getrapidity(p,y)
      implicit none
      real * 8 p(4),y
      real * 8 x
      if (p(4).le.0d0) then
c     should NEVER enter here
         write(*,*) 'WARNING!! pwhg_getrapidity called with energy <= 0'
         y=sign(1.d0,p(3))*1.d8
      else
         x=p(3)/p(4)
         if (x.ge.1d0.or.x.le.-1d0) then
            y = sign(1d0,x)*1d8
         else
            y=0.5*log((1+x)/(1-x))
         endif
      endif
      end

      subroutine pwhg_getpseudorapidity(p,eta)
      implicit none
      real * 8 p(4),eta
      real * 8 mod, costh
      mod = sqrt(p(1)**2+p(2)**2+p(3)**2)
      if (mod.eq.0d0) then
c     should NEVER enter here
         write(*,*)
     $        'WARNING!! pwhg_getpseudorapidity called with vector 0'
         eta=sign(1.d0,p(3))*1.d8
      else
         costh = p(3)/mod
         if (costh.ge.1d0.or.costh.le.-1d0) then
            eta = sign(1d0,costh)*1d8
         else
            eta=0.5*log((1+costh)/(1-costh))
         endif
      endif
      end
