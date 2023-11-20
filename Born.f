      subroutine setborn(p,bflav,born,bornjk,bmunu)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer nlegs
      parameter (nlegs=nlegborn)
      real * 8 p(0:3,nlegs),bornjk(nlegs,nlegs)
      integer bflav(nlegs)
      real * 8 bmunu(0:3,0:3,nlegs),born
      real * 8 amp2

      amp2 = 0d0
      call compborn(p,bflav,amp2,bornjk,bmunu)
      if(amp2.eq.0d0) then
        write(*,*) 'WARNING setborn: returning 0 amplitude!'
        write(*,*) bflav(1),' ',bflav(2),'->',bflav(3),' ',bflav(4)
      endif
      born=amp2
      end


      subroutine compborn(p,bflav,born,bornjk,bmunu)
      use Born_amplitudes, only: compPolVectors
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_kn.h'
      include 'pwhg_st.h'
      include 'PhysPars.h'
      integer nlegs
      parameter (nlegs=nlegborn)
      real * 8 p(0:3,nlegs),bornjk(nlegs,nlegs)
      real * 8 pTH2,shat,that,uhat
      integer bflav(nlegs)
      real * 8 bmunu(0:3,0:3,nlegs),born,born1,born2
      real * 8 modg1, modg2
      integer mu,nu,i,j,k,l,lp
      real* 8 amp2_inf,amp2_top,colcf,prefac,check
      real * 8 gtens(0:3,0:3)
      complex * 16 mpol(-1:1,-1:1)
      complex * 16 eps(2,0:3,-1:1)
      real * 8 NA
      parameter(NA=8.0d0)
      real * 8 square_Lorentz_vector
      external square_Lorentz_vector
      real * 8 scalar_product
      external scalar_product
      data gtens/1d0, 0d0, 0d0, 0d0,
     #           0d0,-1d0, 0d0, 0d0,
     #           0d0, 0d0,-1d0, 0d0,
     #           0d0, 0d0, 0d0,-1d0/
      save gtens
      integer mtdep
      common/cmtdep/mtdep
      logical ini
      data ini/.true./
      save ini
      real *8 powheginput
      external powheginput

      if (ini) then
         ini=.false.
      endif

c******************************************************
c     Choose the type of computation:
c     mtdep = 0: heavy top limit (htl)
c     mtdep = 1: Born improved heavy top limit (full theory born)
c     mtdep = 2: full theory approx (with htl virtual)
c     mtdep = 3: full theory
c     mtdep = 4,5: checks
c******************************************************

      if(mtdep.eq.0) then
C--   Born Cross section in the large top mass limit
         call ME2born_htl(p,born)

      else if(mtdep.eq.1) then
C--   Born Cross section in the large top mass limit and rescale afterwards
         call ME2born_htl(p,born)

      else if(mtdep.eq.2) then
C--   Born Cross section with finite top mass
         call ME2born_full(p,born,mpol)

      else if(mtdep.eq.3) then
C--   Born Cross section with finite top mass
         call ME2born_full(p,born,mpol)

      else if(mtdep.eq.4) then
C--   Test full Born Glover/Bij expressions with GoSam:
         call ME2born_full(p,born1,mpol)
         call ME2born_top(p,born2)
         write(*,*) "--> Born Glover/Bij:", born1
         write(*,*) "--> Born GoSam     :", born2
         write(*,*) "--> Born: should be one:", born2/born1
         born=born1
c     pause

      else if(mtdep.eq.5) then
c   compare Born Glover/Bij expressions to values from grid
         call ME2born_full(p,born,mpol)
         call ME2born_grid(p,born1)
         write(*,*) "Born_full = ", born
         write(*,*) "Born_grid = ", born1
         write(*,*) "--> Born_full/Born_grid: ", born/born1
         born=born1
c      pause

c$$$      else if(mtdep.eq.4) then
c$$$C--   Test full GoSam expressions with large top mass limit:
c$$$         call ME2born_htl(p,born1)
c$$$         call ME2born_full(p,born,mpol)
c$$$         write(*,*) "--> Born: should be one (Full/HEFT):", born2/born1
c$$$         write(*,*) "--> Born: should be one (GS/G.vd.B):", born2/born
      else
         write(*,*) "Unknown value of 'largemtlim', abort!"
         call exit(1)
      endif

C     Born spin- and color-correlated:
      do j=1,nlegs
         do mu=0,3
            do nu=0,3
               bmunu(mu,nu,j)=0d0
            enddo
         enddo
         do k=j+1,nlegs
            bornjk(j,k)=0d0
            bornjk(k,j)=bornjk(j,k)
         enddo
         if(abs(bflav(j)).le.6) then
C     Spin-correlated Born amplitudes
            if(bflav(j).eq.0) then
               if(j.eq.1) then
                  i=2
               elseif(j.eq.2) then
                  i=1
               else
                  write(*,*) 'Error in bmunu'
                  call exit(1)
               endif
               if((mtdep.eq.0).or.(mtdep.eq.1)) then
                  do mu=0,3
                     do nu=0,3
                        bmunu(mu,nu,j) = born * ((p(mu,j)*p(nu,i)+
     $                       p(nu,j)*p(mu,i))/kn_sborn-gtens(mu,nu)/2d0)
                     enddo
                  enddo
               else if(mtdep.ge.2) then
                  call compPolVectors(p,eps)
                  prefac = NA*ph_GF**2/2.0d0/256.0d0
                  do mu=0,3
                     do nu=0,3
                        do l=-1,1,2
                           do lp=-1,1,2
                              do k=-1,1,2
                                 bmunu(mu,nu,j) = bmunu(mu,nu,j) +
     $                                mpol(k,l)*conjg(mpol(k,lp))/2d0*
     $                                eps(j,mu,l)*conjg(eps(j,nu,lp))*
     $                                prefac*(st_alpha/(2d0*pi))**2
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
                  do mu=0,3
                     check=check+gtens(mu,mu)*bmunu(mu,mu,j)
                  enddo

C--   Print for checking purposes:
c                  print *, j, "check/born=", check/born
c                  print *, "bmunu",j,":"
c                  print *, bmunu(:,0,j)
c$$$                  print *, bmunu(:,1,j)
c$$$                  print *, bmunu(:,2,j)
c$$$                  print *, bmunu(:,3,j)
c$$$                  print *, "bmunu2:"
c$$$                  print *, bmunu(:,0,2)
c$$$                  print *, bmunu(:,1,2)
c$$$                  print *, bmunu(:,2,2)
c$$$                  print *, bmunu(:,3,2)
               endif
            endif

C--   Colour factors for colour-correlated Born amplitudes:
C--   Rule from 2.98 in FNO2007, leads to B_ij=Cj * B, where i#j
            do k=j+1,nlegs
               if(bflav(k).eq.0) then
                  colcf=CA
                  bornjk(j,k)=born*colcf
                  bornjk(k,j)=bornjk(j,k)
               endif
            enddo
         endif
      enddo
      end

      subroutine ME2born_htl(p, amp2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'PhysPars.h'
      include 'pwhg_math.h'
      integer nlegs
      parameter (nlegs=nlegborn)
c     input:
      real * 8 p(0:3,nlegs)
c     output:
      real * 8 amp2
c     invariants, abbreviations:
      real * 8 square_Lorentz_vector
      external square_Lorentz_vector
      real * 8 s,denH,MH2
      real * 8 NA,TR
      parameter(NA=8.0d0,TR=0.5d0)
      real * 8 modg1,modg2
      real * 8 ao2pi
c      real * 8 tmp
      real * 8 ct,chhh,ctt,cg,cgg

      MH2   = ph_Hmass2
      s     = square_Lorentz_vector(p(:,1)+p(:,2))
      denH  = 1/(s-MH2)*ph_mdlchhh
      ao2pi = st_alpha/(2d0*pi)

      ct = ph_mdlct
      chhh = ph_mdlchhh
      ctt = 2d0*ph_mdlcthh
      cg = ph_mdlcgg/8d0
      cgg = ph_mdlcgghh/4d0

c      modg1 = 4d0/9d0
c      modg1 = modg1 - 8d0/3d0*MH2*ph_mdlchhh*denH
c      modg1 = modg1 + 4d0*(MH2*ph_mdlchhh*denH)**2
c      write(*,*) 'modg1 = ', modg1
c      tmp = modg1

      modg1 = 2*(
     $      + 32*cgg**2
     $     + 16./3.*ctt*cgg
     $     + 2./9.*ctt**2
     $     - 16./3.*ct**2*cgg
     $     - 4./9.*ct**2*ctt
     $     + 2./9.*ct**4
     $     + 192*MH2*denH*cg*cgg
     $     + 16*MH2*denH*ctt*cg
     $     + 16*MH2*denH*ct*cgg
     $     + 4./3.*MH2*denH*ct*ctt
     $     - 16*MH2*denH*ct**2*cg
     $     - 4./3.*MH2*denH*ct**3
     $     + 288*MH2**2*denH**2*cg**2
     4     + 48*MH2**2*denH**2*ct*cg
     $     + 2*MH2**2*denH**2*ct**2
     $     )

c      write(*,*) 'modg1_new = ', modg1
c      write(*,*) 'Ratio (should be one) = ', tmp/modg1

      modg2 = 0

      amp2 = (modg1 + modg2) * NA*s**2
     $      *ao2pi**2*ph_GF**2
     $      /2.0d0               ! Symmetry factor
     $      /256.0d0             ! spin/colour average

      end subroutine ME2born_htl

      subroutine ME2born_full(p, amp2, mpol)
      use Born_amplitudes, only: ME2born_gbij
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'PhysPars.h'
      include 'pwhg_math.h'
      integer nlegs
      parameter (nlegs=nlegborn)
c     input:
      real * 8 p(0:3,nlegs)
c     output:
      real * 8 amp2
      complex * 16 mpol(-1:1,-1:1)
c     invariants, abbreviations:
      real * 8 MH2, MT2, ao2pi
      real * 8 NA
      parameter(NA=8.0d0)

      MH2   = ph_Hmass2
      MT2   = ph_topmass**2
      ao2pi = st_alpha/(2d0*pi)

      call ME2born_gbij(p, MH2, MT2,st_muren2, ph_mdlchhh, ph_mdlct,
     $                  ph_mdlcthh, ph_mdlcgg, ph_mdlcgghh, mpol)

c     sum over all helicity configurations:
      amp2 = mpol(1,1)*conjg(mpol(1,1)) + mpol(1,-1)*conjg(mpol(1,-1))

      amp2 = amp2 * NA
     $      *ao2pi**2*ph_GF**2
     $      /2.0d0               ! Symmetry factor
     $      /256.0d0             ! spin/colour average

      end subroutine ME2born_full


      function BI_rescaling()
      use Born_amplitudes, only: ME2born_gbij
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'PhysPars.h'
      include 'pwhg_math.h'
      integer nlegs
      parameter (nlegs=nlegborn)
      integer mtdep
      common/cmtdep/mtdep
c     output:
      real * 8 BI_rescaling
c     invariants, abbreviations:
      real * 8 square_Lorentz_vector
      external square_Lorentz_vector
      real * 8 amp2full, amp2heft
      complex * 16 mpol(-1:1,-1:1)
      real * 8 modg1,modg2
      real * 8 s, denH, MH2, MT2
      real * 8 NA
      parameter(NA=8.0d0)
      real *8 dotp
      external dotp

      if (mtdep.eq.1) then

         modg1 = 0d0
         modg2 = 0d0

         MH2   = ph_Hmass2
         MT2   = ph_topmass**2
         s     = kn_sborn
         denH  = 1/(s-MH2)

         modg1 = 4d0/9d0
         modg1 = modg1 - 8d0/3d0*MH2*ph_mdlchhh*denH
         modg1 = modg1 + 4d0*(MH2*ph_mdlchhh*denH)**2
         modg2 = 0

         call ME2born_gbij(kn_cmpborn, MH2, MT2,st_muren2,
     $                     ph_mdlchhh, ph_mdlct, ph_mdlcthh,
     $                     ph_mdlcgg, ph_mdlcgghh, mpol)

c     sum over all helicity configurations:
         amp2full = mpol(1,1)*conjg(mpol(1,1)) +
     $              mpol(1,-1)*conjg(mpol(1,-1))

         BI_rescaling = amp2full / ((modg1 + modg2) *s**2)
      else
         BI_rescaling = 1d0
      endif

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ME2born_grid(p,born)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'PhysPars.h'
      include 'pwhg_math.h'
      integer nlegs
      parameter (nlegs=nlegborn)
c     input:
      real * 8 p(0:3,nlegs)
c     output:
      real * 8 born
c     invariants, abbreviations:
      real * 8 square_Lorentz_vector
      external square_Lorentz_vector
      real * 8 mh,ao2pi
      real * 8 s,t,MH2
      real * 8 NA,nf
      parameter(NA=8.0d0)
      real * 8 grid_born

      MH2   = ph_Hmass2
      s     = square_Lorentz_vector(p(:,1)+p(:,2))
      t     = square_Lorentz_vector(p(:,1)-p(:,3)) ! Conventions: t=(p1-p3)^2
      ao2pi = st_alpha/(2d0*pi)

c  note that grid.in should be symlink to Born.grid in this case!
      born = grid_born(s,t)

c     born = ! Prefactors here

      end subroutine ME2born_grid

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Wrapper function for calling the python grid function
      function grid_born(s,t)
      implicit none
      include 'pwhg_rnd.h' ! rndiwhichseed
      real * 8 grid_born
      real * 8 s,t
      real * 8 gridvirt
c      character(len=500) :: res  ! Buffer for function result
c      character(len=500) :: arg  ! Buffer for function argument
c      character(len=16) :: pyin
c      character(len=17) :: pyout
      integer parallelstage,rndiwhichseed
      common/cpwhg_info/parallelstage,rndiwhichseed
      logical verbose

      verbose = .false.

c     Use input seed to determine which FIFOs to use, e.g. seed = 1 => FIFOs: pyInputPipe-0001, pyOutputPipe-0001
c      pyin = "pyInputPipe-"
c      pyout = "pyOutputPipe-"
c      write(pyin,'(A12,I0.4)') pyin,rndiwhichseed
c      write(pyout,'(A13,I0.4)') pyout,rndiwhichseed
      if (verbose) then
c         write(*,*) "Using FIFOs:"
c         write(*,*) pyin
c         write(*,*) pyout

c         Build input to python grid
         write(*,*) "Input to grid_born:"
         write(*,*) "s = ", s
         write(*,*) "t = ", t
      endif
c      write(arg,'(ES50.40E4,A,ES50.40E4)') s,',',t
c      if (verbose) then
c         write(*,*) "Will send the following char(len=500) to python:"
c         write(*,*) arg
c      endif

c     Send input to python script
c      open(1,file=pyin,position='asis',action='write')
c      write(1,'(A)',advance='no') arg
c      close(1)

c     Receive result from python script
c      open(2,file=pyout,position='asis',action='read')
c      read(2,'(A)') res
c      close(2)

c     Parse result of python grid

c      if (verbose) then
c         write(*,*) "Got the following character(len=500) from python:"
c         write(*,*) res
c      endif
c      read(res,*) grid_born

      grid_born = gridvirt(s,t)

      if (verbose) then
         write(*,*) "Output of grid_born:"
         write(*,*) grid_born
      endif
      end function grid_born


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ME2born_top(p, amp2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'PhysPars.h'
      include 'pwhg_math.h'
      integer nlegs
      parameter (nlegs=nlegborn)
c     input:
      real * 8 p(0:3,nlegs)
      real * 8 ao2pi
c     output:
      real * 8 amp2
c     additional variables
      integer dim_mom_array,processid,i
      parameter (dim_mom_array=50)
      real * 8 pgosam(dim_mom_array)
      real * 8 params(10),muren,res(4)
      logical print_momenta

      processid=0
      
      call gosam_momenta(p,pgosam)
      muren=sqrt(st_muren2)
      params(1)=1d0
      call OLP_EvalSubProcess(processid,pgosam,muren,params,res)
      amp2=res(3)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     GOSAM returns this result with NO gs factor ==>
C     virt_gosam ->  virt_gosam * (gs^2)^AlphasPower =
C     virt_gosam * (4*pi*st_alpha)^AlphasPower
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c         amp2=amp2 * (4d0*pi*ph_alphaem)**2
         amp2=amp2 * (4d0*pi*st_alpha)**2

c     The libraries for the  scalar integrals need an extra factor of 1/(16*pi^2)
c     Since this contribution is the square of a one-loop calculation, this factor has to be squared
     $             * 1d0/(16d0*pi**2)**2
c     There is an extra factor of 8*pi from the instrisic GoSam normalization (4*2*pi)
     $             * 8d0*pi

      processid=1
      
      call gosam_momenta(p,pgosam)
      muren=sqrt(st_muren2)
      params(1)=1d0
      call OLP_EvalSubProcess(processid,pgosam,muren,params,res)

      amp2 = amp2 + res(4) * (4d0*pi*st_alpha)**2 
     $                     * 1d0/(16d0*pi**2)
     $            + res(3) * (4d0*pi*st_alpha)**2
     $                     * 1d0/(16d0*pi**2)
c     $            + res(3) * (4d0*pi*st_alpha)**2 
c     $                     * 1d0/(16d0*pi**2)**2
c     $                     * 8d0 * pi


      end subroutine ME2born_top


      function square_Lorentz_vector(p)
      implicit none
      real * 8 square_Lorentz_vector
      real * 8 p(0:3)
      square_Lorentz_vector = p(0)**2 - p(1)**2 - p(2)**2 - p(3)**2
      return
      end function square_Lorentz_vector

      function scalar_product(p1,p2)
      implicit none
      real * 8 scalar_product
      real * 8 p1(0:3),p2(0:3)
      scalar_product = p1(0)*p2(0) - p1(1)*p2(1)
     $     - p1(2)*p2(2) - p1(3)*p2(3)
      return
      end function scalar_product

      subroutine gosam_momenta(p,pgosam)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      real * 8 p(0:3,nlegborn)
c     real * 8 pgosam(5*nlegborn)
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

      do i=1,nlegborn
         pgosam(1+5*(i-1))=p(0,i)
         pgosam(2+5*(i-1))=p(1,i)
         pgosam(3+5*(i-1))=p(2,i)
         pgosam(4+5*(i-1))=p(3,i)
c         write(*,*) i,p(0,i)**2-p(1,i)**2-p(2,i)**2-p(3,i)**2
         pgosam(5+5*(i-1))=kn_masses(i)
      enddo
      end

      subroutine borncolour_lh
c     Sets up the colour for the given flavour configuration already
c     filled in the Les Houches interface.  In case there are several
c     colour structure, one should pick one with a probability
c     proportional to the value of the corresponding cross section, for
c     the kinematics defined in the Les Houches interface
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
c colours of incoming gluons
      integer icolgi(2)
      data icolgi/ 502, 501 /
      save icolgi

c     neutral particles
      icolup(1,3)=0
      icolup(2,3)=0
      icolup(1,4)=0
      icolup(2,4)=0
c     colored particles
      icolup(1,1)=icolgi(1)
      icolup(2,1)=icolgi(2)
      icolup(1,2)=icolgi(2)
      icolup(2,2)=icolgi(1)

      end

      subroutine finalize_lh
c     Set up the resonances whose mass must be preserved
c     on the Les Houches interface.
c
c     Since the Higgses decay is always handled by the SMC
c     no resonance is present here. Higgs bosons are final state
c     particles.
c     call add_resonance(25,3,4)
      end
