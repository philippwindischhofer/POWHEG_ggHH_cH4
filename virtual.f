c     returns 2 Re(M_B * M_V)/(as/(2pi)),
c     where M_B is the Born amplitude and
c     M_V is the finite part of the virtual amplitude
c     The as/(2pi) factor is attached at a later point
      subroutine setvirtual(p,vflav,virtual)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      real * 8 p(0:3,nlegborn)
      integer vflav(nlegborn)
      real * 8 virtual,virtualfull,bornhtl,bornfull
      complex * 16 mpol(-1:1,-1:1)
      integer mtdep
      common/cmtdep/mtdep

c******************************************************
c     Choose the type of computation:
c     mtdep = 0: heavy top limit (htl)
c     mtdep = 1: Born improved heavy top limit (full theory born)
c     mtdep = 2: full theory approx (with htl virtual)
c     mtdep = 3: full theory
c     mtdep = 4,5: checks
c******************************************************
      virtual = 0d0

      if (mtdep.eq.0) then
         call ME2virt_htl(p, virtual)

      else if (mtdep.eq.1) then
         call ME2virt_htl(p, virtual)

      else if (mtdep.eq.2) then
         call ME2born_htl(p, bornhtl)
         call ME2born_full(p, bornfull, mpol)
         call ME2virt_htl(p, virtual)
         virtual=virtual*bornfull/bornhtl

      else if (mtdep.eq.3) then
         call ME2virt_full(p, virtual)

      else if (mtdep.eq.4) then
         call ME2virt_htl(p, virtual)
         call ME2virt_full(p, virtualfull)

         print *, "virtual HEFT=", virtual
         print *, "virtual FULL=", virtualfull
         print *, "--->Must be 1: ",virtual/virtualfull
c         pause
         virtual=virtualfull
      else
         write(*,*) "Unknown value of 'largemtlim' for virtual, abort!"
         call exit(1)
      endif
      end


      subroutine ME2virt_htl(p, virtual)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'PhysPars.h'
      include 'pwhg_math.h'
      integer nlegs
      parameter (nlegs=nlegborn)
c     input:
      real * 8 p(0:3,nlegs)
      real * 8 mh,ao2pi
c     output:
      real * 8 born,virtual,grid_virt
c     invariants, abbreviations:
      real * 8 square_Lorentz_vector
      external square_Lorentz_vector
      real * 8 s,t,u,denH,MH2
      real * 8 NA,nf,TR,st_alpha_old
      parameter(NA=8.0d0,TR=0.5d0,nf=5.d0)
      real * 8 pwhg_alphas
      external pwhg_alphas
      real * 8 muf,mur
      real * 8 tmp
      real * 8 ct,chhh,ctt,cg,cgg

      MH2   = ph_Hmass2
      s     = square_Lorentz_vector(p(:,1)+p(:,2))
      t     = square_Lorentz_vector(p(:,1)-p(:,3)) ! Conventions: t=(p1-p3)^2
      denH  = 1/(s-MH2)*ph_mdlchhh
      ao2pi = st_alpha/(2d0*pi)

c      ct = ph_mdlct
c      chhh = ph_mdlchhh
c      ctt = 2d0*ph_mdlcthh
c      cg = ph_mdlcgg/8d0
c      cgg = ph_mdlcgghh/4d0

      call ME2born_htl(p, born)

      virtual = born*(11d0+CA*pi**2-CA*log(st_muren2/s)**2)

      virtual = virtual +
     $     (ph_GF**2*ao2pi**2*NA*s**2
     $     *(-16d0/27d0+16d0/9d0*MH2*denH)
     $     /2.0d0               ! Symmetry factor
     $     /256.0d0)            ! spin/colour average


c     for test only! !!!!!!!!
c$$$      virtual = grid_virt(s,t)
c$$$
c$$$      call set_fac_ren_scales(muf,mur)
c$$$
c$$$      st_alpha_old = pwhg_alphas(mur**2,st_lambda5MSB,st_nlight)
c$$$
c$$$      virtual = virtual * st_alpha**2/st_alpha_old**2 +
c$$$     $     CA*born*(log(mur**2/s)**2 - log(st_muren2/s)**2)

c assuming Powheg will multiply by alphas/2/pi

c     for test only! !!!!!!!!

      end subroutine ME2virt_htl


      subroutine ME2virt_full(p, virtual)
c     use Born_amplitudes, only: ME2born_gbij
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'PhysPars.h'
      include 'pwhg_math.h'
      include 'pwhg_flg.h'
      integer nlegs
      parameter (nlegs=nlegborn)
c     input:
      real * 8 p(0:3,nlegs)
c     output:
      real * 8 virtual,bornfull
c     invariants, abbreviations:
      real * 8 square_Lorentz_vector
      external square_Lorentz_vector
      real * 8 mh,ao2pi,st_alpha_old
      real * 8 s,t,MH2
      real * 8 NA,NF,TR
      parameter(NA=8.d0,TR=0.5d0,NF=5.d0)
      real * 8 gridvirt
      complex * 16 mpol(-1:1,-1:1)
      real * 8 pwhg_alphas
      external pwhg_alphas
      real * 8 muf,mur
      
      MH2   = ph_Hmass2
      s     = square_Lorentz_vector(p(:,1)+p(:,2))
      t     = square_Lorentz_vector(p(:,1)-p(:,3)) ! Conventions: t=(p1-p3)^2
      ao2pi = st_alpha/(2d0*pi)

      virtual = gridvirt(s,t)

c    GH: added 23.2.17, warning: NF=5 hardcoded here!
      call ME2born_full(p, bornfull, mpol)
      call set_fac_ren_scales(muf,mur)

c    including the terms from transformation to Vfin:

      st_alpha_old = pwhg_alphas(mur**2,st_lambda5MSB,st_nlight)

      virtual = virtual * st_alpha**2/st_alpha_old**2 +
     $     CA*bornfull*(log(mur**2/s)**2 - log(st_muren2/s)**2)

      end subroutine ME2virt_full


c Wrapper function for calling the python grid function
c$$$      function grid_virt(s,t)
c$$$      implicit none
c$$$      include 'pwhg_rnd.h' ! rndiwhichseed
c$$$      real * 8 grid_virt
c$$$      real * 8 s,t
c$$$      character(len=500) :: res  ! Buffer for function result
c$$$      character(len=500) :: arg  ! Buffer for function argument
c$$$      character(len=16) :: pyin
c$$$      character(len=17) :: pyout
c$$$      integer parallelstage,rndiwhichseed
c$$$      common/cpwhg_info/parallelstage,rndiwhichseed
c$$$      logical verbose
c$$$
c$$$      verbose = .false.
c$$$
c$$$c     Use input seed to determine which FIFOs to use, e.g. seed = 1 => FIFOs: pyInputPipe-0001, pyOutputPipe-0001
c$$$      pyin = "pyInputPipe-"
c$$$      pyout = "pyOutputPipe-"
c$$$      write(pyin,'(A12,I0.4)') pyin,rndiwhichseed
c$$$      write(pyout,'(A13,I0.4)') pyout,rndiwhichseed
c$$$      if (verbose) then
c$$$         write(*,*) "Using FIFOs:"
c$$$         write(*,*) pyin
c$$$         write(*,*) pyout
c$$$
c$$$c         Build input to python grid
c$$$         write(*,*) "Input to grid_virt:"
c$$$         write(*,*) s
c$$$         write(*,*) t
c$$$      endif
c$$$      write(arg,'(ES50.40E4,A,ES50.40E4)') s,',',t
c$$$      if (verbose) then
c$$$         write(*,*) "Will send the following char(len=500) to python:"
c$$$         write(*,*) arg
c$$$      endif
c$$$
c$$$c     Send input to python script
c$$$      open(1,file=pyin,position='asis',action='write')
c$$$      write(1,'(A)',advance='no') arg
c$$$      close(1)
c$$$
c$$$c     Receive result from python script
c$$$      open(2,file=pyout,position='asis',action='read')
c$$$      read(2,'(A)') res
c$$$      close(2)
c$$$
c$$$c     Parse result of python grid
c$$$
c$$$      if (verbose) then
c$$$         write(*,*) "Got the following character(len=500) from python:"
c$$$         write(*,*) res
c$$$      endif
c$$$      read(res,*) grid_virt
c$$$
c$$$      if (verbose) then
c$$$         write(*,*) "Output of grid_virt:"
c$$$         write(*,*) grid_virt
c$$$      endif
c$$$      end function grid_virt
c$$$
c$$$
c$$$c Wrapper function for killing the python grid script
c$$$      logical function kill_python()
c$$$      implicit none
c$$$      character(len=16) :: pyin
c$$$      integer parallelstage,rndiwhichseed
c$$$      common/cpwhg_info/parallelstage,rndiwhichseed
c$$$
c$$$      ! Use input seed to determine which FIFOs to use,
c$$$      ! e.g. seed = 1 => FIFOs: pyInputPipe-0001, pyOutputPipe-0001
c$$$      pyin = "pyInputPipe-"
c$$$      write(pyin,'(A12,I0.4)') pyin,rndiwhichseed
c$$$      write(*,*) "Killing Python script with FIFO:"
c$$$      write(*,*) pyin
c$$$
c$$$      open(1,file=pyin,position='asis',action='write')
c$$$      write(1,'(A)',advance='no') 'exit'
c$$$      close(1)
c$$$      end function
