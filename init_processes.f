      subroutine init_processes
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      include 'LesHouches.h'
      integer i1,i2,i3,i4,i5,k,ii(5)
      equivalence (i1,ii(1)),(i2,ii(2)),(i3,ii(3)),
     #  (i4,ii(4)),(i5,ii(5))
      logical debug
      parameter (debug=.true.)
      integer j
      logical condition
      logical flavequiv
      external flavequiv
      real * 8 powheginput
      external powheginput
      integer mtdep
      common/cmtdep/mtdep
      integer hdecaymode
      common/chdecaymode/hdecaymode
      logical ini_test
      logical btildebornon,btildevirton,btildecollon,btilderealon
      data ini_test/.true./
      save ini_test,btildebornon,btildevirton,btildecollon,btilderealon
      real * 8 cHHH,ct,ctt,cg,cgg
c     check nlegborn. This is only a sanity check while we are TESTING
c     the code and we change often from one process to the other
      if (nlegborn.ne.4) then
         write(*,*) ' ERROR: set nlegborn to the appropriate value'
         write(*,*) ' for this process in nlegborn.h'
         stop
      endif

      if(ini_test) then
         btildebornon = .not.(powheginput("#btildeborn").eq.0)
         btildevirton = .not.(powheginput("#btildevirt").eq.0)
         btildecollon = .not.(powheginput("#btildecoll").eq.0)
         btilderealon = .not.(powheginput("#btildereal").eq.0)
         ini_test = .false.
      endif

      mtdep=3

      write(*,*) "***********************************"
      if (mtdep.eq.3) then
         write(*,*) "**  Running in full theory mode  **"
      else
         write(*,*) "**  Running mode not supported!  **"
         stop
      endif
      write(*,*) "***********************************"


CCCCCCCC  Set some flags:  CCCCCCCC

c even up upper bound limits for radiation in equivalent amplitudes
      flg_evenmaxrat = .true.

c******************************************************
c     Choose the process to be implemented
c******************************************************
c     decay products of the Higgs boson
      hdecaymode=powheginput('#hdecaymode')
c     change the LHUPI id of the process according to vector boson id
c     and decay
      lprup(1)=10000+hdecaymode ! 10000+idup of Higgs decay product of the W


*********************************************************************
c     index of the first LIGHT coloured parton in the final state
      flst_lightpart=5
*********************************************************************

*********************************************************************
***********            BORN SUBPROCESSES              ***************
*********************************************************************

      i1=0
      i2=0
      i3=25
      i4=25

      flst_nborn=1
      do k=1,4
         flst_born(k,flst_nborn)=ii(k)
      enddo

      if (debug) then
         write(*,*) ' born processes',flst_nborn
         do j=1,flst_nborn
            write(*,*) (flst_born(k,j),k=1,nlegborn)
         enddo
      endif

*********************************************************************
***********            INITIALIZE VIRTUAL             ***************
*********************************************************************

      if (mtdep.eq.3.or.mtdep.eq.5) then
         cHHH = powheginput('#chhh')
         ct = 1.0d0
         ctt = 0.0d0
         cg = 0.0d0
         cgg = 0.0d0
         if( flg_bornonly .or. .not.btildevirton )then
            write(*,*)'>>> VIRTUAL not needed!'
            write(*,*)'    not loading virt-grids...'
         else
            write(*,*) 'PROC ', ct, ctt, cg, cgg
            call initgrids(cHHH,ct,ctt,cg,cgg)
         endif
      endif

*********************************************************************
***********            REAL SUBPROCESSES              ***************
*********************************************************************
      flst_nreal=0
      condition=.false.
      do i1=-5,5
         do i2=-5,5
            do i5=-5,5
               condition=(((i1+i2).eq.0).and.(i5.eq.0)).or.
     $              (((i1+i2).eq.i5).and.((i1*i2).eq.0))
               if(condition) then
                  flst_nreal=flst_nreal+1
                  if(flst_nreal.gt.maxprocreal) goto 998
                  do k=1,5
                     flst_real(k,flst_nreal)=ii(k)
                  enddo
 10               continue
               endif
            enddo
         enddo
      enddo
      if (debug) then
         write(*,*) ' real processes',flst_nreal
         do j=1,flst_nreal
            write(*,*) (flst_real(k,j),k=1,nlegreal)
         enddo
      endif
      return
 998  write(*,*) 'init_processes: increase maxprocreal'
      stop
 999  write(*,*) 'init_processes: increase maxprocborn'
      stop
      end
