c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine prtxyz  --  output of XYZ atomic coordinates  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prtxyz" writes out a set of Cartesian atomic coordinates
c     to an external file in Tinker XYZ format
c
c
      subroutine prtxyz (ixyz)
      use atomid
      use atoms
      use bound
      use boxes
      use couple
      use files
      use inform
      use titles
      implicit none
      integer i,j,k,ixyz
      integer size,crdsiz
      real*8 crdmin,crdmax
      logical opened
      character*2 atmc
      character*2 crdc
      character*2 digc
      character*25 fstr
      character*240 xyzfile
c
c
c     open the output unit if not already done
c
      inquire (unit=ixyz,opened=opened)
      if (.not. opened) then
         xyzfile = filename(1:leng)//'.xyz'
         call version (xyzfile,'new')
         open (unit=ixyz,file=xyzfile,status='new')
      end if
c
c     check for large systems needing extended formatting
c
      atmc = 'i6'
      if (n .ge. 100000)  atmc = 'i7'
      if (n .ge. 1000000)  atmc = 'i8'
      crdmin = 0.0d0
      crdmax = 0.0d0
      do i = 1, n
         crdmin = min(crdmin,x(i),y(i),z(i))
         crdmax = max(crdmax,x(i),y(i),z(i))
      end do
      crdsiz = 6
      if (crdmin .le. -1000.0d0)  crdsiz = 7
      if (crdmax .ge. 10000.0d0)  crdsiz = 7
      if (crdmin .le. -10000.0d0)  crdsiz = 8
      if (crdmax .ge. 100000.0d0)  crdsiz = 8
      crdsiz = crdsiz + max(6,digits)
      size = 0
      call numeral (crdsiz,crdc,size)
      if (digits .le. 6) then
         digc = '6 '
      else if (digits .le. 8) then
         digc = '8'
      else
         digc = '10'
      end if
c
c     write out the number of atoms and the title
c
      if (ltitle .eq. 0) then
         fstr = '('//atmc//')'
         write (ixyz,fstr(1:4))  n
      else
         fstr = '('//atmc//',2x,a)'
         write (ixyz,fstr(1:9))  n,title(1:ltitle)
      end if
c
c     write out the periodic cell lengths and angles
c
      if (use_bounds) then
         fstr = '(1x,6f'//crdc//'.'//digc//')'
         write (ixyz,fstr)  xbox,ybox,zbox,alpha,beta,gamma
      end if
c
c     write out the atomic coordinates for each atom
c
      fstr = '('//atmc//',2x,a3,3f'//crdc//
     &          '.'//digc//',i6,8'//atmc//')'
      do i = 1, n
         k = n12(i)
         if (k .eq. 0) then
            write (ixyz,fstr)  i,name(i),x(i),y(i),z(i),type(i)
         else
            write (ixyz,fstr)  i,name(i),x(i),y(i),z(i),type(i),
     &                         (i12(j,i),j=1,k)
         end if
      end do
c
c     close the output unit if opened by this routine
c
      if (.not. opened)  close (unit=ixyz)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine prtdcd  --  output of DCD atomic coordinates  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prtdcd" writes out a set of Cartesian atomic coordinates to
c     a file in CHARMM DCD binary format compatible with the VMD
c     visualization software and other packages
c
c     note the format used is based on the "dcdplugin.c" code from
c     the NAMD and VMD programs, and tutorial 4.1 from the software
c     package GENESIS: Generalized-Ensemble Simulation System
c
c     variables and parameters:
c
c     header     type of data (CORD=coordinates, VELD=velocities)
c     nframe     number of frames stored in the DCD file
c     nprev      number of previous integration steps
c     ncrdsav    frequency in steps for saving coordinate frames
c     nstep      number of integration steps in the total run
c     nvelsav    frequency of coordinate saves with velocity data
c     ndfree     number of degrees of freedom for the system
c     nfixat     number of fixed atoms for the system
c     usebox     flag for periodic boundaries (1=true, 0=false)
c     use4d      flag for 4D trajectory (1=true, 0=false)
c     usefq      flag for fluctuating charges (1=true, 0=false)
c     merged     result of merge without checks (1=true, 0=false)
c     vcharmm    version of CHARMM software for compatibility
c
c     in general a value of zero for any of the above indicates that
c     the particular feature is unused
c
c
      subroutine prtdcd (idcd,first)
      use atoms
      use bound
      use boxes
      use files
      use titles
      implicit none
      integer i,idcd
      integer zero,one
      integer nframe,nprev
      integer ncrdsav,nstep
      integer nvelsav,ndfree
      integer nfixat,usebox
      integer use4d,usefq
      integer merged,vcharmm
      integer ntitle
      real*4 tdelta
      logical opened,first
      character*4 header
      character*240 dcdfile
c
c
c     open the output unit if not already done
c
      inquire (unit=idcd,opened=opened)
      if (.not. opened) then
         dcdfile = filename(1:leng)//'.dcd'
         call version (dcdfile,'new')
         open (unit=idcd,file=dcdfile,form='unformatted',status='new')
      end if
c
c     write header info along with title and number of atoms
c
      if (first) then
         first = .false.
         zero = 0
         one = 1
         header = 'CORD'
         nframe = zero
         nprev = zero
         ncrdsav = one
         nstep = zero
         nvelsav = zero
         ndfree = zero
         nfixat = zero
         tdelta = 0.0
         usebox = zero
         if (use_bounds)  usebox = one
         use4d = zero
         usefq = zero
         merged = zero
         vcharmm = 24
         ntitle = one
         write (idcd)  header,nframe,nprev,ncrdsav,nstep,
     &                 nvelsav,zero,zero,ndfree,nfixat,
     &                 tdelta,usebox,use4d,usefq,merged,
     &                 zero,zero,zero,zero,zero,vcharmm
         write (idcd)  ntitle,title(1:80)
         write (idcd)  n
      end if
c
c     append the lattice values based on header flag value;
c     using angle values is NAMD style, cosine values is CHARMM
c
      if (use_bounds) then
c        write (idcd)  xbox,gamma_cos,ybox,beta_cos,alpha_cos,zbox
         write (idcd)  xbox,gamma,ybox,beta,alpha,zbox
      end if
c
c     append the atomic coordinates along each axis in turn
c
      write (idcd)  (real(x(i)),i=1,n)
      write (idcd)  (real(y(i)),i=1,n)
      write (idcd)  (real(z(i)),i=1,n)
c
c     close the output unit if opened by this routine
c
      if (.not. opened)  close (unit=idcd)
      return
      end
