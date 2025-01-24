c
c
c     ###########################################################
c     ##  COPYRIGHT (C) 2022 by Zhi Wang & Jay William Ponder  ##
c     ##                  All Rights Reserved                  ##
c     ###########################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine readdcd  --  input of DCD coordinate archive  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "readdcd" reads in a set of Cartesian coordinates from an
c     external file in CHARMM DCD binary format
c
c
      subroutine readdcd (idcd,first)
      use atoms
      use bound
      use boxes
      use files
      use inform
      use iounit
      use math
      use titles
      implicit none
      integer i,idcd
      integer blank
      integer nframe,nprev
      integer ncrdsav,nstep
      integer nvelsav,ndfree
      integer nfixat,usebox
      integer use4d,usefq
      integer merged,vcharmm
      integer ntitle
      real*4 tdelta
      real*4, allocatable :: xs(:)
      real*4, allocatable :: ys(:)
      real*4, allocatable :: zs(:)
      logical exist,opened
      logical first
      character*4 header
      character*80 info(10)
      character*240 dcdfile
c
c
c     open the input unit if it has not already been done
c
      inquire (unit=idcd,opened=opened)
      if (.not. opened) then
         dcdfile = filename(1:leng)//'.dcd'
         call version (dcdfile,'old')
         inquire (file=dcdfile,exist=exist)
         if (exist) then
            open (unit=idcd,file=dcdfile,form='unformatted',
     &               status='old')
            rewind (unit=idcd)
         else
            write (iout,10)
   10       format (/,' READDCD  --  Unable to Find the DCD',
     &                 ' Binary Archive File')
            call fatal
         end if
      end if
c
c     read header info along with title and number of atoms
c
      abort = .true.
      if (first) then
         first = .false.
         read (idcd,err=20,end=20)  header,nframe,nprev,ncrdsav,
     &                              nstep,nvelsav,blank,blank,ndfree,
     &                              nfixat,tdelta,usebox,use4d,usefq,
     &                              merged,blank,blank,blank,blank,
     &                              blank,vcharmm
         read (idcd,err=20,end=20)  ntitle,(info(i),i=1,ntitle)
         read (idcd,err=20,end=20)  n
         if (usebox .eq. 1)  use_bounds = .true.
         title(1:80) = info(1)
      end if
c
c     quit if the binary DCD file was not parsed correctly
c
      abort = .false.
   20 continue
      if (abort) then
         write (iout,30)
   30    format (/,' READDCD  --  Error Reading Header from',
     &              ' Binary DCD File')
         call fatal
      end if
c
c     read the lattice values based on header flag value;
c     using angle values is NAMD style, cosine values is CHARMM
c
      abort = .true.
      if (use_bounds) then
         call unitcell
         read (idcd,err=40,end=60)  xbox,gamma,ybox,beta,alpha,zbox
         if (abs(alpha) .le. 1.0d0)  alpha = radian * acos(alpha)
         if (abs(beta) .le. 1.0d0)  beta = radian * acos(beta)
         if (abs(gamma) .le. 1.0d0)  gamma = radian * acos(gamma)
         call lattice
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (xs(n))
      allocate (ys(n))
      allocate (zs(n))
c
c     read the atomic coordinates along each axis in turn
c
      abort = .true.
      read (idcd,err=40,end=60)  (xs(i),i=1,n)
      read (idcd,err=40,end=40)  (ys(i),i=1,n)
      read (idcd,err=40,end=40)  (zs(i),i=1,n)
c
c     quit if the binary DCD file was not parsed correctly
c
      abort = .false.
   40 continue
      if (abort) then
         write (iout,50)
   50    format (/,' READDCD  --  Error Reading Coordinates from',
     &              ' Binary DCD File')
         call fatal
      end if
c
c     copy the atomic coordinates into the current structure
c
      do i = 1, n
         x(i) = dble(xs(i))
         y(i) = dble(ys(i))
         z(i) = dble(zs(i))
      end do
c
c     perform deallocation of some local arrays
c
   60 continue
      if (allocated(xs))  deallocate (xs)
      if (allocated(ys))  deallocate (ys)
      if (allocated(zs))  deallocate (zs)
c
c     close the input unit if opened by this routine
c
      if (.not. opened)  close (unit=idcd)
      return
      end
