c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2024  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine chksymm  --  check for 1D, 2D & mirror plane  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "chksymm" examines the current coordinates for linearity,
c     planarity or an internal mirror plane of symmetry
c
c
      subroutine chksymm (symmtyp)
      use atoms
      implicit none
      integer i
      real*8 eps
      real*8 xave,yave,zave
      logical xnul,ynul,znul
      character*6 symmtyp
c
c
c     copy current coordinates into a reference storage area
c
      call makeref (1)
c
c     move the atomic coordinates into the inertial frame
c
      call inertia (2)
c
c     test maximal coordinate value along each inertial axis
c
      eps = 0.001d0
      symmtyp = 'NONE'
      xnul = .true.
      ynul = .true.
      znul = .true.
      do i = 1, n
         if (abs(x(i)) .gt. eps)  xnul = .false.
         if (abs(y(i)) .gt. eps)  ynul = .false.
         if (abs(z(i)) .gt. eps)  znul = .false.
      end do
      if (xnul)  symmtyp = 'PLANAR'
      if (ynul)  symmtyp = 'PLANAR'
      if (znul)  symmtyp = 'PLANAR'
      if (xnul .and. ynul)  symmtyp = 'LINEAR'
      if (xnul .and. znul)  symmtyp = 'LINEAR'
      if (ynul .and. znul)  symmtyp = 'LINEAR'
      if (xnul .and. ynul .and. znul)  symmtyp = 'ATOMIC'
c
c     test average coordinate value along each inertial axis
c
      if (symmtyp .eq. 'NONE') then
         xave = 0.0d0
         yave = 0.0d0
         zave = 0.0d0
         do i = 1, n
            xave = xave + x(i)
            yave = yave + y(i)
            zave = zave + z(i)
         end do
         xave = abs(xave) / dble(n)
         yave = abs(yave) / dble(n)
         zave = abs(zave) / dble(n)
         if (xave .lt. eps)  symmtyp = 'MIRROR'
         if (yave .lt. eps)  symmtyp = 'MIRROR'
         if (zave .lt. eps)  symmtyp = 'MIRROR'
      end if
c
c     move original coordinates into the current structure
c
      call getref (1)
      return
      end
