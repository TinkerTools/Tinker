c
c
c     ######################################################
c     ##  COPYRIGHT (C) 2023 by Zhi Wang & Jay W. Ponder  ##
c     ##                All Rights Reserved               ##
c     ######################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine exfield  --  external electric field energy  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "exfield" calculates the electrostatic energy due to an
c     external electric field applied to the system
c
c
      subroutine exfield (mode,exf)
      use atoms
      use charge
      use chgpot
      use energi
      use extfld
      use mpole
      use usage
      implicit none
      integer i,ii
      real*8 exf,e,f,phi
      real*8 xi,yi,zi
      real*8 ci,dix,diy,diz
      character*6 mode
c
c
c     zero out the external electric field energy
c
      exf = 0.0d0
      f = electric / dielec
c
c     calculate external field energy over partial charges
c
      if (mode .eq. 'CHARGE') then
!$OMP    PARALLEL default(private) shared(nion,iion,use,
!$OMP&    x,y,z,f,pchg,texfld,exf)
!$OMP    DO reduction(+:exf) schedule(guided)
         do ii = 1, nion
            i = iion(ii)
            if (use(i)) then
               xi = x(i)
               yi = y(i)
               zi = z(i)
               ci = pchg(i)
               phi = xi*texfld(1) + yi*texfld(2) + zi*texfld(3)
               e = -f * ci * phi
               exf = exf + e
            end if
         end do
!$OMP    END DO
!$OMP    END PARALLEL
      end if
c
c     calculate external field energy over atomic multipoles
c
      if (mode .eq. 'MPOLE') then
!$OMP    PARALLEL default(private) shared(npole,ipole,use,
!$OMP&    x,y,z,f,rpole,texfld,exf)
!$OMP    DO reduction(+:exf) schedule(guided)
         do ii = 1, npole
            i = ipole(ii)
            if (use(i)) then
               xi = x(i)
               yi = y(i)
               zi = z(i)
               ci = rpole(1,i)
               phi = xi*texfld(1) + yi*texfld(2) + zi*texfld(3)
               dix = rpole(2,i)
               diy = rpole(3,i)
               diz = rpole(4,i)
               e = -f * (ci*phi + dix*texfld(1)
     &                      + diy*texfld(2) + diz*texfld(3))
               exf = exf + e
            end if
         end do
!$OMP    END DO
!$OMP    END PARALLEL
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine exfield1  --  external field energy & gradient  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "exfield1" calculates the electrostatic energy, gradient and
c     virial due to an external electric field applied to the system
c
c
      subroutine exfield1 (mode,exf)
      use atoms
      use charge
      use chgpot
      use deriv
      use energi
      use extfld
      use mpole
      use usage
      use virial
      implicit none
      integer i,ii
      integer ix,iy,iz
      real*8 exf,e,f,phi
      real*8 xi,yi,zi
      real*8 ci,dix,diy,diz
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 frx,fry,frz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 fix(3),fiy(3)
      real*8 fiz(3),tem(3)
      character*6 mode
c
c
c     zero out the external electric field energy
c
      exf = 0.0d0
      f = electric / dielec
c
c     calculate energy and derivatives over partial charges
c
      if (mode .eq. 'CHARGE') then
!$OMP    PARALLEL default(private) shared(nion,iion,use,
!$OMP&    x,y,z,f,pchg,texfld,exf,dec,vir)
!$OMP    DO reduction(+:exf,dec,vir) schedule(guided)
         do ii = 1, nion
            i = iion(ii)
            if (use(i)) then
               xi = x(i)
               yi = y(i)
               zi = z(i)
               ci = pchg(i)
               phi = xi*texfld(1) + yi*texfld(2) + zi*texfld(3)
               e = -f * ci * phi
               exf = exf + e
c
c     gradient and virial components from charge interactions
c
               frx = -f * texfld(1) * ci
               fry = -f * texfld(2) * ci
               frz = -f * texfld(3) * ci
               dec(1,i) = dec(1,i) + frx
               dec(2,i) = dec(2,i) + fry
               dec(3,i) = dec(3,i) + frz
               vxx = xi * frx
               vyy = yi * fry
               vzz = zi * frz
               vxy = 0.5d0 * (yi*frx+xi*fry)
               vxz = 0.5d0 * (zi*frx+xi*frz)
               vyz = 0.5d0 * (zi*fry+yi*frz)
c
c     increment the total internal virial tensor components
c
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
            end if
         end do
!$OMP    END DO
!$OMP    END PARALLEL
      end if
c
c     calculate energy and derivatives over atomic multipoles
c
      if (mode .eq. 'MPOLE') then
!$OMP    PARALLEL default(private) shared(npole,ipole,use,
!$OMP&    x,y,z,xaxis,yaxis,zaxis,f,rpole,texfld,exf,dem,vir)
!$OMP    DO reduction(+:exf,dem,vir) schedule(guided)
         do ii = 1, npole
            i = ipole(ii)
            if (use(i)) then
               xi = x(i)
               yi = y(i)
               zi = z(i)
               ci = rpole(1,i)
               dix = rpole(2,i)
               diy = rpole(3,i)
               diz = rpole(4,i)
               phi = xi*texfld(1) + yi*texfld(2) + zi*texfld(3)
               e = -f * (ci*phi + dix*texfld(1)
     &                      + diy*texfld(2) + diz*texfld(3))
               exf = exf + e
c
c     gradient and virial components from dipole interactions
c
               tem(1) = f * (diy*texfld(3)-diz*texfld(2))
               tem(2) = f * (diz*texfld(1)-dix*texfld(3))
               tem(3) = f * (dix*texfld(2)-diy*texfld(1))
               call torque (i,tem,fix,fiy,fiz,dem)
               iz = zaxis(i)
               ix = xaxis(i)
               iy = abs(yaxis(i))
               if (iz .eq. 0)  iz = i
               if (ix .eq. 0)  ix = i
               if (iy .eq. 0)  iy = i
               xiz = x(iz) - x(i)
               yiz = y(iz) - y(i)
               ziz = z(iz) - z(i)
               xix = x(ix) - x(i)
               yix = y(ix) - y(i)
               zix = z(ix) - z(i)
               xiy = x(iy) - x(i)
               yiy = y(iy) - y(i)
               ziy = z(iy) - z(i)
               vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
               vxy = 0.5d0 * (yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                           + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
               vxz = 0.5d0 * (zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                           + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3))
               vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
               vyz = 0.5d0 * (zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                           + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
               vzz = zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
c
c     gradient and virial components from monopole interactions
c
               frx = -f * texfld(1) * ci
               fry = -f * texfld(2) * ci
               frz = -f * texfld(3) * ci
               dem(1,i) = dem(1,i) + frx
               dem(2,i) = dem(2,i) + fry
               dem(3,i) = dem(3,i) + frz
               vxx = vxx + xi*frx
               vyy = vyy + yi*fry
               vzz = vzz + zi*frz
               vxy = vxy + 0.5d0*(yi*frx+xi*fry)
               vxz = vxz + 0.5d0*(zi*frx+xi*frz)
               vyz = vyz + 0.5d0*(zi*fry+yi*frz)
c
c     increment the total internal virial tensor components
c
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
            end if
         end do
!$OMP    END DO
!$OMP    END PARALLEL
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine exfield3  --  electric field energy & analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "exfield3" calculates the electrostatic energy and partitions
c     the energy among the atomsdue to an external electric field
c     applied to the system
c
c
      subroutine exfield3 (mode,exf)
      use action
      use analyz
      use atoms
      use charge
      use chgpot
      use energi
      use extfld
      use mpole
      use usage
      implicit none
      integer i,ii
      real*8 exf,e,f,phi
      real*8 xi,yi,zi
      real*8 ci,dix,diy,diz
      character*6 mode
c
c
c     zero out the external electric field energy
c
      exf = 0.0d0
      f = electric / dielec
c
c     calculate energy and partitioning over partial charges
c
      if (mode .eq. 'CHARGE') then
!$OMP    PARALLEL default(private) shared(nion,iion,use,
!$OMP&    x,y,z,f,pchg,texfld,exf,nec,aec)
!$OMP    DO reduction(+:exf,nec,aec) schedule(guided)
         do ii = 1, nion
            i = iion(ii)
            if (use(i)) then
               xi = x(i)
               yi = y(i)
               zi = z(i)
               ci = pchg(i)
               phi = xi*texfld(1) + yi*texfld(2) + zi*texfld(3)
               e = -f * ci * phi
               exf = exf + e
               nec = nec + 1
               aec(i) = aec(i) + e
            end if
         end do
!$OMP    END DO
!$OMP    END PARALLEL
      end if
c
c     calculate energy and partitioning over atomic multipoles
c
      if (mode .eq. 'MPOLE') then
!$OMP    PARALLEL default(private) shared(npole,ipole,use,
!$OMP&    x,y,z,f,rpole,texfld,exf,nem,aem)
!$OMP    DO reduction(+:exf,nem,aem) schedule(guided)
         do ii = 1, npole
            i = ipole(ii)
            if (use(i)) then
               xi = x(i)
               yi = y(i)
               zi = z(i)
               ci = rpole(1,i)
               phi = xi*texfld(1) + yi*texfld(2) + zi*texfld(3)
               dix = rpole(2,i)
               diy = rpole(3,i)
               diz = rpole(4,i)
               e = -f * (ci*phi + dix*texfld(1)
     &                      + diy*texfld(2) + diz*texfld(3))
               exf = exf + e
               nem = nem + 1
               aem(i) = aem(i) + e
            end if
         end do
!$OMP    END DO
!$OMP    END PARALLEL
      end if
      return
      end
