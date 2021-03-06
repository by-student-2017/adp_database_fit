C author: X. W. Zhou, xzhou@sandia.gov
C updates by: Lucas Hale lucas.hale@nist.gov
C show difference of F boundary value by By Student
c      open(unit=5,file='a.i')
      call inter
c      close(5)
      call writeset
      stop
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c main subroutine.                                                c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine inter
      implicit real*8 (a-h,o-z)
      implicit integer*8 (i-m)
      character*80 atomtype,atommatch,outfile,outelem
      namelist /funccard/ atomtype
      common /pass1/ re(16),fe(16),rhoe(16),alpha(16),
     *   beta(16),beta1(16),A(16),B(16),cai(16),ramda(16),
     *   ramda1(16),Fi0(16),Fi1(16),Fi2(16),Fi3(16),
     *   Fm0(16),Fm1(16),Fm2(16),Fm3(16),Fm4(16),
     *   fnn(16),Fn(16),rhoin(16),rhoout(16),rhol(16),
     *   rhoh(16),rhos(16),
     *   urhoe(16),urhos(16),
     *   ui0(16),ui1(16),ui2(16),ui3(16),
     *   um0(16),um1(16),um21(16),um22(16),um31(16),um32(16),
     *   unn(16),un(16),urhoin(16),urhoout(16),
     *   urhol(16),urhoh(16),
     *   wi0(16),wi1(16),wi2(16),wi3(16),
     *   wrhoe(16),wrhos(16),
     *   wm0(16),wm1(16),wm21(16),wm22(16),wm31(16),wm32(16),
     *   wnn(16),wn(16),wrhoin(16),wrhoout(16),
     *   wrhol(16),wrhoh(16)
      common /pass2/ amass(16),Fr(5000,16),rhor(5000,16),
     *   z2r(5000,16,16),u(5000,16,16),w(5000,16,16),
     *   blat(16),drho,dr,rc,outfile,outelem
      common /pass3/ ielement(16),ntypes,nrho,nr
      ntypes=0
10    continue
      atomtype='none'
      read(5,funccard)
      if (atomtype .eq. 'none') goto 1200
      open(unit=10,file='ADP_code',form='FORMATTED',status='OLD')
11    read(10,9501,end=1210)atommatch
9501  format(a80)
      if (atomtype .eq. atommatch) then
         ntypes=ntypes+1
         length=len_trim(outfile)
         if (length .eq. len(outfile)) then
            outfile = atomtype
         else
            outfile = outfile(1:length)//atomtype
         endif
         length=len_trim(outelem)
         if (length .eq. len(outelem)) then
            outelem = atomtype
         else
            outelem = outelem(1:length)//' '//atomtype
         endif
         read(10,*) re(ntypes)
         read(10,*) fe(ntypes)
         read(10,*) rhoe(ntypes)
         read(10,*) rhos(ntypes)
         read(10,*) alpha(ntypes)
         read(10,*) beta(ntypes)
         read(10,*) A(ntypes)
         read(10,*) B(ntypes)
         read(10,*) cai(ntypes)
         read(10,*) ramda(ntypes)
         read(10,*) Fi0(ntypes)
         read(10,*) Fi1(ntypes)
         read(10,*) Fi2(ntypes)
         read(10,*) Fi3(ntypes)
         read(10,*) Fm0(ntypes)
         read(10,*) Fm1(ntypes)
         read(10,*) Fm2(ntypes)
         read(10,*) Fm3(ntypes)
         read(10,*) fnn(ntypes)
         read(10,*) Fn(ntypes)
         read(10,*) ielement(ntypes)
         read(10,*) amass(ntypes)
         read(10,*) Fm4(ntypes)
         read(10,*) beta1(ntypes)
         read(10,*) ramda1(ntypes)
         read(10,*) rhol(ntypes)
         read(10,*) rhoh(ntypes)
         blat(ntypes)=sqrt(2.0d0)*re(ntypes)
         rhoin(ntypes)=rhol(ntypes)*rhoe(ntypes)
         rhoout(ntypes)=rhoh(ntypes)*rhoe(ntypes)
c u function
         read(10,*) urhoe(ntypes)
c         read(10,*) urhos(ntypes)
         read(10,*) ui0(ntypes)
         read(10,*) ui1(ntypes)
         read(10,*) ui2(ntypes)
         read(10,*) ui3(ntypes)
         read(10,*) um0(ntypes)
         read(10,*) um1(ntypes)
         read(10,*) um21(ntypes)
         read(10,*) um22(ntypes)
         read(10,*) um31(ntypes)
         read(10,*) um32(ntypes)
         read(10,*) unn(ntypes)
         read(10,*) un(ntypes)
         read(10,*) urhol(ntypes)
         read(10,*) urhoh(ntypes)
         urhoin(ntypes)=urhol(ntypes)*urhoe(ntypes)
         urhoout(ntypes)=urhoh(ntypes)*urhoe(ntypes)
         urhos(ntypes)=urhoout(ntypes)
c w function
         read(10,*) wrhoe(ntypes)
c         read(10,*) wrhos(ntypes)
         read(10,*) wi0(ntypes)
         read(10,*) wi1(ntypes)
         read(10,*) wi2(ntypes)
         read(10,*) wi3(ntypes)
         read(10,*) wm0(ntypes)
         read(10,*) wm1(ntypes)
         read(10,*) wm21(ntypes)
         read(10,*) wm22(ntypes)
         read(10,*) wm31(ntypes)
         read(10,*) wm32(ntypes)
         read(10,*) wnn(ntypes)
         read(10,*) wn(ntypes)
         read(10,*) wrhol(ntypes)
         read(10,*) wrhoh(ntypes)
         wrhoin(ntypes)=wrhol(ntypes)*wrhoe(ntypes)
         wrhoout(ntypes)=wrhoh(ntypes)*wrhoe(ntypes)
         wrhos(ntypes)=wrhoout(ntypes)
c
      else
         do 1 i=1,27
1        read(10,*)vtmp
         goto 11
      endif
      close(10)
      goto 10
1210  write(6,*)'error: atom type ',atomtype,' not found'
      stop
1200  continue
      nr=2000
      nrho=2000
      alatmax=blat(1)
      rhoemax=rhoe(1)
      do 2 i=2,ntypes
         if (alatmax .lt. blat(i)) alatmax=blat(i)
         if (rhoemax .lt. rhoe(i)) rhoemax=rhoe(i)
2     continue
      rc=sqrt(10.0d0)/2.0d0*alatmax
      rst=0.5d0
      dr=rc/(nr-1.0d0)
      fmax=-1.0d0
      do 3 i1=1,ntypes
      do 3 i2=1,i1
      if ( i1 .eq. i2) then
         do 4 i=1,nr
            r=(i-1)*dr
            if (r .lt. rst) r=rst
            call prof(i1,r,fvalue)
            if (fmax .lt. fvalue) fmax=fvalue
            rhor(i,i1)=fvalue
            call pair(i1,i2,r,psi)
            z2r(i,i1,i2)=r*psi
c            r=(i-1)*dr
            call uembed(i1,r,emb)
            u(i,i1,i2)=emb
            call wembed(i1,r,emb)
            w(i,i1,i2)=emb
4        continue
      else
         do 5 i=1,nr
            r=(i-1)*dr
            if (r .lt. rst) r=rst
            call pair(i1,i2,r,psi)
            z2r(i,i1,i2)=r*psi
            z2r(i,i2,i1)=z2r(i,i1,i2)
c            r=(i-1)*dr
            call prof(i1,r,f1)
            call prof(i2,r,f2)
            u(i,i1,i2)=0.5d0*(f2/f1*u(i,i1,1)+f1/f2*u(i,i2,1))
            w(i,i1,i2)=0.5d0*(f2/f1*w(i,i1,1)+f1/f2*w(i,i2,1))
5        continue
      endif
3     continue
      rhom=fmax
      if (rhom .lt. 2.0d0*rhoemax) rhom=2.0d0*rhoemax
      if (rhom .lt. 100.0d0) rhom=100.0d0
      drho=rhom/(nrho-1.0d0)
      do 6 it=1,ntypes
      do 7 i=1,nrho
         rhoF=(i-1)*drho
         if (i .eq. 1) rhoF=0.0d0
         call embed(it,rhoF,emb)
         Fr(i,it)=emb
7     continue
      do 20 i=1,nrho
         rhoF=(i-1)*drho
         if (rhoF .lt. rhoin(it)) then
           embb11 = Fr(i,it)
           embb12 = Fr(i+1,it)
         else if (rhoF .lt. rhoout(it)) then
           embb21 = Fr(i,it) 
           embb22 = Fr(i+1,it)
         else
         endif
20    continue
c
      do 21 i=1,nr
         r=(i-1)*dr
c u bounary
         if (r .lt. urhoin(it)) then
           embb11u = u(i,it,1)
           embb12u = u(i+1,it,1)
         else if (r .lt. urhoout(it)) then
           embb21u = u(i,it,1) 
           embb22u = u(i+1,it,1)
         else
         endif
c w boundary
         if (r .lt. wrhoin(it)) then
           embb11w = w(i,it,1)
           embb12w = w(i+1,it,1)
         else if (r .lt. wrhoout(it)) then
           embb21w = w(i,it,1) 
           embb22w = w(i+1,it,1)
         else
         endif
21    continue
c
      diff = diff + abs(embb11 - embb12) + abs(embb21 - embb22)
     *            + abs(embb11u - embb12u) + abs(embb21u - embb22u)
     *            + abs(embb11w - embb12w) + abs(embb21w - embb22w)
6     continue
      open(unit=50,file='diff.dat',form='FORMATTED',status='OLD')
      write(50,*) diff
      close(50)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine calculates the electron density.                c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prof(it,r,f)
      implicit real*8 (a-h,o-z)
      implicit integer*8 (i-m)
      common /pass1/ re(16),fe(16),rhoe(16),alpha(16),
     *   beta(16),beta1(16),A(16),B(16),cai(16),ramda(16),
     *   ramda1(16),Fi0(16),Fi1(16),Fi2(16),Fi3(16),
     *   Fm0(16),Fm1(16),Fm2(16),Fm3(16),Fm4(16),
     *   fnn(16),Fn(16),rhoin(16),rhoout(16),rhol(16),
     *   rhoh(16),rhos(16),
     *   urhoe(16),urhos(16),
     *   ui0(16),ui1(16),ui2(16),ui3(16),
     *   um0(16),um1(16),um21(16),um22(16),um31(16),um32(16),
     *   unn(16),un(16),urhoin(16),urhoout(16),
     *   urhol(16),urhoh(16),
     *   wi0(16),wi1(16),wi2(16),wi3(16),
     *   wrhoe(16),wrhos(16),
     *   wm0(16),wm1(16),wm21(16),wm22(16),wm31(16),wm32(16),
     *   wnn(16),wn(16),wrhoin(16),wrhoout(16),
     *   wrhol(16),wrhoh(16)
      f=fe(it)*exp(-beta1(it)*(r/re(it)-1.0d0))
      f=f/(1.0d0+(r/re(it)-ramda1(it))**20)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine calculates the pair potential.                  c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pair(it1,it2,r,psi)
      implicit real*8 (a-h,o-z)
      implicit integer*8 (i-m)
      common /pass1/ re(16),fe(16),rhoe(16),alpha(16),
     *   beta(16),beta1(16),A(16),B(16),cai(16),ramda(16),
     *   ramda1(16),Fi0(16),Fi1(16),Fi2(16),Fi3(16),
     *   Fm0(16),Fm1(16),Fm2(16),Fm3(16),Fm4(16),
     *   fnn(16),Fn(16),rhoin(16),rhoout(16),rhol(16),
     *   rhoh(16),rhos(16),
     *   urhoe(16),urhos(16),
     *   ui0(16),ui1(16),ui2(16),ui3(16),
     *   um0(16),um1(16),um21(16),um22(16),um31(16),um32(16),
     *   unn(16),un(16),urhoin(16),urhoout(16),
     *   urhol(16),urhoh(16),
     *   wi0(16),wi1(16),wi2(16),wi3(16),
     *   wrhoe(16),wrhos(16),
     *   wm0(16),wm1(16),wm21(16),wm22(16),wm31(16),wm32(16),
     *   wnn(16),wn(16),wrhoin(16),wrhoout(16),
     *   wrhol(16),wrhoh(16)
      if (it1 .eq. it2) then
         psi1=A(it1)*exp(-alpha(it1)*(r/re(it1)-1.0d0))
         psi1=psi1/(1.0d0+(r/re(it1)-cai(it1))**20)
         psi2=B(it1)*exp(-beta(it1)*(r/re(it1)-1.0d0))
         psi2=psi2/(1.0d0+(r/re(it1)-ramda(it1))**20)
         psi=psi1-psi2
      else
         psi1=A(it1)*exp(-alpha(it1)*(r/re(it1)-1.0d0))
         psi1=psi1/(1.0d0+(r/re(it1)-cai(it1))**20)
         psi2=B(it1)*exp(-beta(it1)*(r/re(it1)-1.0d0))
         psi2=psi2/(1.0d0+(r/re(it1)-ramda(it1))**20)
         psia=psi1-psi2
         psi1=A(it2)*exp(-alpha(it2)*(r/re(it2)-1.0d0))
         psi1=psi1/(1.0d0+(r/re(it2)-cai(it2))**20)
         psi2=B(it2)*exp(-beta(it2)*(r/re(it2)-1.0d0))
         psi2=psi2/(1.0d0+(r/re(it2)-ramda(it2))**20)
         psib=psi1-psi2
         call prof(it1,r,f1)
         call prof(it2,r,f2)
         psi=0.5d0*(f2/f1*psia+f1/f2*psib)
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine calculates the embedding energy.                c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine embed(it,rho,emb)
      implicit real*8 (a-h,o-z)
      implicit integer*8 (i-m)
      common /pass1/ re(16),fe(16),rhoe(16),alpha(16),
     *   beta(16),beta1(16),A(16),B(16),cai(16),ramda(16),
     *   ramda1(16),Fi0(16),Fi1(16),Fi2(16),Fi3(16),
     *   Fm0(16),Fm1(16),Fm2(16),Fm3(16),Fm4(16),
     *   fnn(16),Fn(16),rhoin(16),rhoout(16),rhol(16),
     *   rhoh(16),rhos(16),
     *   urhoe(16),urhos(16),
     *   ui0(16),ui1(16),ui2(16),ui3(16),
     *   um0(16),um1(16),um21(16),um22(16),um31(16),um32(16),
     *   unn(16),un(16),urhoin(16),urhoout(16),
     *   urhol(16),urhoh(16),
     *   i0(16),wi1(16),wi2(16),wi3(16),
     *   wrhoe(16),wrhos(16),
     *   wm0(16),wm1(16),wm21(16),wm22(16),wm31(16),wm32(16),
     *   wnn(16),wn(16),wrhoin(16),wrhoout(16),
     *   wrhol(16),wrhoh(16)
      if (rho .lt. rhoe(it)) then
         Fm33=Fm3(it)
      else 
         Fm33=Fm4(it)
      endif
      if (rho .eq. 0.0d0) then
        emb = 0.0d0
      else if (rho .lt. rhoin(it)) then
         emb=Fi0(it)+
     *       Fi1(it)*(rho/rhoin(it)-1.0d0)+
     *       Fi2(it)*(rho/rhoin(it)-1.0d0)**2+
     *       Fi3(it)*(rho/rhoin(it)-1.0d0)**3
      else if (rho .lt. rhoout(it)) then
         emb=Fm0(it)+
     *       Fm1(it)*(rho/rhoe(it)-1.0d0)+
     *       Fm2(it)*(rho/rhoe(it)-1.0d0)**2+
     *          Fm33*(rho/rhoe(it)-1.0d0)**3
      else
         emb=Fn(it)*(1.0d0-fnn(it)*log(rho/rhos(it)))*
     *       (rho/rhos(it))**fnn(it)
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine calculates u function.                          c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine uembed(it,r,emb)
      implicit real*8 (a-h,o-z)
      implicit integer*8 (i-m)
      common /pass1/ re(16),fe(16),rhoe(16),alpha(16),
     *   beta(16),beta1(16),A(16),B(16),cai(16),ramda(16),
     *   ramda1(16),Fi0(16),Fi1(16),Fi2(16),Fi3(16),
     *   Fm0(16),Fm1(16),Fm2(16),Fm3(16),Fm4(16),
     *   fnn(16),Fn(16),rhoin(16),rhoout(16),rhol(16),
     *   rhoh(16),rhos(16),
     *   urhoe(16),urhos(16),
     *   ui0(16),ui1(16),ui2(16),ui3(16),
     *   um0(16),um1(16),um21(16),um22(16),um31(16),um32(16),
     *   unn(16),un(16),urhoin(16),urhoout(16),
     *   urhol(16),urhoh(16),
     *   wi0(16),wi1(16),wi2(16),wi3(16),
     *   wrhoe(16),wrhos(16),
     *   wm0(16),wm1(16),wm21(16),wm22(16),wm31(16),wm32(16),
     *   wnn(16),wn(16),wrhoin(16),wrhoout(16),
     *   wrhol(16),wrhoh(16)
      if (r .lt. urhoe(it)) then
         um2=um21(it)
         um3=um31(it)
      else 
         um2=um22(it)
         um3=um32(it)
      endif
      if (r .lt. urhoin(it)) then
         emb=ui0(it)+
     *       ui1(it)*(r/urhoin(it)-1.0d0)+
     *       ui2(it)*(r/urhoin(it)-1.0d0)**2+
     *       ui3(it)*(r/urhoin(it)-1.0d0)**3
      else if (r .lt. urhoout(it)) then
         emb=um0(it)+
     *       um1(it)*(r/urhoe(it)-1.0d0)+
     *           um2*(r/urhoe(it)-1.0d0)**2+
     *           um3*(r/urhoe(it)-1.0d0)**3
      else
         emb=un(it)*(1.0d0-unn(it)*log(r/urhos(it)))/
     *       (r/urhos(it))**unn(it)
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine calculates w function.                          c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wembed(it,r,emb)
      implicit real*8 (a-h,o-z)
      implicit integer*8 (i-m)
      common /pass1/ re(16),fe(16),rhoe(16),alpha(16),
     *   beta(16),beta1(16),A(16),B(16),cai(16),ramda(16),
     *   ramda1(16),Fi0(16),Fi1(16),Fi2(16),Fi3(16),
     *   Fm0(16),Fm1(16),Fm2(16),Fm3(16),Fm4(16),
     *   fnn(16),Fn(16),rhoin(16),rhoout(16),rhol(16),
     *   rhoh(16),rhos(16),
     *   urhoe(16),urhos(16),
     *   ui0(16),ui1(16),ui2(16),ui3(16),
     *   um0(16),um1(16),um21(16),um22(16),um31(16),um32(16),
     *   unn(16),un(16),urhoin(16),urhoout(16),
     *   urhol(16),urhoh(16),
     *   wi0(16),wi1(16),wi2(16),wi3(16),
     *   wrhoe(16),wrhos(16),
     *   wm0(16),wm1(16),wm21(16),wm22(16),wm31(16),wm32(16),
     *   wnn(16),wn(16),wrhoin(16),wrhoout(16),
     *   wrhol(16),wrhoh(16)
      if (r .lt. wrhoe(it)) then
         wm2=wm21(it)
         wm3=wm31(it)
      else 
         wm2=wm22(it)
         wm3=wm32(it)
      endif
      if (r .lt. wrhoin(it)) then
         emb=wi0(it)+
     *       wi1(it)*(r/wrhoin(it)-1.0d0)+
     *       wi2(it)*(r/wrhoin(it)-1.0d0)**2+
     *       wi3(it)*(r/wrhoin(it)-1.0d0)**3
      else if (r .lt. wrhoout(it)) then
         emb=wm0(it)+
     *       wm1(it)*(r/wrhoe(it)-1.0d0)+
     *           wm2*(r/wrhoe(it)-1.0d0)**2+
     *           wm3*(r/wrhoe(it)-1.0d0)**3
      else
         emb=wn(it)*(1.0d0-wnn(it)*log(r/wrhos(it)))/
     *       (r/wrhos(it))**wnn(it)
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c write out set file.                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine writeset
      implicit real*8 (a-h,o-z)
      implicit integer*8 (i-m)
      character*80 outfile,outelem
      common /pass1/ re(16),fe(16),rhoe(16),alpha(16),
     *   beta(16),beta1(16),A(16),B(16),cai(16),ramda(16),
     *   ramda1(16),Fi0(16),Fi1(16),Fi2(16),Fi3(16),
     *   Fm0(16),Fm1(16),Fm2(16),Fm3(16),Fm4(16),
     *   fnn(16),Fn(16),rhoin(16),rhoout(16),rhol(16),
     *   rhoh(16),rhos(16),
     *   urhoe(16),urhos(16),
     *   ui0(16),ui1(16),ui2(16),ui3(16),
     *   um0(16),um1(16),um21(16),um22(16),um31(16),um32(16),
     *   unn(16),un(16),urhoin(16),urhoout(16),
     *   urhol(16),urhoh(16),
     *   wi0(16),wi1(16),wi2(16),wi3(16),
     *   wrhoe(16),wrhos(16),
     *   wm0(16),wm1(16),wm21(16),wm22(16),wm31(16),wm32(16),
     *   wnn(16),wn(16),wrhoin(16),wrhoout(16),
     *   wrhol(16),wrhoh(16)
      common /pass2/ amass(16),Fr(5000,16),rhor(5000,16),
     *   z2r(5000,16,16),u(5000,16,16),w(5000,16,16),
     *   blat(16),drho,dr,rc,outfile,outelem
      common /pass3/ ielement(16),ntypes,nrho,nr
      character*80 struc
      struc='ZZZ'
      outfile = outfile(1:index(outfile,' ')-1)//'_Zhou04.adp'
      open(unit=1,file=outfile)
      write(1,*) ' DATE: 2018-03-30 ',
     *   'CONTRIBUTOR: Xiaowang Zhou xzhou@sandia.gov and ',
     *   'Lucas Hale lucas.hale@nist.gov ',
     *   'CITATION: X. W. Zhou, R. A. Johnson, ',
     *   'H. N. G. Wadley, Phys. Rev. B, 69, 144113(2004)'
      write(1,*) ' Generated from Zhou04_create_v2.f'
      write(1,*) ' Fixes precision issues with older version'
      write(1,8)ntypes,outelem
8     format(i5,' ',a24)
      write(1,9)nrho,drho,nr,dr,rc
9     format(i5,e24.16,i5,2e24.16)
      do 10 i=1,ntypes
        write(1,11)ielement(i),amass(i),blat(i),struc
        write(1,12)(Fr(j,i),j=1,nrho)
        write(1,12)(rhor(j,i),j=1,nr)
10    continue
11    format(i5,2g15.5,a8)
12    format(5e24.16)
      do 13 i1=1,ntypes
        do 13 i2=1,i1
          write(1,12)(z2r(i,i1,i2),i=1,nr)
13    continue
      do 31 i1=1,ntypes
        do 31 i2=1,i1
          write(1,12) (u(i,i1,i2),i=1,nr)
31    continue
      do 41 i1=1,ntypes
        do 41 i2=1,i1
          write(1,12) (w(i,i1,i2),i=1,nr)
41    continue
      close(1)
      return
      end
