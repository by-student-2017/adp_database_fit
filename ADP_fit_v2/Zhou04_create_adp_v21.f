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
      common /pass1/ re(16),p0(16),p1(16),h(16),
     *   pk1(16),pk2(16),pk3(16),pk4(16),
     *   pk5(16),pk6(16),pk7(16),pk8(16),
     *   rhoe(16),
     *   Fi0(16),Fi1(16),Fi2(16),Fi3(16),Fi4(16),Fi5(16),Fi6(16),
     *   Fm2(16),Fm3(16),Fm4(16),Fm5(16),Fm6(16),
     *   Fn0(16),Fn1(16),Fn2(16),Fn3(16),
     *   rhoin(16),rhoout(16),
     *   rhol(16),rhoh(16),
     *   urhoe(16),
     *   ui0(16),ui1(16),ui2(16),ui3(16),ui4(16),ui5(16),ui6(16),
     *   um2(16),um3(16),um4(16),um5(16),um6(16),
     *   un0(16),un1(16),un2(16),un3(16),
     *   urhoin(16),urhoout(16),
     *   urhol(16),urhoh(16),
     *   wrhoe(16),
     *   wi0(16),wi1(16),wi2(16),wi3(16),wi4(16),wi5(16),wi6(16),
     *   wm2(16),wm3(16),wm4(16),wm5(16),wm6(16),
     *   wn0(16),wn1(16),wn2(16),wn3(16),
     *   wrhoin(16),wrhoout(16),
     *   wrhol(16),wrhoh(16),
     *   rcut(16)
      common /pass2/ amass(16),Fr(5000,16),rhor(5000,16),
     *   z2r(5000,16,16),u(5000,16,16),w(5000,16,16),
     *   blat(16),drho,dr,rc,outfile,outelem
      common /pass3/ ielement(16),ntypes,nrho,nr
      ntypes=0
10    continue
      atomtype='none'
      read(5,funccard)
      if (atomtype .eq. 'none') goto 1200
      open(unit=10,file='ADP_code_v21',form='FORMATTED',status='OLD')
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
         read(10,*) p0(ntypes)
         read(10,*) rhoe(ntypes)
         read(10,*) pk1(ntypes)
         read(10,*) pk2(ntypes)
         read(10,*) pk3(ntypes)
         read(10,*) pk4(ntypes)
         read(10,*) pk5(ntypes)
         read(10,*) pk6(ntypes)
         read(10,*) pk7(ntypes)
         read(10,*) pk8(ntypes)
         read(10,*) Fi0(ntypes)
         read(10,*) Fi1(ntypes)
         read(10,*) Fi2(ntypes)
         read(10,*) Fi3(ntypes)
         read(10,*) Fi4(ntypes)
         read(10,*) Fi5(ntypes)
         read(10,*) Fi6(ntypes)
         read(10,*) Fm2(ntypes)
         read(10,*) Fm3(ntypes)
         read(10,*) Fm4(ntypes)
         read(10,*) Fm5(ntypes)
         read(10,*) Fm6(ntypes)
         read(10,*) Fn0(ntypes)
         read(10,*) Fn1(ntypes)
         read(10,*) Fn2(ntypes)
         read(10,*) Fn3(ntypes)
         read(10,*) ielement(ntypes)
         read(10,*) amass(ntypes)
         read(10,*) p1(ntypes)
         read(10,*) h(ntypes)
         read(10,*) rhol(ntypes)
         read(10,*) rhoh(ntypes)
         blat(ntypes)=sqrt(2.0d0)*re(ntypes)
         rhoin(ntypes)=rhol(ntypes)*rhoe(ntypes)
         rhoout(ntypes)=rhoh(ntypes)*rhoe(ntypes)
         rcut(ntypes)=sqrt(10.0d0)/2.0d0*blat(ntypes)
c u function
         read(10,*) urhoe(ntypes)
         read(10,*) ui0(ntypes)
         read(10,*) ui1(ntypes)
         read(10,*) ui2(ntypes)
         read(10,*) ui3(ntypes)
         read(10,*) ui4(ntypes)
         read(10,*) ui5(ntypes)
         read(10,*) ui6(ntypes)
         read(10,*) um2(ntypes)
         read(10,*) um3(ntypes)
         read(10,*) um4(ntypes)
         read(10,*) um5(ntypes)
         read(10,*) um6(ntypes)
         read(10,*) un0(ntypes)
         read(10,*) un1(ntypes)
         read(10,*) un2(ntypes)
         read(10,*) un3(ntypes)
         read(10,*) urhol(ntypes)
         read(10,*) urhoh(ntypes)
         urhoin(ntypes)=urhol(ntypes)*urhoe(ntypes)
         urhoout(ntypes)=urhoh(ntypes)*urhoe(ntypes)
c         urhos(ntypes)=urhoout(ntypes)
c w function
         read(10,*) wrhoe(ntypes)
         read(10,*) wi0(ntypes)
         read(10,*) wi1(ntypes)
         read(10,*) wi2(ntypes)
         read(10,*) wi3(ntypes)
         read(10,*) wi4(ntypes)
         read(10,*) wi5(ntypes)
         read(10,*) wi6(ntypes)
         read(10,*) wm2(ntypes)
         read(10,*) wm3(ntypes)
         read(10,*) wm4(ntypes)
         read(10,*) wm5(ntypes)
         read(10,*) wm6(ntypes)
         read(10,*) wn0(ntypes)
         read(10,*) wn1(ntypes)
         read(10,*) wn2(ntypes)
         read(10,*) wn3(ntypes)
         read(10,*) wrhol(ntypes)
         read(10,*) wrhoh(ntypes)
         wrhoin(ntypes)=wrhol(ntypes)*wrhoe(ntypes)
         wrhoout(ntypes)=wrhoh(ntypes)*wrhoe(ntypes)
c         wrhos(ntypes)=wrhoout(ntypes)
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
            u(i,i1,i2)=0.0d0
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
         if ( (rhoF .ge. rhoin(it)) .and.
     *        (rhoF .lt. rhoin(it)+drho) ) then
           call embed0(it,rhoF,emb)
           embb11 = emb
           embb12 = Fr(i,it)
         endif
         if ( (rhoF .ge. rhoout(it)) .and. 
     *        (rhoF .lt. rhoout(it)+drho) ) then
           call embed0(it,rhoF,emb)
           embb21 = emb
           embb22 = Fr(i,it)
         endif
20    continue
c
      do 21 i=1,nr
         r=(i-1)*dr
c u bounary
         if ( (r .ge. urhoin(it)) .and.
     *        (r .lt. urhoin(it)+dr) ) then
           call uembed0(it,r,emb)
           embb11u = emb
           embb12u = u(i,it,it)
         endif
         if ( (r .ge. urhoout(it)) .and.
     *        (r .lt. urhoout(it)+dr) ) then
           call uembed0(it,r,emb)
           embb21u = emb
           embb22u = u(i,it,it)
         endif
c w boundary
         if ( (r .ge. wrhoin(it)) .and.
     *        (r .lt. wrhoin(it)+dr) ) then
           call wembed0(it,r,emb)
           embb11w = emb
           embb12w = w(i,it,it)
         endif
         if ( (r .ge. wrhoout(it)) .and.
     *        (r .lt. wrhoout(it)+dr) ) then
           call wembed0(it,r,emb)
           embb21w = emb
           embb22w = w(i+1,it,it)
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
      common /pass1/ re(16),p0(16),p1(16),h(16),
     *   pk1(16),pk2(16),pk3(16),pk4(16),
     *   pk5(16),pk6(16),pk7(16),pk8(16),
     *   rhoe(16),
     *   Fi0(16),Fi1(16),Fi2(16),Fi3(16),Fi4(16),Fi5(16),Fi6(16),
     *   Fm2(16),Fm3(16),Fm4(16),Fm5(16),Fm6(16),
     *   Fn0(16),Fn1(16),Fn2(16),Fn3(16),
     *   rhoin(16),rhoout(16),
     *   rhol(16),rhoh(16),
     *   urhoe(16),
     *   ui0(16),ui1(16),ui2(16),ui3(16),ui4(16),ui5(16),ui6(16),
     *   um2(16),um3(16),um4(16),um5(16),um6(16),
     *   un0(16),un1(16),un2(16),un3(16),
     *   urhoin(16),urhoout(16),
     *   urhol(16),urhoh(16),
     *   wrhoe(16),
     *   wi0(16),wi1(16),wi2(16),wi3(16),wi4(16),wi5(16),wi6(16),
     *   wm2(16),wm3(16),wm4(16),wm5(16),wm6(16),
     *   wn0(16),wn1(16),wn2(16),wn3(16),
     *   wrhoin(16),wrhoout(16),
     *   wrhol(16),wrhoh(16),
     *   rcut(16)
c Zhou04
c      f=fe(it)*exp(-beta1(it)*(r/re(it)-1.0d0))
c      f=f/(1.0d0+(r/re(it)-ramda1(it))**20)
c fe,beta1,re,ramda1
c re = rc
c fe = p0
c beta1 = p1
c ramda1 = h
      x=(r-re(it))/h(it)
      if (x .ge. 0.0) then
        f=0.0d0
      else
        f=p0(it)*exp(-p1(it)*r)*(x**4.0d0/(1.0d0+x**4.0d0))
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine calculates the pair potential.                  c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pair(it1,it2,r,psi)
      implicit real*8 (a-h,o-z)
      implicit integer*8 (i-m)
      common /pass1/ re(16),p0(16),p1(16),h(16),
     *   pk1(16),pk2(16),pk3(16),pk4(16),
     *   pk5(16),pk6(16),pk7(16),pk8(16),
     *   rhoe(16),
     *   Fi0(16),Fi1(16),Fi2(16),Fi3(16),Fi4(16),Fi5(16),Fi6(16),
     *   Fm2(16),Fm3(16),Fm4(16),Fm5(16),Fm6(16),
     *   Fn0(16),Fn1(16),Fn2(16),Fn3(16),
     *   rhoin(16),rhoout(16),
     *   rhol(16),rhoh(16),
     *   urhoe(16),
     *   ui0(16),ui1(16),ui2(16),ui3(16),ui4(16),ui5(16),ui6(16),
     *   um2(16),um3(16),um4(16),um5(16),um6(16),
     *   un0(16),un1(16),un2(16),un3(16),
     *   urhoin(16),urhoout(16),
     *   urhol(16),urhoh(16),
     *   wrhoe(16),
     *   wi0(16),wi1(16),wi2(16),wi3(16),wi4(16),wi5(16),wi6(16),
     *   wm2(16),wm3(16),wm4(16),wm5(16),wm6(16),
     *   wn0(16),wn1(16),wn2(16),wn3(16),
     *   wrhoin(16),wrhoout(16),
     *   wrhol(16),wrhoh(16),
     *   rcut(16)
      if (it1 .eq. it2) then
c Zhou04
c         psi1=A(it1)*exp(-alpha(it1)*(r/re(it1)-1.0d0))
c         psi1=psi1/(1.0d0+(r/re(it1)-cai(it1))**20)
c         psi2=B(it1)*exp(-beta(it1)*(r/re(it1)-1.0d0))
c         psi2=psi2/(1.0d0+(r/re(it1)-ramda(it1))**20)
c         psi=psi1-psi2
c A,alpha,re,cai,B,beta,ramda
c re = rc
c h = ramda1
c alpha = k1
c beta = k2
c A = k3
c B = k4
c cai = k5
c ramda = k6
        x=(r-re(it1))/h(it1)
        if (x .ge. 0.0d0) then
          psi=0.0d0
        else
          psi=(pk1(it1)*r**-12.0d0 + pk2(it1)*r**-6.0d0 + 
     *          pk3(it1)*r**-3.0d0 + pk4(it1)*r**-2.0d0 + 
     *          pk5(it1)*r**-1.0d0 + pk6(it1) + pk7(it1)*r + 
     *          pk8(it1)*r**6.0d0)*(x**4.0d0/(1.0d0+x**4.0d0))
        endif
      else
c Zhou04
c         psi1=A(it1)*exp(-alpha(it1)*(r/re(it1)-1.0d0))
c         psi1=psi1/(1.0d0+(r/re(it1)-cai(it1))**20)
c         psi2=B(it1)*exp(-beta(it1)*(r/re(it1)-1.0d0))
c         psi2=psi2/(1.0d0+(r/re(it1)-ramda(it1))**20)
c         psia=psi1-psi2
c         psi1=A(it2)*exp(-alpha(it2)*(r/re(it2)-1.0d0))
c         psi1=psi1/(1.0d0+(r/re(it2)-cai(it2))**20)
c         psi2=B(it2)*exp(-beta(it2)*(r/re(it2)-1.0d0))
c         psi2=psi2/(1.0d0+(r/re(it2)-ramda(it2))**20)
c         psib=psi1-psi2
c A,alpha,re,cai,B,beta,ramda
         x=(r-re(it1))/h(it1)
         if (x .ge. 0.0d0) then
           psia=0.0d0
         else
           psia=(pk1(it1)*r**-12.0d0 + pk2(it1)*r**-6.0d0 +
     *           pk3(it1)*r**-3.0d0 + pk4(it1)*r**-2.0d0 +
     *           pk5(it1)*r**-1.0d0 + pk6(it1) + pk7(it1)*r +
     *           pk8(it1)*r**6.0d0)*(x**4.0d0/(1.0d0+x**4.0d0))
         endif
         x=(r-re(it2))/h(it2)
         if (x .ge. 0.0d0) then
           psib=0.0d0
         else
           psib=(pk1(it2)*r**-12.0d0 + pk2(it2)*r**-6.0d0 +
     *           pk3(it2)*r**-3.0d0 + pk4(it2)*r**-2.0d0 +
     *           pk5(it2)*r**-1.0d0 + pk6(it2) + pk7(it2)*r +
     *           pk8(it2)*r**6.0d0)*(x**4.0d0/(1.0d0+x**4.0d0))
         endif
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
      common /pass1/ re(16),p0(16),p1(16),h(16),
     *   pk1(16),pk2(16),pk3(16),pk4(16),
     *   pk5(16),pk6(16),pk7(16),pk8(16),
     *   rhoe(16),
     *   Fi0(16),Fi1(16),Fi2(16),Fi3(16),Fi4(16),Fi5(16),Fi6(16),
     *   Fm2(16),Fm3(16),Fm4(16),Fm5(16),Fm6(16),
     *   Fn0(16),Fn1(16),Fn2(16),Fn3(16),
     *   rhoin(16),rhoout(16),
     *   rhol(16),rhoh(16),
     *   urhoe(16),
     *   ui0(16),ui1(16),ui2(16),ui3(16),ui4(16),ui5(16),ui6(16),
     *   um2(16),um3(16),um4(16),um5(16),um6(16),
     *   un0(16),un1(16),un2(16),un3(16),
     *   urhoin(16),urhoout(16),
     *   urhol(16),urhoh(16),
     *   wrhoe(16),
     *   wi0(16),wi1(16),wi2(16),wi3(16),wi4(16),wi5(16),wi6(16),
     *   wm2(16),wm3(16),wm4(16),wm5(16),wm6(16),
     *   wn0(16),wn1(16),wn2(16),wn3(16),
     *   wrhoin(16),wrhoout(16),
     *   wrhol(16),wrhoh(16),
     *   rcut(16)
      if (rho .eq. 0.0d0) then
        emb = 0.0d0
      else if (rho .lt. rhoin(it)) then
         emb=Fi0(it)+
     *       Fi1(it)*(rho/rhoin(it)-1.0d0)+
     *       Fi2(it)*(rho/rhoin(it)-1.0d0)**2+
     *       Fi3(it)*(rho/rhoin(it)-1.0d0)**3+
     *       Fi4(it)*(rho/rhoin(it)-1.0d0)**4+
     *       Fi5(it)*(rho/rhoin(it)-1.0d0)**5+
     *       Fi6(it)*(rho/rhoin(it)-1.0d0)**6
      else if (rho .lt. rhoout(it)) then
         emb=Fi0(it)+
     *       Fi1(it)*(rho/rhoin(it)-1.0d0)+
     *       Fm2(it)*(rho/rhoin(it)-1.0d0)**2+
     *       Fm3(it)*(rho/rhoin(it)-1.0d0)**3+
     *       Fm4(it)*(rho/rhoin(it)-1.0d0)**4+
     *       Fm5(it)*(rho/rhoin(it)-1.0d0)**5+
     *       Fm6(it)*(rho/rhoin(it)-1.0d0)**6
      else
         emb=Fn0(it)+
     *       Fn1(it)*(rho/rhoout(it)-1.0d0)+
     *       Fn2(it)*(rho/rhoout(it)-1.0d0)**2+
     *       Fn3(it)*(rho/rhoout(it)-1.0d0)**3
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine calculates the embedding energy for boundary.   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine embed0(it,rho,emb)
      implicit real*8 (a-h,o-z)
      implicit integer*8 (i-m)
      common /pass1/ re(16),p0(16),p1(16),h(16),
     *   pk1(16),pk2(16),pk3(16),pk4(16),
     *   pk5(16),pk6(16),pk7(16),pk8(16),
     *   rhoe(16),
     *   Fi0(16),Fi1(16),Fi2(16),Fi3(16),Fi4(16),Fi5(16),Fi6(16),
     *   Fm2(16),Fm3(16),Fm4(16),Fm5(16),Fm6(16),
     *   Fn0(16),Fn1(16),Fn2(16),Fn3(16),
     *   rhoin(16),rhoout(16),
     *   rhol(16),rhoh(16),
     *   urhoe(16),
     *   ui0(16),ui1(16),ui2(16),ui3(16),ui4(16),ui5(16),ui6(16),
     *   um2(16),um3(16),um4(16),um5(16),um6(16),
     *   un0(16),un1(16),un2(16),un3(16),
     *   urhoin(16),urhoout(16),
     *   urhol(16),urhoh(16),
     *   wrhoe(16),
     *   wi0(16),wi1(16),wi2(16),wi3(16),wi4(16),wi5(16),wi6(16),
     *   wm2(16),wm3(16),wm4(16),wm5(16),wm6(16),
     *   wn0(16),wn1(16),wn2(16),wn3(16),
     *   wrhoin(16),wrhoout(16),
     *   wrhol(16),wrhoh(16),
     *   rcut(16)
      if (rho .eq. 0.0d0) then
        emb = 0.0d0
      else if (rho .lt. rhoin(it)) then
         emb=Fi0(it)+
     *       Fi1(it)*(rho/rhoin(it)-1.0d0)+
     *       Fm2(it)*(rho/rhoin(it)-1.0d0)**2+
     *       Fm3(it)*(rho/rhoin(it)-1.0d0)**3+
     *       Fm4(it)*(rho/rhoin(it)-1.0d0)**4+
     *       Fm5(it)*(rho/rhoin(it)-1.0d0)**5+
     *       Fm6(it)*(rho/rhoin(it)-1.0d0)**6
      else if (rho .lt. rhoout(it)) then
         emb=Fn0(it)+
     *       Fn1(it)*(rho/rhoout(it)-1.0d0)+
     *       Fn2(it)*(rho/rhoout(it)-1.0d0)**2+
     *       Fn3(it)*(rho/rhoout(it)-1.0d0)**3
      else
         emb=Fi0(it)+
     *       Fi1(it)*(rho/rhoin(it)-1.0d0)+
     *       Fm2(it)*(rho/rhoin(it)-1.0d0)**2+
     *       Fm3(it)*(rho/rhoin(it)-1.0d0)**3+
     *       Fm4(it)*(rho/rhoin(it)-1.0d0)**4+
     *       Fm5(it)*(rho/rhoin(it)-1.0d0)**5+
     *       Fm6(it)*(rho/rhoin(it)-1.0d0)**6
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine calculates u function.                          c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine uembed(it,r,emb)
      implicit real*8 (a-h,o-z)
      implicit integer*8 (i-m)
      common /pass1/ re(16),p0(16),p1(16),h(16),
     *   pk1(16),pk2(16),pk3(16),pk4(16),
     *   pk5(16),pk6(16),pk7(16),pk8(16),
     *   rhoe(16),
     *   Fi0(16),Fi1(16),Fi2(16),Fi3(16),Fi4(16),Fi5(16),Fi6(16),
     *   Fm2(16),Fm3(16),Fm4(16),Fm5(16),Fm6(16),
     *   Fn0(16),Fn1(16),Fn2(16),Fn3(16),
     *   rhoin(16),rhoout(16),
     *   rhol(16),rhoh(16),
     *   urhoe(16),
     *   ui0(16),ui1(16),ui2(16),ui3(16),ui4(16),ui5(16),ui6(16),
     *   um2(16),um3(16),um4(16),um5(16),um6(16),
     *   un0(16),un1(16),un2(16),un3(16),
     *   urhoin(16),urhoout(16),
     *   urhol(16),urhoh(16),
     *   wrhoe(16),
     *   wi0(16),wi1(16),wi2(16),wi3(16),wi4(16),wi5(16),wi6(16),
     *   wm2(16),wm3(16),wm4(16),wm5(16),wm6(16),
     *   wn0(16),wn1(16),wn2(16),wn3(16),
     *   wrhoin(16),wrhoout(16),
     *   wrhol(16),wrhoh(16),
     *   rcut(16)
      if (r .lt. urhoin(it)) then
         emb=ui0(it)+
     *       ui1(it)*(r/urhoin(it)-1.0d0)+
     *       ui2(it)*(r/urhoin(it)-1.0d0)**2+
     *       ui3(it)*(r/urhoin(it)-1.0d0)**3+
     *       ui4(it)*(r/urhoin(it)-1.0d0)**4+
     *       ui5(it)*(r/urhoin(it)-1.0d0)**5+
     *       ui6(it)*(r/urhoin(it)-1.0d0)**6
      else if (r .lt. urhoout(it)) then
         emb=ui0(it)+
     *       ui1(it)*(r/urhoin(it)-1.0d0)+
     *       um2(it)*(r/urhoin(it)-1.0d0)**2+
     *       um3(it)*(r/urhoin(it)-1.0d0)**3+
     *       um4(it)*(r/urhoin(it)-1.0d0)**4+
     *       um5(it)*(r/urhoin(it)-1.0d0)**5+
     *       um6(it)*(r/urhoin(it)-1.0d0)**6
      else if (r .lt. rcut(it)) then
         emb=un0(it)+
     *       un1(it)*(r/rcut(it)-1.0d0)+
     *       un2(it)*(r/rcut(it)-1.0d0)**2+
     *       un3(it)*(r/rcut(it)-1.0d0)**3
      else
         emb=0.0d0
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine calculates u function for boundary.             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine uembed0(it,r,emb)
      implicit real*8 (a-h,o-z)
      implicit integer*8 (i-m)
      common /pass1/ re(16),p0(16),p1(16),h(16),
     *   pk1(16),pk2(16),pk3(16),pk4(16),
     *   pk5(16),pk6(16),pk7(16),pk8(16),
     *   rhoe(16),
     *   Fi0(16),Fi1(16),Fi2(16),Fi3(16),Fi4(16),Fi5(16),Fi6(16),
     *   Fm2(16),Fm3(16),Fm4(16),Fm5(16),Fm6(16),
     *   Fn0(16),Fn1(16),Fn2(16),Fn3(16),
     *   rhoin(16),rhoout(16),
     *   rhol(16),rhoh(16),
     *   urhoe(16),
     *   ui0(16),ui1(16),ui2(16),ui3(16),ui4(16),ui5(16),ui6(16),
     *   um2(16),um3(16),um4(16),um5(16),um6(16),
     *   un0(16),un1(16),un2(16),un3(16),
     *   urhoin(16),urhoout(16),
     *   urhol(16),urhoh(16),
     *   wrhoe(16),
     *   wi0(16),wi1(16),wi2(16),wi3(16),wi4(16),wi5(16),wi6(16),
     *   wm2(16),wm3(16),wm4(16),wm5(16),wm6(16),
     *   wn0(16),wn1(16),wn2(16),wn3(16),
     *   wrhoin(16),wrhoout(16),
     *   wrhol(16),wrhoh(16),
     *   rcut(16)
      if (r .lt. urhoin(it)) then
         emb=ui0(it)+
     *       ui1(it)*(r/urhoin(it)-1.0d0)+
     *       um2(it)*(r/urhoin(it)-1.0d0)**2+
     *       um3(it)*(r/urhoin(it)-1.0d0)**3+
     *       um4(it)*(r/urhoin(it)-1.0d0)**4+
     *       um5(it)*(r/urhoin(it)-1.0d0)**5+
     *       um6(it)*(r/urhoin(it)-1.0d0)**6
      else if (r .lt. urhoout(it)) then
         emb=un0(it)+
     *       un1(it)*(r/rcut(it)-1.0d0)+
     *       un2(it)*(r/rcut(it)-1.0d0)**2+
     *       un3(it)*(r/rcut(it)-1.0d0)**3
      else if (r .lt. rcut(it)) then
         emb=ui0(it)+
     *       ui1(it)*(r/urhoin(it)-1.0d0)+
     *       um2(it)*(r/urhoin(it)-1.0d0)**2+
     *       um3(it)*(r/urhoin(it)-1.0d0)**3+
     *       um4(it)*(r/urhoin(it)-1.0d0)**4+
     *       um5(it)*(r/urhoin(it)-1.0d0)**5+
     *       um6(it)*(r/urhoin(it)-1.0d0)**6
      else
         emb=0.0d0
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine calculates w function.                          c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wembed(it,r,emb)
      implicit real*8 (a-h,o-z)
      implicit integer*8 (i-m)
      common /pass1/ re(16),p0(16),p1(16),h(16),
     *   pk1(16),pk2(16),pk3(16),pk4(16),
     *   pk5(16),pk6(16),pk7(16),pk8(16),
     *   rhoe(16),
     *   Fi0(16),Fi1(16),Fi2(16),Fi3(16),Fi4(16),Fi5(16),Fi6(16),
     *   Fm2(16),Fm3(16),Fm4(16),Fm5(16),Fm6(16),
     *   Fn0(16),Fn1(16),Fn2(16),Fn3(16),
     *   rhoin(16),rhoout(16),
     *   rhol(16),rhoh(16),
     *   urhoe(16),
     *   ui0(16),ui1(16),ui2(16),ui3(16),ui4(16),ui5(16),ui6(16),
     *   um2(16),um3(16),um4(16),um5(16),um6(16),
     *   un0(16),un1(16),un2(16),un3(16),
     *   urhoin(16),urhoout(16),
     *   urhol(16),urhoh(16),
     *   wrhoe(16),
     *   wi0(16),wi1(16),wi2(16),wi3(16),wi4(16),wi5(16),wi6(16),
     *   wm2(16),wm3(16),wm4(16),wm5(16),wm6(16),
     *   wn0(16),wn1(16),wn2(16),wn3(16),
     *   wrhoin(16),wrhoout(16),
     *   wrhol(16),wrhoh(16),
     *   rcut(16)
      if (r .lt. wrhoin(it)) then
         emb=wi0(it)+
     *       wi1(it)*(r/wrhoin(it)-1.0d0)+
     *       wi2(it)*(r/wrhoin(it)-1.0d0)**2+
     *       wi3(it)*(r/wrhoin(it)-1.0d0)**3+
     *       wi3(it)*(r/wrhoin(it)-1.0d0)**4+
     *       wi3(it)*(r/wrhoin(it)-1.0d0)**5+
     *       wi3(it)*(r/wrhoin(it)-1.0d0)**6
      else if (r .lt. wrhoout(it)) then
         emb=wi0(it)+
     *       wi1(it)*(r/wrhoin(it)-1.0d0)+
     *       wm2(it)*(r/wrhoin(it)-1.0d0)**2+
     *       wm3(it)*(r/wrhoin(it)-1.0d0)**3+
     *       wm3(it)*(r/wrhoin(it)-1.0d0)**4+
     *       wm3(it)*(r/wrhoin(it)-1.0d0)**5+
     *       wm3(it)*(r/wrhoin(it)-1.0d0)**6
      else if (r .lt. rcut(it)) then
         emb=wn0(it)+
     *       wn1(it)*(r/wrhoin(it)-1.0d0)+
     *       wn2(it)*(r/wrhoin(it)-1.0d0)**2+
     *       wn3(it)*(r/wrhoin(it)-1.0d0)**3
      else
         emb=0.0d0
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine calculates w function for boundary.             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wembed0(it,r,emb)
      implicit real*8 (a-h,o-z)
      implicit integer*8 (i-m)
      common /pass1/ re(16),p0(16),p1(16),h(16),
     *   pk1(16),pk2(16),pk3(16),pk4(16),
     *   pk5(16),pk6(16),pk7(16),pk8(16),
     *   rhoe(16),
     *   Fi0(16),Fi1(16),Fi2(16),Fi3(16),Fi4(16),Fi5(16),Fi6(16),
     *   Fm2(16),Fm3(16),Fm4(16),Fm5(16),Fm6(16),
     *   Fn0(16),Fn1(16),Fn2(16),Fn3(16),
     *   rhoin(16),rhoout(16),
     *   rhol(16),rhoh(16),
     *   urhoe(16),
     *   ui0(16),ui1(16),ui2(16),ui3(16),ui4(16),ui5(16),ui6(16),
     *   um2(16),um3(16),um4(16),um5(16),um6(16),
     *   un0(16),un1(16),un2(16),un3(16),
     *   urhoin(16),urhoout(16),
     *   urhol(16),urhoh(16),
     *   wrhoe(16),
     *   wi0(16),wi1(16),wi2(16),wi3(16),wi4(16),wi5(16),wi6(16),
     *   wm2(16),wm3(16),wm4(16),wm5(16),wm6(16),
     *   wn0(16),wn1(16),wn2(16),wn3(16),
     *   wrhoin(16),wrhoout(16),
     *   wrhol(16),wrhoh(16),
     *   rcut(16)
      if (r .lt. wrhoin(it)) then
         emb=wi0(it)+
     *       wi1(it)*(r/wrhoin(it)-1.0d0)+
     *       wm2(it)*(r/wrhoin(it)-1.0d0)**2+
     *       wm3(it)*(r/wrhoin(it)-1.0d0)**3+
     *       wm3(it)*(r/wrhoin(it)-1.0d0)**4+
     *       wm3(it)*(r/wrhoin(it)-1.0d0)**5+
     *       wm3(it)*(r/wrhoin(it)-1.0d0)**6
      else if (r .lt. wrhoout(it)) then
         emb=wn0(it)+
     *       wn1(it)*(r/wrhoin(it)-1.0d0)+
     *       wn2(it)*(r/wrhoin(it)-1.0d0)**2+
     *       wn3(it)*(r/wrhoin(it)-1.0d0)**3
      else if (r .lt. rcut(it)) then
         emb=wi0(it)+
     *       wi1(it)*(r/wrhoin(it)-1.0d0)+
     *       wm2(it)*(r/wrhoin(it)-1.0d0)**2+
     *       wm3(it)*(r/wrhoin(it)-1.0d0)**3+
     *       wm3(it)*(r/wrhoin(it)-1.0d0)**4+
     *       wm3(it)*(r/wrhoin(it)-1.0d0)**5+
     *       wm3(it)*(r/wrhoin(it)-1.0d0)**6
      else
         emb=0.0d0
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
      common /pass1/ re(16),p0(16),p1(16),h(16),
     *   pk1(16),pk2(16),pk3(16),pk4(16),
     *   pk5(16),pk6(16),pk7(16),pk8(16),
     *   rhoe(16),
     *   Fi0(16),Fi1(16),Fi2(16),Fi3(16),Fi4(16),Fi5(16),Fi6(16),
     *   Fm2(16),Fm3(16),Fm4(16),Fm5(16),Fm6(16),
     *   Fn0(16),Fn1(16),Fn2(16),Fn3(16),
     *   rhoin(16),rhoout(16),
     *   rhol(16),rhoh(16),
     *   urhoe(16),
     *   ui0(16),ui1(16),ui2(16),ui3(16),ui4(16),ui5(16),ui6(16),
     *   um2(16),um3(16),um4(16),um5(16),um6(16),
     *   un0(16),un1(16),un2(16),un3(16),
     *   urhoin(16),urhoout(16),
     *   urhol(16),urhoh(16),
     *   wrhoe(16),
     *   wi0(16),wi1(16),wi2(16),wi3(16),wi4(16),wi5(16),wi6(16),
     *   wm2(16),wm3(16),wm4(16),wm5(16),wm6(16),
     *   wn0(16),wn1(16),wn2(16),wn3(16),
     *   wrhoin(16),wrhoout(16),
     *   wrhol(16),wrhoh(16),
     *   rcut(16)
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
c      do 311 i=1,nr
c        do 311 i1=1,ntypes
c          do 311 i2=1,i1
c            u(i,i1,i2) = 0.0d0
c311   continue
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
