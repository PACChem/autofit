      program lsfit
c
      implicit double precision (a-h,o-z)
c
      include 'param.inc'
      parameter(npp=maxcoef)
      parameter(mpp=maxcoef)
      parameter(maxzone=10)
      parameter(autocmi=219474.63067d0)
      dimension vvmin(maxzone)
      dimension coef(maxcoef),sig(maxdata),mzone(maxdata),nzone(maxzone)
      dimension basis(maxterm),sc(maxzone),sx(maxzone)

      dimension vv(maxdata),rrr(maxdata,maxpair)
      dimension wzone(maxzone),wzone2(maxzone)
      dimension rcom2(maxdata),xprint(50,20)
      dimension rcom(maxdata)
      dimension x(maxatom),y(maxatom),z(maxatom)

      character*2 dum

      common/foox/rrr,nncoef

      open(56,file="coef.dat")
      cut0=500.
      cut1=100.
      cut2=25.
      epsilon=100.d0
      read(5,*)cut0,cut1,cut2,cut3
      read(5,*)epsilon
      read(5,*)sx(1),sx(2),sx(3)
      read(5,*)sc(1),sc(2),sc(3)

      open(7,file='ai.all')
      read(7,*)ndat2,natom
      write(6,*)'Reading ',ndat2,' data '
      if (ndat2.gt.maxdata) then
        write(6,*)"ndat = ",ndat2," while maxdata = ",maxdata,
     &        ". Change in LSS-PIP."
        stop
      endif

      ndat=0
      izone=0
      do i=1,maxzone
        vvmin(i)=1.d10
        nzone(i)=0
      enddo
      do i=1,ndat2
        read(7,*)
        read(7,*)iz,dum,vvx
        if (iz.eq.-1) izone=izone+1
        do j=1,natom
          read(7,*)dum,x(j),y(j),z(j)
        enddo

        if (vvx.lt.cut0.and.vvx.ne.0.) then
          if (vvx.lt.vvmin(izone)) vvmin(izone)=vvx
          ndat=ndat+1
          nzone(izone)=nzone(izone)+1
          vv(ndat)=vvx
c        if (vvx.lt.0.1) print *,ndat,vvx
          sig(ndat)=1.d0/(epsilon/(vv(ndat)+epsilon))
          mzone(ndat)=izone
          ii=0
          do j=1,natom
            do k=j+1,natom
              ii=ii+1
        rrr(ndat,ii)=dsqrt((x(j)-x(k))**2+(y(j)-y(k))**2+(z(j)-z(k))**2)
            enddo  
          enddo  
        endif
      enddo  

      print *,ndat,ndat2

      do i=1,izone
        wzone(i)=0.d0
        wzone2(i)=0.d0
      enddo
      do i=1,ndat
        sig(i)=1.d0/(epsilon/(dabs(vv(i)-sx(mzone(i)))+epsilon))
c     &      *(dble(nzone(mzone(i)))/dble(ndat))
        wzone(mzone(i))=wzone(mzone(i))+1.d0/sig(i)
        wzone2(mzone(i))=wzone2(mzone(i))+1.d0/sig(i)**2
      enddo
      do i=1,ndat
        sig(i)=sig(i)*dsqrt(wzone2(mzone(i))/wzone2(1))
        sig(i)=sig(i)/sc(mzone(i))
      enddo
      do i=1,izone
        print *," ZONE ",i," Emin ",vvmin(i)," # ",nzone(i),
     &    ndat,wzone(i),dsqrt(wzone(i))
      enddo

      write(6,*)'Using   ',ndat,' data '

      call prepot
      ncoef=nncoef
      print *,ncoef," coefficients"

      call svdfit(vv,sig,ndat,coef,ncoef,ndat,ncoef)

      err=0.d0
      err2=0.d0 
      errx=0.d0
      errx2=0.d0
      erry=0.d0 
      erry2=0.d0 
      errz=0.d0 
      errz2=0.d0 
      errb=0.d0 
      errb2=0.d0 
      nn1=0
      nn2=0
      nn3=0
      nnb=0
      wn=0.d0
      wnx=0.d0
      wny=0.d0
      wnz=0.d0
      wnb=0.d0
      vvxm=10.d0
      vvim=10.d0
      do i=1,ndat
        call funcs1(i,basis,ncoef) 
        vvx=0.d0
        do j=1,ncoef
          vvx=vvx+coef(j)*basis(j)
          if (i.eq.1) write(6,'(a,i5,a,e20.10)')
     &  '       coef(',j,') = ',coef(j)
          if (i.eq.1) write(56,'(i5,e20.10)')j,coef(j)
c          write(21,'(2i10,5e18.8)')i,j,basis(j)
        enddo
        if (i.le.200) write(6,'(i7,99e20.10)')i,vvx,vv(i),1.d0/sig(i)
c        write(22,'(i10,10f18.8)')i,(rrr(i,k),k=1,6),
c     &    vvx,vv(i),1.d0/sig(i)
        err=err+dabs(vvx-vv(i))
        err2=err2+(vvx-vv(i))**2/sig(i)**2
        wn=wn+1.d0/sig(i)**2
        if (vv(i).lt.cut1) then
          errb=errb+dabs(vvx-vv(i))
          errb2=errb2+(vvx-vv(i))**2/sig(i)**2
          wnb=wnb+1.d0/sig(i)**2
          nnb=nnb+1
        endif
        if (vv(i).lt.cut1.and.vv(i).gt.cut2) then
          errx=errx+dabs(vvx-vv(i))
          errx2=errx2+(vvx-vv(i))**2/sig(i)**2
          wnx=wnx+1.d0/sig(i)**2
          nn1=nn1+1
        endif
        if (vv(i).lt.cut2.and.vv(i).gt.cut3) then
          erry=erry+dabs(vvx-vv(i))
          erry2=erry2+(vvx-vv(i))**2/sig(i)**2
          wny=wny+1.d0/sig(i)**2
          nn2=nn2+1
        endif
        if (vv(i).lt.cut3) then
          errz=errz+dabs(vvx-vv(i))
          errz2=errz2+(vvx-vv(i))**2/sig(i)**2
          wnz=wnz+1.d0/sig(i)**2
          nn3=nn3+1
        endif
        if (vv(i).lt.vvim) then
          vvim2=vvx
          vvim=vv(i)
        endif
        if (vvx.lt.vvxm) then
          vvxm=vvx
          vvxm2=vv(i)
        endif
      enddo
      close(56)
c      err=err/dble(ndat)
c      err2=dsqrt(err2/dble(ndat))
c      errx=errx/dble(nn1)
c      errx2=dsqrt(errx2/dble(nn1))
c      erry=erry/dble(nn2)
c      erry2=dsqrt(erry2/dble(nn2))
      err=err/dble(ndat)
      err2=dsqrt(err2/wn)
      errx=errx/dble(nn1)
      errx2=dsqrt(errx2/wnx)
      erry=erry/dble(nn2)
      erry2=dsqrt(erry2/wny)
      errz=errz/dble(nn3)
      errz2=dsqrt(errz2/wnz)
      errb=errb/dble(nnb)
      errb2=dsqrt(errb2/wnb)
      print *,'ERRORS'
      print *,' < ',cut0,ndat,err,err2
      print *,' < ',cut1,nn1,errx,errx2
      print *,' < ',cut2,nn2,erry,erry2
      print *,' < ',cut3,nn3,errz,errz2
      print *,"         fit  ai"
      print *,"fit min",vvxm,vvxm2,vvxm-vvxm2
      print *,"ai  min",vvim2,vvim,vvim2-vvim

c test
      close(7)
      open(7,file='ai.mep')
      rewind(7)
      read(7,*)ndat2,natom
      write(6,*)'Reading ',ndat2,' data '
      if (ndat2.gt.maxdata) then
        write(6,*)"ndat = ",ndat2," while maxdat = ",maxdata,
     &        ". Change in FUNC."
        stop
      endif

      do i=1,ndat2
        read(7,*)natom
        read(7,*)dum,rcom2(i),vv(i)
        sig(i)=1.d0/(epsilon/(vv(i)+epsilon))
        do j=1,natom
          read(7,*)dum,x(j),y(j),z(j)
        enddo
        ii=0
        do j=1,natom
          do k=j+1,natom
            ii=ii+1
           rrr(i,ii)=dsqrt((x(j)-x(k))**2+(y(j)-y(k))**2+(z(j)-z(k))**2)
          enddo
        enddo
      enddo

      rlast=0.d0
      k=2
      j=0
      err2=0.d0
      ndat=0
      wn=0.d0
      do i=1,ndat2
        j=j+1
        call funcs1(i,basis,ncoef)
        vvx=0.d0
        do l=1,ncoef
          vvx=vvx+coef(l)*basis(l)
        enddo
        if (vv(i).lt.cut1.and.vv(i).ne.0.d0) then
          err2=err2+(vvx-vv(i))**2/sig(i)**2
          wn=wn+1.d0/sig(i)**2
          ndat=ndat+1
        endif
        print *,i,rcom2(i),vv(i),vvx
      enddo
      err2=dsqrt(err2/wn)
      print *,' test set < ',cut1,ndat,err2

c all
      close(7)
      open(7,file='ai.test')
      rewind(7)
      read(7,*)ndat2,natom
      write(6,*)'Reading ',ndat2,' data '
      if (ndat2.gt.maxdata) then
        write(6,*)"ndat = ",ndat2," while maxdata = ",maxdata,
     &        ". Change in LSS-PIP."
        stop
      endif

      ndat=0
      izone=0
      do i=1,maxzone
        vvmin(i)=1.d10
        nzone(i)=0
      enddo
      do i=1,ndat2
        read(7,*)
        read(7,*)iz,rdum,vvx
        if (iz.eq.-1) izone=izone+1
        do j=1,natom
          read(7,*)dum,x(j),y(j),z(j)
        enddo

        if (vvx.lt.cut0.and.vvx.ne.0.) then
          if (vvx.lt.vvmin(izone)) vvmin(izone)=vvx
          ndat=ndat+1
          nzone(izone)=nzone(izone)+1
          vv(ndat)=vvx
          sig(ndat)=1.d0/(epsilon/(vv(ndat)+epsilon))
          mzone(ndat)=izone
          rcom(ndat)=rdum
          ii=0
          do j=1,natom
            do k=j+1,natom
              ii=ii+1
      rrr(ndat,ii)=dsqrt((x(j)-x(k))**2+(y(j)-y(k))**2+(z(j)-z(k))**2)
            enddo  
          enddo  
        endif
      enddo  

      do i=1,izone
        wzone(i)=0.d0
        wzone2(i)=0.d0
      enddo
      do i=1,ndat
        sig(i)=1.d0/(epsilon/(dabs(vv(i)-sx(mzone(i)))+epsilon))
c     &      *(dble(nzone(mzone(i)))/dble(ndat))
        wzone(mzone(i))=wzone(mzone(i))+1.d0/sig(i)
        wzone2(mzone(i))=wzone2(mzone(i))+1.d0/sig(i)**2
      enddo
      do i=1,ndat
        sig(i)=sig(i)*dsqrt(wzone2(mzone(i))/wzone2(1))
        sig(i)=sig(i)/sc(mzone(i))
      enddo

      rlast=0.d0
      k=2
      j=0
      err3=0.d0
      ndat3=0
      wn=0.d0
      do i=1,ndat
        j=j+1
        call funcs1(i,basis,ncoef)
        vvx=0.d0
        do l=1,ncoef
           vvx=vvx+coef(l)*basis(l)
        enddo
        if (vv(i).lt.cut1.and.vv(i).ne.0.d0) then
          err3=err3+(vvx-vv(i))**2/sig(i)**2
          if (dabs(vvx-vv(i))/sig(i).gt.20.d0) 
     & print *,rcom(i),mzone(i),vv(i),vvx
          wn=wn+1.d0/sig(i)**2
          ndat3=ndat3+1
        endif
c        print *,i,rcom2(i),vv(i),vvx
      enddo
      err3=dsqrt(err3/wn)
      print *,' test set < ',cut1,ndat3,err3

      print *,"errors",errb2,err2,err3
 
      end
