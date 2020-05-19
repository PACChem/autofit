********************************************

      subroutine prepot

      implicit double precision(a-h,o-z)
      include 'param.inc'
      dimension iagroup(maxatom),
     &  ind(maxterm,maxpair),
     &  iatom(maxperm,maxatom),
     &  idum(maxatom),nngroup(maxatom),idisc(maxterm),
     &  idum2(maxperm,maxatom),
     &  idum3(maxperm,maxatom),nperm0(maxperm),nperm1(maxperm),
     &  basis(maxterm),ibasis(maxterm),r(maxpair),
     &  rrr(maxdata,maxpair),index(maxatom,maxatom),ix(maxperm,maxpair)
      character*2 symb(maxatom),dum
      logical lreadbasis,lreaddisc
 
      common/foox/rrr,nncoef

      save npairs,nterms,ind,ibasis

      read(5,*)natom
      if (natom.gt.maxatom) then
        print *,"natom (",natom,") > maxatom (",maxatom,")"
        stop
      endif
      read(5,*)(iagroup(k),k=1,natom)
      read(5,*)(symb(k),k=1,natom)
      read(5,*)ipow,ipowt
      read(5,*)lreadbasis
      read(5,*)lreaddisc
      npairs=natom*(natom-1)/2            ! Number of interatomic pairs
      if (npairs.gt.maxpair) then
        print *,"npairs (",npairs,") > maxpair (",maxpair,")"
        stop
      endif

      write(6,*)
      write(6,*)"Atoms"
      write(6,'(1x,a9,100a5)')"Symbol",(symb(i),i=1,natom)
      write(6,'(a10,100i5)')"Index",(i,i=1,natom)
      write(6,'(a10,100i5)')"Group",(iagroup(i),i=1,natom)


ccc GENERATE BASIS ccc
      if (.not.lreadbasis) then
ccc GENERATE BASIS ccc


c     generate atom permutation lists
      do i=1,natom
        nngroup(i)=0
      enddo
      ngroup=1
      do i=1,natom
        if (iagroup(i).gt.ngroup) ngroup=iagroup(i)
        nngroup(iagroup(i))=nngroup(iagroup(i))+1
      enddo

      nn=0

      do i=1,ngroup
        n=0
        do k=1,natom
          if (iagroup(k).eq.i) then
            n=n+1
            idum(n)=k
          endif
        enddo
      
        npermute=0
        call heapp(idum,n,n,idum2,npermute)
        nperm0(i)=nn+1
        nperm1(i)=nn+npermute
        do k=1,npermute
          nn=nn+1
          m=0
          do j=1,natom
            idum3(nn,j)=0
            if (iagroup(j).eq.i) then
              m=m+1
              idum3(nn,j)=idum2(k,m)
            endif
          enddo
        enddo
      enddo ! i=1,ngroup

      ntmp=1
      do i=1,ngroup
        idum(i)=nperm0(i)
        print *,"Group ",i," has ",(nperm1(i)-nperm0(i)+1),
     & " permutations"
        ntmp=ntmp*(nperm1(i)-nperm0(i)+1)
      enddo
      print *,"For a total of ",ntmp," permutations"

      npermute=0
      do while (.true.)
        npermute=npermute+1
        if (npermute.gt.maxperm) then
          print *,"npermute (",npermute,") > maxperm (",maxperm,")"
          print *,"NOTE: maxperm needs to be at least npermute + 1"
          stop
        endif

        do i=1,natom
          iatom(npermute,i)=0
          do j=1,ngroup
            iatom(npermute,i)=iatom(npermute,i)+idum3(idum(j),i)
          enddo
        enddo

        idum(ngroup)=idum(ngroup)+1
 777    continue

        do i=1,ngroup
          if (idum(i).gt.nperm1(i)) then
            if (i.eq.1) go to 778
            idum(i)=nperm0(i)
            idum(i-1)=idum(i-1)+1
            go to 777
          endif
        enddo 

      enddo ! while (.true.)
 778  continue

      print *
      print *,'Atom permutations',npermute
      do i=1,min(npermute,100)
        print *,i,":",(iatom(i,j),j=1,natom)
      enddo

      ii=0
      do i=1,natom
        do j=i+1,natom
          ii=ii+1
          index(i,j)=ii
          index(j,i)=ii
        enddo
      enddo

      write(6,*)
      write(6,*)"Pair permutations"
      write(6,'(22x,100(a3,"- ",a3,4x))')
     &   ((symb(i),symb(j),j=1+i,natom),i=1,natom) 
      write(6,'(21x,100(i3," -",i3,4x))')((i,j,j=1+i,natom),
     &   i=1,natom) 
      do ii=1,npermute
        iix=0
        do i=1,natom
          do j=i+1,natom
            iix=iix+1
            ix(ii,iix)=index(iatom(ii,i),iatom(ii,j))
          enddo
        enddo
        if (ii.le.100) print *,ii,":",(ix(ii,iix),iix=1,npairs)
      enddo

c generate terms using individual power constraints
      ii=1
      do i=1,npairs
        ind(ii,i)=0
      enddo
      do while (.true.)
        ii=ii+1
        if (ii.gt.maxterm) then
          print *,"number of terms (",ii,") > maxterm (",maxterm,")"
          stop
        endif

        do i=1,npairs
          ind(ii,i)=ind(ii-1,i)
        enddo
        ind(ii,npairs)=ind(ii,npairs)+1
 300    continue

        indtot=0
        do i=1,npairs
          indtot=indtot+ind(ii,i)
          if (ind(ii,i).gt.ipow.or.indtot.gt.ipowt) then ! ipow(i) would allow atom-atom-type-dependent limits
            if (i.eq.1) go to 400
            ind(ii,i)=0
            ind(ii,i-1)=ind(ii,i-1)+1
            go to 300
          endif
        enddo
      enddo

 400  continue
      nterms=ii-1

c      print *
c      print *,"Basis # (Group):  Powers"

c symmetrize
      nbasis=0
      DO ii=1,nterms
        if (mod(ii,100).eq.0) print *,"symmetrizing ",ii," / ",nterms
        ifail=0
        do i=1,ii-1
          do j=1,npermute
            ifail=1
            do k=1,npairs
              if (ind(i,k).ne.ind(ii,ix(j,k))) ifail=0
            enddo
            if (ifail.eq.1) go to 1010
          enddo
        enddo
 1010 continue

        if (ifail.eq.0) then
          nbasis=nbasis+1
          ibasis(ii)=nbasis
        else
          ibasis(ii)=ibasis(i)
        endif
c      write(6,'(i5,"  (",i5,"):",100i8)')
c     &   ii,ibasis(ii),(ind(ii,j),j=1,npairs)
      ENDDO

      nncoef=nbasis

c remove disconnected terms
c      IF (.FALSE.) THEN   
      IF (lreaddisc) THEN
        ndisc=0
        do ii=1,nterms
          npx=0
          ipx=0
          do i=1,npairs
            if (ind(ii,i).ne.0) then
              npx=npx+1
              ipx=ipx+i
            endif
          enddo
c        !!! THIS TEST IS SPECIFIC TO OHHH !!!
          if (npx.eq.2.and.ipx.eq.7) then
c          print *,"unconnected"
c      write(6,'(i5,"  (",i5,"):",100i8)')
c     &   ii,ibasis(ii),(ind(ii,j),j=1,npairs)
            ndupe=0
            do j=1,ndisc
              if (idisc(j).eq.ibasis(ii)) ndupe=1
            enddo
            if (ndupe.eq.0) then
              ndisc=ndisc+1
              idisc(ndisc)=ibasis(ii)
            endif
          endif
        enddo
        print *
        print *,"disconncted groups: ",(idisc(j),j=1,ndisc)
      ENDIF ! TURN OFF
c end remove disconnected terms, can delete this whole section to turn this off

      print *
      print *,"Basis # (Group):  Powers"
 
      do ii=nterms,1,-1
        ibad=0
        do j=1,ndisc
          if (ibasis(ii).eq.idisc(j)) ibad=1
        enddo
        if (ibad.eq.1) then
          nterms=nterms-1
          do jj=ii,nterms
            ibasis(jj)=ibasis(jj+1)
            do k=1,npairs
              ind(jj,k)=ind(jj+1,k)
            enddo
          enddo
        endif
      enddo
      do ii=1,nterms
        nx=ibasis(ii)
        do i=1,ndisc
          if (nx.gt.idisc(i)) ibasis(ii)=ibasis(ii)-1
        enddo
      enddo
      do ii=1,nterms
        write(6,'(i5,"  (",i5,"):",100i8)')
     &   ii,ibasis(ii),(ind(ii,j),j=1,npairs)
      enddo
      nbasis=nbasis-ndisc
      nncoef=nbasis
      print *,'nncoef = ',nncoef
      open(55,file="basis.dat")
      write(55,*)natom,npairs,nncoef,nterms,
     & " ! atom pairs, coefficients, terms"
      write(55,*)" TERM GROUP :     EXPONENTS"
      do ii=1,nterms
        write(55,'(2i6," : ",1000i5)')
     &   ii,ibasis(ii),(ind(ii,j),j=1,npairs)
      enddo
      close(55)

ccc READ BASIS ccc
      else
      open(55,file="basis.dat")
      read(55,*)natom,npairs,nncoef,nterms
      read(55,*)
      do i=1,nterms
        read(55,*)k,ibasis(k),dum,(ind(k,j),j=1,npairs)
      enddo
      close(55)

      endif

      return

      entry funcs1(iii,basis,ncoef)

      do j=1,ncoef
        basis(j)=0.d0
      enddo

      do j=1,npairs
        r(j)=dexp(-rrr(iii,j))
      enddo

      do i=1,nterms
        arg=1.d0
        do j=1,npairs
          arg=arg*(r(j)**ind(i,j))
c          if (iii.eq.1) print *,"lll",i,j,arg
        enddo
        basis(ibasis(i))=basis(ibasis(i))+arg
      enddo

      return 
      end

***************************************************

      recursive subroutine heapp(ia,size,n,iia,ii)

      include 'param.inc'
      integer i,n,size,ii
      integer ia(maxatom)
      integer iia(maxperm,maxatom)
      integer iagroup(maxatom)

      if (size.eq.1) then
        ii=ii+1
        do i=1,n
          iia(ii,i)=ia(i)
        enddo
        return
      endif

      do i=1,size
        call heapp(ia,size-1,n,iia,ii)
        if (mod(size,2).eq.1) then
          tmp=ia(1)
          ia(1)=ia(size)
          ia(size)=tmp
        else
          tmp=ia(i)
          ia(i)=ia(size)
          ia(size)=tmp
        endif
      enddo

      end subroutine

***************************************************
