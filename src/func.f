********************************************

      subroutine prepot

      implicit double precision(a-h,p-z)
      include 'param.inc'
      integer c,dtot,dpcheck,discind,pairpair,npow
      integer dgroup,dg1,dg2,ndg,dgtot
!      integer OMP_GET_MAX_THREADS
!      integer OMP_GET_NUM_THREADS
!      integer OMP_GET_THREAD_NUM
      dimension iagroup(maxatom),
     &  ind(maxterm,maxpair),
     &  iatom(maxperm,maxatom),
     &  idum(maxatom),nngroup(maxatom),idisc(maxterm),
     &  idum2(maxperm,maxatom),
     &  idum3(maxperm,maxatom),nperm0(maxperm),nperm1(maxperm),
     &  basis(maxterm),ibasis(maxterm),r(maxpair),
     &  dbasisdr(maxterm,maxpair),
     &  rrr(maxdata,maxpair),index(maxatom,maxatom),ix(maxperm,maxpair),
     &  dpcheck(maxpair),discind(maxatom,maxatom),
     &  pairpair(maxpair,maxpair),npow(maxpair),dgroup(maxatom),
     &  ndg(maxatom)
      character*2 symb(maxatom),dum!,dsymb1,dsymb2
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
c DISC
      if (lreaddisc) then
        read(5,*)(dgroup(k),k=1,natom) ! disconnected groups of products  (e.g. 1 1 2 2) 
      endif

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
      IF (.not.lreadbasis) THEN
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
          enddo ! j=1,natom
        enddo ! k=1,npermute
      enddo ! i=1,ngroup

      ntmp=1
      do i=1,ngroup
        idum(i)=nperm0(i)
      print *,"Group ",i," has ",(nperm1(i)-nperm0(i)+1)," permutations"
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

      enddo
 778  continue

      print *
      print *,'Atom permutations',npermute
      do i=1,min(npermute,100) ! prints only first 100 permutations
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

c symmetrize
      nbasis=0
c!$OMP PARALLEL
c      write(6,"(2(a,i3))") " OpenMP: N_threads = ", 
c     &   OMP_GET_NUM_THREADS()," thread = ", OMP_GET_THREAD_NUM()
c!$OMP END PARALLEL
c!$OMP PARALLEL DO
      DO ii=1,nterms
        if (mod(ii,100).eq.0) print *,"symmetrizing ",ii," / ",nterms
c!$OMP THREADPRIVATE(ifail)
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
      ENDDO
c!$OMP END PARALLEL DO

      nncoef=nbasis

c DISC remove disconnected terms
c 1: Identify atom pairs that satisfy possible requirements
c 2: record atom indices of pairs
c 3: find other pair(s) (P2) sharing 0 or 1 (not both) atom types as pair P1
c 4: if neither atom of pair P2 has the same index as those in pair P1,
c     classify as disconnected term if only includes powers of these pairs

      if (lreaddisc) then
        dgtot=0         ! number of DGs
        do i=1,natom
          ndg(i)=0      ! number of atoms in each DG
        enddo
        do i=1,npairs      
          dpcheck(i)=0   ! becomes 1 if pair matches DG conditions
          discind(i,1)=0
          discind(i,2)=0
        enddo
        ndp=0        ! number of disconnected pairs

c       steps 1 and 2
        do i=1,natom
          ndg(dgroup(i))=ndg(dgroup(i))+1
          if (dgroup(i).gt.dgtot) dgtot=dgroup(i)
        enddo
        print *,(ndg(i),i=1,dgtot)

        DO dg1=1,dgtot-1
        if (ndg(dg1).le.1) then
          print *,"No disconnected groups possible"
        else
        DO dg2=dg1+1,dgtot
          print *,"Checking connections between DGs ",dg1," and ",dg2
          if (ndg(dg2).le.1) then
            print *,"No disconnected groups possible"
          else
            do i=1,natom
              do j=i+1,natom
                ag=0
                dg=0
                if (iagroup(i).eq.iagroup(j)) ag=1
                if (dgroup(i).eq.dgroup(j)) dg=1
                if (dg.eq.1) then ! if dg = 1, all good
                  ndp=ndp+1
                  dpcheck(index(i,j))=1
                  discind(ndp,1)=i
                  discind(ndp,2)=j
                  do k=1,natom ! inner loop over atoms, find matching AGs
                    if (iagroup(k).eq.iagroup(i)) then
                      do l=k+1,natom
                        if (iagroup(l).eq.iagroup(j)) then
                          if (dpcheck(index(k,l)).eq.0) then
                            ndp=ndp+1
                            dpcheck(index(k,l))=1
                            discind(ndp,1)=k
                            discind(ndp,2)=l
                          endif
                        endif
                      enddo
                    elseif (iagroup(k).eq.iagroup(j)) then
                      do l=k+1,natom
                        if (iagroup(l).eq.iagroup(i)) then
                          if (dpcheck(index(k,l)).eq.0) then
                            ndp=ndp+1
                            dpcheck(index(k,l))=1
                            discind(ndp,1)=k
                            discind(ndp,2)=l
                          endif
                        endif
                      enddo
                    endif
                  enddo
                elseif (ag.eq.1) then ! dg = 0, ag = 1, also all good
                  ndp=ndp+1
                  dpcheck(index(i,j))=1
                  discind(ndp,1)=i
                  discind(ndp,2)=j
                endif
              enddo
            enddo
          endif
        ENDDO
        endif
        ENDDO

        print *,"Molecular Groups"
        write(6,'(15x,100(2x,a2,"- ",a2))')
     &   ((symb(i),symb(j),j=1+i,natom),i=1,natom) 
        write(6,'(15x,100i8)')
     &   (dpcheck(i),i=1,npairs)
      ! should give dpcheck(1,2,3)=1 (O-H bond pairs)
c
c       step 3 - create pair-pair matrix -- PARALLELIZE
        do i=1,npairs
          npow(i)=0
          do j=1,npairs
            pairpair(i,j)=0
          enddo
        enddo

        do c=1,ndp ! disconnected pairs
          i=discind(c,1) ! pair 1, atom 1
          j=discind(c,2) ! pair 1, atom 2
          do k=i+1,natom ! pair 2, atom 1
            do l=k+1,natom ! pair 2, atom 2
              if (j.ne.k.and.j.ne.l) then
                pairpair(index(i,j),index(k,l))=1
                pairpair(index(k,l),index(i,j))=1
              endif
            enddo
          enddo
        enddo

        ndisc=0
        print *,"disconnected terms:"
        do ii=1,nterms ! basis terms
          npx=0 ! number of pairs
          ipx=0 ! total power of pairs
          do i=1,npairs ! permutation pairs
            if (ind(ii,i).ne.0) then ! if there are powers of ii
              npx=npx+1
              ipx=ipx+i
              npow(npx)=i
            endif
          enddo

c DISC  new general test
          if (npx.eq.2.and.pairpair(npow(1),npow(2)).eq.1) then ! pair is a monitored pair
            write(6,'(i5,"  (",i5,"):",100i8)')
     &      ii,ibasis(ii),(ind(ii,j),j=1,npairs)
!            write(6,'(i5,"  (",i5,"):",100i8)')
!     &      ii,ibasis(ii),(ind(ii,j),j=1,npairs)

            ndupe=0
            do j=1,ndisc
              if (idisc(j).eq.ibasis(ii)) ndupe=1
            enddo
            if (ndupe.eq.0) then
              ndisc=ndisc+1
              idisc(ndisc)=ibasis(ii)
            endif
          endif
        enddo ! ii=1,nterms
 
        print *
        print *,"Disconnected groups: ",(idisc(j),j=1,ndisc)
c remove disconnected terms
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
      endif ! lreaddisc
c end remove disconnected terms

      print *
      print *,"Basis # (Group):  Powers"

      do ii=1,nterms
        write(6,'(i5,"  (",i5,"):",100i8)')
     &   ii,ibasis(ii),(ind(ii,j),j=1,npairs)
      enddo

      nbasis=nbasis-ndisc
      nncoef=nbasis
      print *,'nncoef = ',nncoef

      open(55,file="basis.dat")
      write(55,*)natom,npairs,nncoef,nterms,
     & " ! atoms, atom pairs, coefficients, terms"
      write(55,*)" TERM GROUP :     EXPONENTS"
      do ii=1,nterms
      write(55,'(2i6," : ",1000i5)')
     &   ii,ibasis(ii),(ind(ii,j),j=1,npairs)
      enddo
      close(55)

ccc READ BASIS ccc
      ELSE
        open(55,file="basis.dat")
        read(55,*)natom,npairs,nncoef,nterms
        read(55,*)
        do i=1,nterms
          read(55,*)k,ibasis(k),dum,(ind(k,j),j=1,npairs)
        enddo
        close(55)
      ENDIF

      return

      entry funcs1(iii,basis,ncoef)
!      entry funcs1(iii,basis,ncoef,dbasisdr)

      !write(6,*)iii,basis,ncoef

      do j=1,ncoef
        basis(j)=0.d0
!        do j=1,npairs
!          dbasisdr(i,j)=0.d0
!        enddo
      enddo

      do j=1,npairs
        r(j)=dexp(-rrr(iii,j))
!        r(j)=dexp(-rrr(iii,j)*autoang)
      enddo

      do i=1,nterms
        arg=1.d0
        do j=1,npairs
          arg=arg*(r(j)**ind(i,j))
!          if (ind(i,j).ne.0) arg=arg*(r(j)**ind(i,j))
        enddo
        basis(ibasis(i))=basis(ibasis(i))+arg
!        do j=1,npairs
!          dbasisdr(ibasis(i),j)=dbasisdr(ibasis(i),j)   ! dV/dy * dy/dr
!     &                         -arg*dble(ind(i,j))*autoang
!        enddo
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
