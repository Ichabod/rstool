      subroutine chlsco(n,m,sb,rb,A,B,err,pinfo)

      implicit none
      integer n,m,sb,rb,err,pinfo
      double precision A,B
      dimension A(*),B(*)
      dimension pinfo(*)

      call chloco(n,m,sb,rb,A,B,err,
     /  pinfo(pinfo(2)),pinfo(pinfo(3)),pinfo(pinfo(4)),
     /  pinfo(pinfo(1)),pinfo(pinfo(5)),pinfo(pinfo(6)),
     /  pinfo(pinfo(7)))

      end


      subroutine chloco(n,m,sb,rb,A,B,err,
     /            ns2,ns3,nbs,sdiag,ss2,ss3,sbs)

      implicit none
      integer n,m,sb,rb,i,j,k,err
      integer ns2,ns3,nbs,sdiag,ss2,ss3,sbs
      integer cns2,cns3,cbs,i1,j1,i2,j2
      double precision A,B,diag
      dimension A(n,n),B(sb,m)
      dimension ns2(n-1),ns3(n-1),nbs(n-1)
      dimension sdiag(*),ss2(*),ss3(2,*),sbs(*)

      err=3
      cns2=1
      cns3=1
      cbs=1

      do 1 i=1,n-1
         i1 = sdiag(i)
         diag = A(i1,i1);
         if (diag.le.0.0D0) return
         diag = dsqrt(diag)
         A(i1,i1) = diag
         do 2 k=1,m
            B(rb+i1,k) = B(rb+i1,k)/diag;
2        continue

         do 3 j=1,ns2(i)
            i2 = ss2(cns2)
            cns2=cns2+1
            A(i2,i1) = A(i2,i1)/diag
            do 3 k=1,m
               B(rb+i2,k)=B(rb+i2,k)-B(rb+i1,k)*A(i2,i1)
3        continue

         do 1 j=1,ns3(i)
            i2=ss3(1,cns3)
            j2=ss3(2,cns3)
            cns3=cns3+1
            A(i2,j2)=A(i2,j2)-A(i2,i1)*A(j2,i1)
1     continue

      i1 = sdiag(n)
      diag = A(i1,i1)
      do 7 k=1,m
         B(rb+i1,k) = B(rb+i1,k)/diag;
7     continue
      A(i1,i1) = sqrt(diag);

      i1=sdiag(n)
      do 8 i=n,2,-1
         do 9 j=1,nbs(n-i+1)
            i2=sbs(cbs)
            cbs=cbs+1
            do 9 k=1,m
               B(rb+i2,k)=B(rb+i2,k)-B(rb+i1,k)*A(i1,i2);
9        continue
         i1=sdiag(i-1)
         diag=A(i1,i1)
         do 8 k=1,m
            B(rb+i1,k)=B(rb+i1,k)/diag;
8     continue

      err = 0
      end


      subroutine rrslvo(n,m,sb,rb,A,B,pinfo)

      implicit none
      integer n,m,sb,rb,pinfo
      double precision A,B
      dimension A(*),B(*)
      dimension pinfo(*)

      call rroco(n,m,sb,rb,A,B,
     /  pinfo(pinfo(2)),pinfo(pinfo(4)),
     /  pinfo(pinfo(1)),pinfo(pinfo(5)),
     /  pinfo(pinfo(7)))

      end

      subroutine rroco(n,m,sb,rb,A,B,
     /            ns2,nbs,sdiag,ss2,sbs)

      implicit none
      integer n,m,sb,rb,i,j,k
      integer ns2,nbs,sdiag,ss2,sbs
      integer cns2,cbs,i1,j1,i2
      double precision A,B,diag
      dimension A(n,n),B(sb,m)
      dimension ns2(n-1),nbs(n-1)
      dimension sdiag(*),ss2(*),sbs(*)

      cns2=1
      cbs=1

      do 1 i=1,n-1
         i1 = sdiag(i)
         diag = A(i1,i1);
         do 2 k=1,m
            B(rb+i1,k) = B(rb+i1,k)/diag;
2        continue

         do 3 j=1,ns2(i)
            i2 = ss2(cns2)
            cns2=cns2+1
            do 3 k=1,m
               B(rb+i2,k)=B(rb+i2,k)-B(rb+i1,k)*A(i2,i1)
3        continue
1     continue

      i1 = sdiag(n)
      diag = A(i1,i1)
      do 7 k=1,m
         B(rb+i1,k) = B(rb+i1,k)/(diag*diag);
7     continue

      i1=sdiag(n)
      do 8 i=n,2,-1
         do 9 j=1,nbs(n-i+1)
            i2=sbs(cbs)
            cbs=cbs+1
            do 9 k=1,m
               B(rb+i2,k)=B(rb+i2,k)-B(rb+i1,k)*A(i1,i2);
9        continue
         i1=sdiag(i-1)
         diag=A(i1,i1)
         do 8 k=1,m
            B(rb+i1,k)=B(rb+i1,k)/diag;
8     continue

      end
