C  Compute C = s*A*B
C  where
C     C in n x m
C     A in n x k
C     B in k x m
C     s scalar
C
C  ***** Size Parameters *****
C  Matrix Dimensions
C     n  Rows of C and A
C     m  Columns of C and B
C     k  Columns of A, Rows of B
C  Memory Layout
C     sc Columnsize of C in memory
C     rc Row of C in memory
C     cc Column of C in memory
C     sa Columnsize of A in memory
C     ra Row of A in memory
C     ca Column of A in memory
C     sb Columnsize of B in memory
C     rb Row of B in memory
C     cb Column of B in memory

      subroutine musmm(n,m,k,sc,rc,cc,sa,ra,ca,sb,rb,cb,C,A,B,s)

      implicit none
      integer n,m,k,sc,rc,cc,sa,ra,ca,sb,rb,cb,i,j,l
      double precision C,A,B,s,temp
      dimension C(sc,*),A(sa,*),B(sb,*)

      do 1 i=1,n
         do 2 j=1,m
            temp=0.0D0
            do 3 l=1,k
               temp=temp+A(i+ra,l+ca)*B(l+rb,j+cb)
3           continue
            C(i+rc,j+cc)=s*temp
2        continue
1     continue

      end


C  Compute C += s*A*B
C  where
C     C in n x m
C     A in n x k
C     B in k x m
C     s scalar
C
C  ***** Size Parameters *****
C  Matrix Dimensions
C     n  Rows of C and A
C     m  Columns of C and B
C     k  Columns of A, Rows of B
C  Memory Layout
C     sc Columnsize of C in memory
C     rc Row of C in memory
C     cc Column of C in memory
C     sa Columnsize of A in memory
C     ra Row of A in memory
C     ca Column of A in memory
C     sb Columnsize of B in memory
C     rb Row of B in memory
C     cb Column of B in memory

      subroutine masmm(n,m,k,sc,rc,cc,sa,ra,ca,sb,rb,cb,C,A,B,s)

      implicit none
      integer n,m,k,sc,rc,cc,sa,ra,ca,sb,rb,cb,i,j,l
      double precision C,A,B,s,temp
      dimension C(sc,*),A(sa,*),B(sb,*)

      do 1 i=1,n
         do 2 j=1,m
            temp=0.0D0
            do 3 l=1,k
               temp=temp+A(i+ra,l+ca)*B(l+rb,j+cb)
3           continue
            C(i+rc,j+cc)=C(i+rc,j+cc)+s*temp
2        continue
1     continue

      end


C  Compute C = s*A*B'
C  where
C     C in n x m
C     A in n x k
C     B in m x k
C     s scalar
C
C  ***** Size Parameters *****
C  Matrix Dimensions
C     n  Rows of C and A
C     m  Columns of C and Rows of B
C     k  Columns of A and B
C  Memory Layout
C     sc Columnsize of C in memory
C     rc Row of C in memory
C     cc Column of C in memory
C     sa Columnsize of A in memory
C     ra Row of A in memory
C     ca Column of A in memory
C     sb Columnsize of B in memory
C     rb Row of B in memory
C     cb Column of B in memory

      subroutine musmmt(n,m,k,sc,rc,cc,sa,ra,ca,sb,rb,cb,C,A,B,s)

      implicit none
      integer n,m,k,sc,rc,cc,sa,ra,ca,sb,rb,cb,i,j,l
      double precision C,A,B,s,temp
      dimension C(sc,*),A(sa,*),B(sb,*)

      do 1 i=1,n
         do 2 j=1,m
            temp=0.0D0
            do 3 l=1,k
               temp=temp+A(i+ra,l+ca)*B(j+rb,l+cb)
3           continue
            C(i+rc,j+cc)=s*temp
2        continue
1     continue

      end



C  Compute C += s*A*B'
C  where
C     C in n x m
C     A in n x k
C     B in m x k
C     s scalar
C
C  ***** Size Parameters *****
C  Matrix Dimensions
C     n  Rows of C and A
C     m  Columns of C and Rows of B
C     k  Columns of A and B
C  Memory Layout
C     sc Columnsize of C in memory
C     rc Row of C in memory
C     cc Column of C in memory
C     sa Columnsize of A in memory
C     ra Row of A in memory
C     ca Column of A in memory
C     sb Columnsize of B in memory
C     rb Row of B in memory
C     cb Column of B in memory

      subroutine masmmt(n,m,k,sc,rc,cc,sa,ra,ca,sb,rb,cb,C,A,B,s)

      implicit none
      integer n,m,k,sc,rc,cc,sa,ra,ca,sb,rb,cb,i,j,l
      double precision C,A,B,s,temp
      dimension C(sc,*),A(sa,*),B(sb,*)

      do 1 i=1,n
         do 2 j=1,m
            temp=0.0D0
            do 3 l=1,k
               temp=temp+A(i+ra,l+ca)*B(j+rb,l+cb)
3           continue
            C(i+rc,j+cc)=C(i+rc,j+cc)+s*temp
2        continue
1     continue

      end




C  Compute C = s*A*A'
C  where
C     C in n x n
C     A in n x m
C     s scalar
C
C  ***** Size Parameters *****
C  Matrix Dimensions
C     n  Rows and Columns of C, Rows of A
C     m  Columns of A
C  Memory Layout
C     sc Columnsize of C in memory
C     rc Row of C in memory
C     cc Column of C in memory
C     sa Columnsize of A in memory
C     ra Row of A in memory
C     ca Column of A in memory

      subroutine musmt(n,m,sc,rc,cc,sa,ra,ca,C,A,s)

      implicit none
      integer n,m,k,sc,rc,cc,sa,ra,ca,i,j,l
      double precision C,A,s,temp
      dimension C(sc,*),A(sa,*)

      do 1 i=1,n
         temp=0.0D0
         do 3 l=1,m
            temp=temp+A(i+ra,l+ca)**2
3        continue
         C(i+rc,i+cc)=s*temp
         do 2 j=i+1,n
            temp=0.0D0
            do 4 l=1,m
               temp=temp+A(i+ra,l+ca)*A(j+ra,l+ca)
4           continue
            temp = s*temp;
            C(i+rc,j+cc)=temp
            C(j+rc,i+cc)=temp
2        continue
1     continue

      end


C  Compute C += s*A*A'
C  where
C     C in n x n and originally symmetric
C     A in n x m
C     s scalar
C
C  ***** Size Parameters *****
C  Matrix Dimensions
C     n  Rows and Columns of C, Rows of A
C     m  Columns of A
C  Memory Layout
C     sc Columnsize of C in memory
C     rc Row of C in memory
C     cc Column of C in memory
C     sa Columnsize of A in memory
C     ra Row of A in memory
C     ca Column of A in memory

      subroutine masmt(n,m,sc,rc,cc,sa,ra,ca,C,A,s)

      implicit none
      integer n,m,k,sc,rc,cc,sa,ra,ca,i,j,l
      double precision C,A,s,temp
      dimension C(sc,*),A(sa,*)

      do 1 i=1,n
         temp=0.0D0
         do 3 l=1,m
            temp=temp+A(i+ra,l+ca)**2
3        continue
         C(i+rc,i+cc)=C(i+rc,i+cc)+s*temp
         do 2 j=i+1,n
            temp=0.0D0
            do 4 l=1,m
               temp=temp+A(i+ra,l+ca)*A(j+ra,l+ca)
4           continue
            temp = C(i+rc,j+cc) + s*temp;
            C(i+rc,j+cc)=temp
            C(j+rc,i+cc)=temp
2        continue
1     continue

      end



C  Compute C = s*A*D*A'
C  where
C     C in n x n
C     A in n x m
C     D in m x m diagonal
C     s scalar
C
C  ***** Size Parameters *****
C  Matrix Dimensions
C     n  Rows and Columns of C, Rows of A
C     m  Rows and Columns of D, Columns of A
C  Memory Layout
C     sc Columnsize of C in memory
C     rc Row of C in memory
C     cc Column of C in memory
C     sa Columnsize of A in memory
C     ra Row of A in memory
C     ca Column of A in memory
C     rd Index of D in memory (D is vector with diagonal elements)

      subroutine musmtd(n,m,sc,rc,cc,sa,ra,ca,rd,C,A,D,s)

      implicit none
      integer n,m,k,sc,rc,cc,sa,ra,ca,rd,i,j,l
      double precision C,A,D,s,temp
      dimension C(sc,*),A(sa,*),D(*)

      do 1 i=1,n
         temp=0.0D0
         do 3 l=1,m
            temp=temp+D(rd+l)*A(i+ra,l+ca)**2
3        continue
         C(i+rc,i+cc)=s*temp
         do 2 j=i+1,n
            temp=0.0D0
            do 4 l=1,m
               temp=temp+A(i+ra,l+ca)*D(rd+l)*A(j+ra,l+ca)
4           continue
            temp = s*temp;
            C(i+rc,j+cc)=temp
            C(j+rc,i+cc)=temp
2        continue
1     continue

      end



C  Compute C += s*A*D*A'
C  where
C     C in n x n and originally symmetric
C     A in n x m
C     D in m x m diagonal
C     s scalar
C
C  ***** Size Parameters *****
C  Matrix Dimensions
C     n  Rows and Columns of C, Rows of A
C     m  Rows and Columns of D, Columns of A
C  Memory Layout
C     sc Columnsize of C in memory
C     rc Row of C in memory
C     cc Column of C in memory
C     sa Columnsize of A in memory
C     ra Row of A in memory
C     ca Column of A in memory
C     rd Index of D in memory (D is vector with diagonal elements)

      subroutine masmtd(n,m,sc,rc,cc,sa,ra,ca,rd,C,A,D,s)

      implicit none
      integer n,m,k,sc,rc,cc,sa,ra,ca,rd,i,j,l
      double precision C,A,D,s,temp
      dimension C(sc,*),A(sa,*),D(*)

      do 1 i=1,n
         temp=0.0D0
         do 3 l=1,m
            temp=temp+D(rd+l)*A(i+ra,l+ca)**2
3        continue
         C(i+rc,i+cc)=C(i+rc,i+cc)+s*temp
         do 2 j=i+1,n
            temp=0.0D0
            do 4 l=1,m
               temp=temp+A(i+ra,l+ca)*D(rd+l)*A(j+ra,l+ca)
4           continue
            temp = C(i+rc,j+cc) + s*temp;
            C(i+rc,j+cc)=temp
            C(j+rc,i+cc)=temp
2        continue
1     continue

      end



C  Compute C = s*A*D
C  where
C     C in n x m
C     A in n x m
C     D in m x m diagonal
C     s scalar
C
C  ***** Size Parameters *****
C  Matrix Dimensions
C     n  Rows of C, Rows of A
C     m  Rows and Columns of D, Columns of A and C
C  Memory Layout
C     sc Columnsize of C in memory
C     rc Row of C in memory
C     cc Column of C in memory
C     sa Columnsize of A in memory
C     ra Row of A in memory
C     ca Column of A in memory
C     rd Index of D in memory (D is vector with diagonal elements)

      subroutine musmd(n,m,sc,rc,cc,sa,ra,ca,rd,C,A,D,s)

      implicit none
      integer n,m,sc,rc,cc,sa,ra,ca,rd,i,j,l
      double precision C,A,D,s,temp
      dimension C(sc,*),A(sa,*),D(*)

      do 1 i=1,n
         do 2 j=1,m
            C(i+rc,j+cc)=s*D(rd+j)*A(i+ra,j+ca);
2        continue
1     continue

      end



C  Compute C += s*A*D
C  where
C     C in n x m
C     A in n x m
C     D in m x m diagonal
C     s scalar
C
C  ***** Size Parameters *****
C  Matrix Dimensions
C     n  Rows of C, Rows of A
C     m  Rows and Columns of D, Columns of A and C
C  Memory Layout
C     sc Columnsize of C in memory
C     rc Row of C in memory
C     cc Column of C in memory
C     sa Columnsize of A in memory
C     ra Row of A in memory
C     ca Column of A in memory
C     rd Index of D in memory (D is vector with diagonal elements)

      subroutine masmd(n,m,sc,rc,cc,sa,ra,ca,rd,C,A,D,s)

      implicit none
      integer n,m,sc,rc,cc,sa,ra,ca,rd,i,j
      double precision C,A,D,s,temp
      dimension C(sc,*),A(sa,*),D(*)

      do 1 i=1,n
         do 2 j=1,m
            C(i+rc,j+cc)=C(i+rc,j+cc) + s*D(rd+j)*A(i+ra,j+ca);
2        continue
1     continue

      end



C  Compute C = s*A*D*B
C  where
C     C in n x m
C     A in n x k
C     B in k x m
C     D in k x k diagonal
C     s scalar
C
C  ***** Size Parameters *****
C  Memory Layout
C     sc Columnsize of C in memory
C     rc Row of C in memory
C     cc Column of C in memory
C     sa Columnsize of A in memory
C     ra Row of A in memory
C     ca Column of A in memory
C     sb Columnsize of B in memory
C     rb Row of B in memory
C     cb Column of B in memory
C     rd Index of D in memory (D is vector with diagonal elements)

      subroutine musmdm(n,m,k,sc,rc,cc,sa,ra,ca,sb,rb,cb,rd,C,A,B,D,s)

      implicit none
      integer n,m,k,sc,rc,cc,sa,ra,ca,sb,rb,cb,rd,i,j,l
      double precision C,A,B,D,s,temp
      dimension C(sc,*),A(sa,*),B(sb,*),D(*)

      do 1 i=1,n
         do 2 j=1,m
            temp=0.0D0
            do 4 l=1,k
               temp=temp+A(i+ra,l+ca)*D(rd+l)*B(l+rb,j+cb)
4           continue
            C(i+rc,j+cc)=s*temp
2        continue
1     continue

      end



C  Compute C += s*A*D*B
C  where
C     C in n x m
C     A in n x k
C     B in k x m
C     D in k x k diagonal
C     s scalar
C
C  ***** Size Parameters *****
C  Memory Layout
C     sc Columnsize of C in memory
C     rc Row of C in memory
C     cc Column of C in memory
C     sa Columnsize of A in memory
C     ra Row of A in memory
C     ca Column of A in memory
C     sb Columnsize of B in memory
C     rb Row of B in memory
C     cb Column of B in memory
C     rd Index of D in memory (D is vector with diagonal elements)

      subroutine masmdm(n,m,k,sc,rc,cc,sa,ra,ca,sb,rb,cb,rd,C,A,B,D,s)

      implicit none
      integer n,m,k,sc,rc,cc,sa,ra,ca,sb,rb,cb,rd,i,j,l
      double precision C,A,B,D,s,temp
      dimension C(sc,*),A(sa,*),B(sb,*),D(*)

      do 1 i=1,n
         do 2 j=1,m
            temp=0.0D0
            do 4 l=1,k
               temp=temp+A(i+ra,l+ca)*D(rd+l)*B(l+rb,j+cb)
4           continue
            C(i+rc,j+cc)=C(i+rc,j+cc)+s*temp
2        continue
1     continue

      end



C  Compute C = s*A
C  where
C     C in n x m
C     A in n x m
C     s scalar
C
C  ***** Size Parameters *****
C  Matrix Dimensions
C     n  Rows of C and A
C     m  Columns of A and C
C  Memory Layout
C     sc Columnsize of C in memory
C     rc Row of C in memory
C     cc Column of C in memory
C     sa Columnsize of A in memory
C     ra Row of A in memory
C     ca Column of A in memory

      subroutine eqsm(n,m,sc,rc,cc,sa,ra,ca,C,A,s)

      implicit none
      integer n,m,k,sc,rc,cc,sa,ra,ca,rd,i,j
      double precision C,A,s
      dimension C(sc,*),A(sa,*)

      do 1 i=1,n
         do 2 j=1,m
            C(i+rc,j+cc)=s*A(i+ra,j+ca);
2        continue
1     continue

      end




C  Compute C = s*A
C  where
C     C in n x m
C     A in n x m
C     s scalar
C
C  ***** Size Parameters *****
C  Matrix Dimensions
C     n  Rows of C and A
C     m  Columns of A and C
C  Memory Layout
C     sc Columnsize of C in memory
C     rc Row of C in memory
C     cc Column of C in memory
C     sa Columnsize of A in memory
C     ra Row of A in memory
C     ca Column of A in memory

      subroutine adsm(n,m,sc,rc,cc,sa,ra,ca,C,A,s)

      implicit none
      integer n,m,k,sc,rc,cc,sa,ra,ca,rd,i,j
      double precision C,A,s,temp
      dimension C(sc,*),A(sa,*)

      do 1 i=1,n
         do 2 j=1,m
            C(i+rc,j+cc)=C(i+rc,j+cc) + s*A(i+ra,j+ca);
2        continue
1     continue

      end



C  Solve A*X=B
C  where
C     A in n x n, symmetrical, positiv definite
C     X,B in n x m
C  The solution X will be stored in B
C
C  ***** Size Parameters *****
C  Matrix Dimensions
C     n  Rows and Columns of A, Rows of X and B
C     m Columns of X and B
C  Memory Layout
C     sb Column size of X and B
C     rb Row of X and B in memory
C     cb Column of X and B in memory

      subroutine choslv(n,m,sb,rb,A,B,err)

      implicit none
      integer n,m,sb,rb,i,j,k,err
      double precision A,B,diag
      dimension A(n,n),B(sb,m)

      err=3
      do 1 i=1,n-1
         diag = A(i,i);
         if (diag.le.0.0D0) return
         diag = dsqrt(diag)
         A(i,i) = diag
         do 2 k=1,m
            B(rb+i,k) = B(rb+i,k)/diag;
2        continue

         do 3 j=i+1,n
            A(j,i) = A(j,i)/diag
            do 4 k=1,m
               B(rb+j,k) = B(rb+j,k) - B(rb+i,k)*A(j,i)
4           continue
3        continue

         do 5 j=i+1,n
            do 6 k=i+1,j
               A(j,k) = A(j,k) - A(j,i)*A(k,i)
6           continue
5        continue
1     continue

      diag = A(n,n)
      do 7 k=1,m
         B(rb+n,k) = B(rb+n,k)/diag;
7     continue
      A(n,n)=sqrt(diag)

      do 8 i=n,2,-1
         do 9 j=1,i-1
            do 10 k=1,m
               B(rb+j,k) = B(rb+j,k) - B(rb+i,k)*A(i,j);
10          continue
9        continue
         do 11 k=1,m
               B(rb+i-1,k) = B(rb+i-1,k)/A(i-1,i-1);
11       continue
8     continue

      err = 0

      end


C  Solve A*X=B
C  where
C     A in n x n, result of RR decomposition
C     X,B in n x m
C  The solution X will be stored in B
C
C  ***** Size Parameters *****
C  Matrix Dimensions
C     n  Rows and Columns of A, Rows of X and B
C     m Columns of X and B
C  Memory Layout
C     sb Column size of X and B
C     rb Row of X and B in memory
C     cb Column of X and B in memory

      subroutine rrslv(n,m,sb,rb,A,B)

      implicit none
      integer n,m,sb,rb,i,j,k
      double precision A,B
      dimension A(n,n),B(sb,m)

      do i=1,n
         do j=1,i-1
            do k=1,m
               B(rb+i,k) = B(rb+i,k) - B(rb+j,k)*A(i,j);
            end do
         end do
         do k=1,m
               B(rb+i,k) = B(rb+i,k)/A(i,i);
         end do
      end do

      do i=n,1,-1
         do j=i+1,n
            do k=1,m
               B(rb+i,k) = B(rb+i,k) - B(rb+j,k)*A(j,i);
            end do
         end do
         do k=1,m
               B(rb+i,k) = B(rb+i,k)/A(i,i);
         end do
      end do

      end


      subroutine chlsch(n,subn,m,sb,rb,A,B,err,
     /            ns3,ns6,nbs,sdiag,ss3,ss6,sbs)

      implicit none
      integer n,m,sb,rb,i,j,k,err,totaln,subn
      integer ns3,ns6,nbs,sdiag,ss3,ss6,sbs
      integer cns3,cns6,cbs,i1,j1,i2,j2,i3,j3
      double precision A,B,diag
      dimension A(subn,(2*n-1)*subn),B(sb,m)
      dimension ns3(n*subn-1),ns6(n*subn-1),nbs(n*subn-1)
      dimension sdiag(2,*),ss3(3,*),ss6(6,*),sbs(3,*)

      totaln=n*subn
      err=3
      cns3=1
      cns6=1
      cbs=1

      do 1 i=1,totaln-1
         i1 = sdiag(1,i)
         j1 = sdiag(2,i);
         diag = A(i1,j1);
         if (diag.le.0.0D0) return
         diag = dsqrt(diag)
         A(i1,j1) = diag
         do 2 k=1,m
            B(rb+i,k) = B(rb+i,k)/diag;
2        continue

         do 3 j=1,ns3(i)
            i3 = ss3(1,cns3)
            i2 = ss3(2,cns3)
            j2 = ss3(3,cns3)
            cns3=cns3+1
            A(i2,j2) = A(i2,j2)/diag
            do 3 k=1,m
               B(rb+i3,k)=B(rb+i3,k)-B(rb+i,k)*A(i2,j2)
3        continue

         do 1 j=1,ns6(i)
            i1=ss6(1,cns6)
            j1=ss6(2,cns6)
            i2=ss6(3,cns6)
            j2=ss6(4,cns6)
            i3=ss6(5,cns6)
            j3=ss6(6,cns6)
            cns6=cns6+1
            A(i1,j1)=A(i1,j1)-A(i2,j2)*A(i3,j3)
1     continue

      diag = A(sdiag(1,totaln),sdiag(2,totaln))
*      if (diag.le.0.0D0) return
      do 7 k=1,m
         B(rb+totaln,k) = B(rb+totaln,k)/diag;
7     continue


      j3=1
      do 8 i=totaln,2,-1
         do 9 j=1,nbs(j3)
            i1=sbs(1,cbs)
            i2=sbs(2,cbs)
            j2=sbs(3,cbs)
            cbs=cbs+1
            do 9 k=1,m
               B(rb+i1,k)=B(rb+i1,k)-B(rb+i,k)*A(i2,j2);
9        continue
         j3=j3+1
         diag=A(sdiag(1,i-1),sdiag(2,i-1))
*         if (diag.le.0.0D0) return
         do 8 k=1,m
            B(rb+i-1,k)=B(rb+i-1,k)/diag;
8     continue

      err = 0
      end



      subroutine gauss(n,m,A,B,P)

      implicit none
      integer n,m,P,i,j,k,pivot,r1,r2
      double precision A,B,alpha
      dimension A(n,n),B(n,m)
      dimension P(n)

      do 1 i=1,n
         P(i)=i
1     continue

      do 2 i=1,n-1
         pivot = i
         do 3 j=i+1,n
            if (dabs(A(P(j),i)) .gt. dabs(A(P(pivot),i))) then
               pivot=j
            end if
3        continue

         if (pivot .ne. i) then
            j=P(i)
            P(i)=P(pivot)
            P(pivot)=j
         end if

         do 4 j=i+1,n
            r1=P(i)
            r2=P(j)
            alpha=A(r2,i)/A(r1,i)
            do 5 k=i+1,n
               A(r2,k)=A(r2,k)-alpha*A(r1,k)
5           continue
            do 6 k=1,m
               B(r2,k)=B(r2,k)-alpha*B(r1,k)
6           continue
4        continue
2     continue

      do 7 i=n,1,-1
         r1=P(i)
         do 8 j=i+1,n
            do 9 k=1,m
               B(r1,k)=B(r1,k)-B(P(j),k)*A(r1,j)
9           continue
8        continue
         do 10 k=1,m
            B(r1,k)=B(r1,k)/A(r1,i)
10        continue
7     continue

      end
