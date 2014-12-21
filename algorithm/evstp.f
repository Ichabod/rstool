      subroutine evstu1(dimu,dimx,dimk,
     \   df,nq,dq,nr,dr,nlb,nub,
     \   x,px,s,ps,z,pz,y1,py1,y2,py2,as,az,
     \   sk,yk,err)

      implicit none

      integer dimu,dimx,dimk
      integer nq,nr,nlb,nub,err

      integer i,j,k

      double precision df,dq,dr
      double precision x,px,s,ps,z,pz,y1,py1,y2,py2,as,az
      double precision sk,yk

      double precision d,d2

      dimension df(dimk,dimu),dq(nq,dimu),dr(nr,dimu)
      dimension x(dimu),px(dimu),s(nub+nlb+nr),ps(nub+nlb+nr)
      dimension z(nub+nlb+nr),pz(nub+nlb+nr)
      dimension y1(dimk+dimx),py1(dimk+dimx),y2(nq),py2(nq)
      dimension sk(dimu),yk(dimu)

      err=0

      do i=1,dimk+dimx
         y1(i)=y1(i)+az*py1(i)
         if (isnan(y1(i)).or.abs(y1(i)).gt.huge(0D0)) then
            err=1
            return
         end if
      end do

      do i=1,nq
         y2(i)=y2(i)+az*py2(i)
         if (isnan(y2(i)).or.abs(y2(i)).gt.huge(0D0)) then
            err=1
            return
         end if
      end do

      do i=1,nub+nlb+nr
         s(i)=s(i)+as*ps(i)
         if (isnan(s(i)).or.abs(s(i)).gt.huge(0D0)) then
            err=1
            return
         end if
         z(i)=z(i)+az*pz(i)
         if (isnan(z(i)).or.abs(z(i)).gt.huge(0D0)) then
            err=1
            return
         end if
      end do

      do i=1,dimu
         d2=as*px(i)
         sk(i)=d2
         x(i)=x(i)+d2
         if (isnan(x(i)).or.abs(x(i)).gt.huge(0D0)) then
            err=1
            return
         end if

         d=0D0
         do j=1,dimk
            d=d+df(j,i)*y1(j)
         end do
         do j=1,nq
            d=d+dq(j,i)*y2(j)
         end do
         do j=1,nr
            d=d+dr(j,i)*z(nub+nlb+j)
         end do
         yk(i)=d
      end do

      end


      subroutine evstu2(dimu,dimx,dimk,
     \   df,nq,dq,nr,dr,nlb,nub,
     \   z,y1,y2,
     \   sk,yk,sum)

      implicit none

      integer dimu,dimx,dimk
      integer nq,nr,nlb,nub,err

      integer i,j,k

      double precision df,dq,dr
      double precision z,y1,y2
      double precision sk,yk,sum

      double precision d

      dimension df(dimk,dimu),dq(nq,dimu),dr(nr,dimu)
      dimension z(nub+nlb+nr)
      dimension y1(dimk+dimx),y2(nq)
      dimension sk(dimu),yk(dimu)

      do i=1,dimu
         d=0D0
         do j=1,dimk
            d=d+df(j,i)*y1(j)
         end do
         do j=1,nq
            d=d+dq(j,i)*y2(j)
         end do
         do j=1,nr
            d=d+dr(j,i)*z(nub+nlb+j)
         end do
         yk(i)=yk(i)-d
         sum=sum+yk(i)*sk(i)
      end do

      end



      subroutine evstx1(dimx,dimk,
     \   df,nq,dq,nr,dr,nlb,nub,
     \   x,px,s,ps,z,pz,y1,y2,py2,as,az,
     \   sk,yk,idf,err)

      implicit none

      integer dimx,dimk
      integer nq,nr,nlb,nub,idf,err

      integer i,j,k

      double precision df,dq,dr
      double precision x,px,s,ps,z,pz,y1,y2,py2,as,az
      double precision sk,yk

      double precision d,d2

      dimension df(dimk,dimx),dq(nq,dimx),dr(nr,dimx)
      dimension x(dimx),px(dimx),s(nub+nlb+nr),ps(nub+nlb+nr)
      dimension z(nub+nlb+nr),pz(nub+nlb+nr)
      dimension y1(dimk+dimx),y2(nq),py2(nq)
      dimension sk(dimx),yk(dimx)

      do i=1,nq
         y2(i)=y2(i)+az*py2(i)
         if (isnan(y2(i)).or.abs(y2(i)).gt.huge(0D0)) then
            err=1
            return
         end if
      end do

      do i=1,nub+nlb+nr
         s(i)=s(i)+as*ps(i)
         if (isnan(s(i)).or.abs(s(i)).gt.huge(0D0)) then
            err=1
            return
         end if
         z(i)=z(i)+az*pz(i)
         if (isnan(z(i)).or.abs(z(i)).gt.huge(0D0)) then
            err=1
            return
         end if
      end do

      do i=1,dimx
         d2=as*px(i)
         sk(i)=d2
         x(i)=x(i)+d2
         if (isnan(x(i)).or.abs(x(i)).gt.huge(0D0)) then
            err=1
            return
         end if

         d=0D0
         if (idf.eq.1) then
            do j=1,dimk
               d=d+df(j,i)*y1(j)
            end do
         end if
         do j=1,nq
            d=d+dq(j,i)*y2(j)
         end do
         do j=1,nr
            d=d+dr(j,i)*z(nub+nlb+j)
         end do
         yk(i)=d
      end do

      end


      subroutine evstx2(dimx,dimk,
     \   df,nq,dq,nr,dr,nlb,nub,
     \   z,y1,y2,
     \   sk,yk,idf,sum)

      implicit none

      integer dimx,dimk
      integer nq,nr,nlb,nub,idf

      integer i,j,k

      double precision df,dq,dr
      double precision z,y1,y2
      double precision sk,yk,sum

      double precision d,d2

      dimension df(dimk,dimx),dq(nq,dimx),dr(nr,dimx)
      dimension z(nub+nlb+nr)
      dimension y1(dimk+dimx),y2(nq)
      dimension sk(dimx),yk(dimx)

      do i=1,dimx
         d=0D0
         if (idf.ne.0) then
            do j=1,dimk
               d=d+df(j,i)*y1(j)
            end do
         end if
         do j=1,nq
            d=d+dq(j,i)*y2(j)
         end do
         do j=1,nr
            d=d+dr(j,i)*z(nub+nlb+j)
         end do
         yk(i)=yk(i)-d
         sum=sum+yk(i)*sk(i)
      end do

      end


      subroutine evstk1(dimx,dimk,stages,
     \   df,
     \   x,px,y,as,
     \   sk,yk,err)

      implicit none

      integer dimx,dimk,stages,err

      integer i,j,k

      double precision df,rkb,hstep
      double precision x,px,y,as
      double precision sk,yk

      double precision d,d2,d3

      dimension df(dimk,dimx,stages),rkb(stages)
      dimension x(dimx,stages),px(dimx,stages)
      dimension y(dimk+dimx)
      dimension sk(dimx,stages),yk(dimx,stages)

      do k=1,stages
         do i=1,dimx
            d2=as*px(i,k)
            sk(i,k)=d2
            x(i,k)=x(i,k)+d2
            if (isnan(x(i,k)).or.abs(x(i,k)).gt.huge(0D0)) then
               err=1
               return
            end if

            d=0D0
            do j=1,dimk
               d=d+df(j,i,k)*y(j)
            end do
            yk(i,k)=d
         end do
      end do

      end

      subroutine evstk2(dimx,dimk,stages,
     \   df,
     \   y,
     \   sk,yk,sum)

      implicit none

      integer dimx,dimk,stages

      integer i,j,k

      double precision df
      double precision y
      double precision sk,yk,sum

      double precision d

      dimension df(dimk,dimx,stages)
      dimension y(dimk+dimx)
      dimension sk(dimx,stages),yk(dimx,stages)

      do k=1,stages
         do i=1,dimx
            d=0D0
            do j=1,dimk
               d=d+df(j,i,k)*y(j)
            end do
            yk(i,k)=yk(i,k)-d
            sum=sum+yk(i,k)*sk(i,k)
         end do
      end do

      end
