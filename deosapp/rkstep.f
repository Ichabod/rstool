C     work: s*dim*s*dim + s*dim + s*dim + dim = (s*s*dim + s + s + 1)*dim = (s*(s*dim + 2) + 1)*dim
C     iwork: s*dim

      subroutine rkstep(dim,h,t,x,xout,u,p,s,k,rkA,rkb,rkc,eps,
     \ work,iwork,ierr)

      implicit none

      integer dim,s,iwork,ierr,i,j
      double precision d
      double precision t,h,x,xout,u,p,k,rkA,rkb,rkc,work,eps
      integer tdim

      integer ifres,ifres2,idf,irktmp

      dimension x(dim),xout(dim),u(*),k(dim,s),rkA(s,s),rkb(s),rkc(s)
      dimension iwork(*),work(*),p(*)

      tdim=s*dim
      idf=1
      ifres=idf+tdim*tdim
      ifres2=ifres+tdim
      irktmp=ifres2+tdim

      call newslv(dim,h,t,x,u,p,k,s,eps,rkA,rkc,
     \  work(idf),work(ifres),work(ifres2),iwork,work(irktmp),ierr)

      if (ierr.ne.0) then
         return
      end if

      do i=1,dim
         d=0D0
         do j=1,s
            d=d+rkb(j)*k(i,j)
         end do
         xout(i)=x(i)+h*d
      end do


      end



      subroutine newslv(dim,h,t,x,u,p,k,s,eps,rkA,rkc,
     \ df,fres,fres2,ip,rktmp,ier)

      implicit none

      integer dim,s,i,j,tdim,ip,ier,cnt
      double precision t,h,x,u,p,k,rkA,rkc,eps
      double precision fres,fres2,rktmp
      double precision df,dx,d

      dimension x(dim),u(*),k(*),p(*),rkA(*),rkc(*),ip(*)
      dimension fres(dim*s),fres2(dim*s),df(s*dim,s*dim),rktmp(*)

      external dec,sol

      tdim=dim*s
      cnt=0
      eps=eps**2

      do
C        Compute Jacobian      
         call evrkf(dim,h,t,x,u,p,k,s,rkA,rkc,fres,rktmp)
         do i=1,tdim
            d=k(i)
            k(i)=k(i)+1D-6
            call evrkf(dim,h,t,x,u,p,k,s,rkA,rkc,fres2,rktmp)
            k(i)=d
            do j=1,tdim
               df(j,i)=-(fres2(j)-fres(j))/1D-6
            end do
         end do

C        Solve LGS
         call dec(tdim,tdim,df,ip,ier)
         if (ier.ne.0) then
            return
         end if
         call sol(tdim,tdim,df,fres,ip)

C        Newton Step
         d=0D0
         do i=1,tdim
            k(i) = k(i) + fres(i)
            d=d+fres(i)**2
         end do
         cnt=cnt+1
         if (d.lt.eps) then
            ier=0
            return
         end if

         if (cnt.ge.100) then
            ier=100
            return
         end if

      end do

      end


      subroutine evrkf(dim,h,t,x,u,p,k,s,rkA,rkc,res,tmp)

      implicit none

      integer dim,i,j,j2,s
      double precision t,h,x,u,p,k,rkA,rkc,res,tmp,d,d2,d3

      external rhs

      dimension x(dim),u(*),k(dim,s),p(*),rkA(s,s),rkc(s)
      dimension res(dim,s),tmp(dim)

      do i=1,s
         d=t+h*rkc(i)
         tmp=0D0
         do j=1,dim
            d2=0D0
            do j2=1,s
               d2=d2+rkA(i,j2)*k(j,j2)
            end do
            tmp(j)=x(j)+h*d2
         end do
         call rhs(d,tmp,u,p,res(1,i))
         do j=1,dim
            res(j,i)=res(j,i)-k(j,i)
         end do
      end do

      end

      subroutine evrkfu(dim,dimu,h,t,x,u,p,k,s,rkA,rkc,res,tmp,tmp2)

      implicit none

      integer dim,dimu,i,j,j2,s
      double precision t,h,x,u,p,k,rkA,rkc,res,tmp,tmp2,d,d2,d3

      external rhs

      dimension x(dim),u(*),k(dim,s),p(*),rkA(s,s),rkc(s)
      dimension res(dim,s,dimu),tmp(dim),tmp2(dim,dimu)

      do i=1,s
         d=t+h*rkc(i)
         tmp=0D0
         do j=1,dim
            d2=0D0
            do j2=1,s
               d2=d2+rkA(i,j2)*k(j,j2)
            end do
            tmp(j)=x(j)+h*d2
         end do
         call rhsu(d,tmp,u,p,tmp2)
         res(1:dim,i,1:dimu) = tmp2(1:dim,1:dimu)         
      end do

      end

      subroutine evrkfx(dim,h,t,x,u,p,k,s,rkA,rkc,res,tmp,tmp2)

      implicit none

      integer dim,i,j,j2,s
      double precision t,h,x,u,p,k,rkA,rkc,res,tmp,tmp2,d,d2,d3

      external rhs

      dimension x(dim),u(*),k(dim,s),p(*),rkA(s,s),rkc(s)
      dimension res(dim,s,dim),tmp(dim),tmp2(dim,dim)

      do i=1,s
         d=t+h*rkc(i)
         tmp=0D0
         do j=1,dim
            d2=0D0
            do j2=1,s
               d2=d2+rkA(i,j2)*k(j,j2)
            end do
            tmp(j)=x(j)+h*d2
         end do
         call rhsx(d,tmp,u,p,tmp2)
         res(1:dim,i,1:dim) = tmp2(1:dim,1:dim)         
      end do

      end