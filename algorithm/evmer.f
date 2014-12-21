      subroutine merstp(xsize,gsize,x,px,s,ps,alpha,xout,sout)

      implicit none

      integer xsize,gsize
      integer i

      double precision x,px,s,ps,alpha,xout,sout

      dimension x(xsize),px(xsize),xout(xsize)
      dimension s(gsize),ps(gsize),sout(gsize)

      do i=1,xsize
         xout(i)=x(i)+alpha*px(i)
      end do
      do i=1,gsize
         sout(i)=s(i)+alpha*ps(i)
      end do

      end


      subroutine evmer(dim,g,h,
     / nlb,nub,plb,pub,nr,dr,nq,dq,
     / px,ps,s,mu,phi0,dphi0)

      implicit none

      integer dim
      integer nlb,nub,plb,pub,nr,nq
      integer i,j,gidx

      double precision g,h,dr,dq,px,ps,s,mu,phi0,dphi0
      double precision d1,d2,d3,d4,logsi

      dimension pub(nub),plb(nlb)
      dimension g(nub+nlb+nr),h(nq)
      dimension dr(nr,dim),dq(nq,dim)
      dimension px(dim),ps(nub+nlb+nr),s(nub+nlb+nr)

      d1=0D0
      d2=0D0
      d3=0D0
      logsi=0D0
      gidx=0
      do i=1,nub
         d1=d1+log(s(i))
         if (g(i).lt.0D0) then
            d2=d2+abs(g(i))
            d3=d3+px(pub(i))
         end if
         logsi=logsi+ps(i)/s(i)
      end do
      do i=1,nlb
         gidx=nub+i
         d1=d1+log(s(gidx))
         if (g(gidx).lt.0D0) then
            d2=d2+abs(g(gidx))
            d3=d3-px(plb(i))
         end if
         logsi=logsi+ps(gidx)/s(gidx)
      end do
      do i=1,nr
         gidx=nub+nlb+i
         d1=d1+log(s(gidx))
         if (g(gidx).lt.0D0) then
            d2=d2+abs(g(gidx))
            do j=1,dim
               d3=d3-dr(i,j)*px(j)
            end do
         end if
         logsi=logsi+ps(gidx)/s(gidx)
      end do

      do i=1,nq
         d2=d2+abs(h(i))
         d4=0D0
         do j=1,dim
            d4=d4+dq(i,j)*px(j)
         end do
         call mysign(d4,h(i),d3)
      end do

      phi0=phi0 - mu*d1 + d2
      dphi0=dphi0 - mu*logsi + d3

      end

      subroutine evme2(dim,g,h,
     / nlb,nub,plb,pub,nr,nq,
     / s,mu,phi0)

      implicit none

      integer dim
      integer nlb,nub,plb,pub,nr,nq
      integer i,j,gidx

      double precision g,h,s,mu,phi0
      double precision d1,d2

      dimension pub(nub),plb(nlb)
      dimension g(nub+nlb+nr),h(nq)
      dimension s(nub+nlb+nr)

      d1=0D0
      d2=0D0
      gidx=0
      do i=1,nub
         d1=d1+log(s(i))
         if (g(i).lt.0D0) then
            d2=d2+abs(g(i))
         end if
      end do
      do i=1,nlb
         gidx=nub+i
         d1=d1+log(s(gidx))
         if (g(gidx).lt.0D0) then
            d2=d2+abs(g(gidx))
         end if
      end do
      do i=1,nr
         gidx=nub+nlb+i
         d1=d1+log(s(gidx))
         if (g(gidx).lt.0D0) then
            d2=d2+abs(g(gidx))
         end if
      end do

      do i=1,nq
         d2=d2+abs(h(i))
      end do

      phi0=phi0 - mu*d1 + d2

      end


      subroutine evmeru(dimu,dimx,dimk,stages,h,px,dfu,dphi0)

      implicit none

      integer dimu,dimx,dimk,stages
      integer i,j

      double precision h,px,dfu,dphi0
      double precision d1,d2

      dimension px(dimu),h(dimk+dimx),dfu(dimk,dimu)

      d1=0D0
      do i=1,dimk
         d2=0D0
         do j=1,dimu
            d2=d2+dfu(i,j)*px(j)
         end do
         call mysign(d2,h(i),d1)
      end do
      dphi0 = dphi0+d1

      end


      subroutine evmerx(dimx,dimk,stages,h,px,pxA,dfx,dfxA,
     \ rkb,hstep,phi0,dphi0)

      implicit none

      integer dimx,dimk,stages
      integer i,j,k

      double precision h,px,pxA,dfx,dfxA,phi0,dphi0,rkb,hstep
      double precision d1,d2,d3,d4

      dimension px(2*dimx),pxA(dimx,stages),h(dimk+dimx)
      dimension dfx(dimk,dimx),dfxA(dimk,dimx,stages),rkb(stages)

      d1=0D0
      d4=0D0
      do i=1,dimx
         call mysign(px(i)-px(i+dimx), h(dimk+i), d1)
         d4=d4+abs(h(dimk+i))
      end do

      do i=1,dimk
         d2=0D0
         do j=1,dimx
            d2=d2+dfx(i,j)*px(j)
         end do
         do k=1,stages
            do j=1,dimx
               d2=d2+dfxA(i,j,k)*pxA(j,k)
            end do
         end do
         call mysign(d2,h(i),d1)
         d4=d4+abs(h(i))
      end do

      do k=1,stages
         d3=hstep*rkb(k)
         do i=1,dimx
            call mysign(pxA(i,k)*d3, h(dimk+i), d1)
         end do
      end do
      dphi0 = dphi0+d1
      phi0 = phi0+d4
      end

      subroutine evmex2(dimx,dimk,h,phi0)

      implicit none

      integer dimx,dimk
      integer i

      double precision h,phi0
      double precision d4

      dimension h(dimk+dimx)

      d4=0D0
      do i=1,dimx
         d4=d4+abs(h(dimk+i))
      end do

      do i=1,dimk
         d4=d4+abs(h(i))
      end do

      phi0 = phi0+d4
      end


      subroutine mysign(x,y,res)

      implicit none

      double precision x,y,res

      if (y.ne.0D0) then
C      if (abs(y).gt.1e-14) then
         if (y.lt.0D0) then
            res=res-x
         else
            res=res+x
         end if
      else
         res=0D0
      end if

      end

      subroutine stepsz(xsize,gsize,hsize,as,az,dphi0,
     \ px,ps,pz,py,mtau,s,z)

      implicit none

      integer xsize,gsize,hsize
      integer i

      double precision as,az,dphi0,px,ps,pz,py,mtau,s,z
      double precision dt

      dimension px(xsize),ps(gsize),pz(gsize),py(hsize)
      dimension s(gsize),z(gsize)

      dt=1D0-mtau

      if (dphi0.gt.0D0) then
         dphi0=-dphi0
         do i=1,gsize
            ps(i)=-ps(i)
            pz(i)=-pz(i)
         end do
         do i=1,hsize
            py(i)=-py(i)
         end do
         do i=1,xsize
            px(i)=-px(i)
         end do
      end if

      as=1D0
      az=1D0

      do i=1,gsize
         if (ps(i).lt.0D0) then
            as=min(-dt*s(i)/ps(i),as)
         end if
         if (pz(i).lt.0D0) then
            az=min(-dt*z(i)/pz(i),az)
         end if
      end do

      end

      subroutine lnsrch(phi0,phia,dphi0,alpha,eps,res)

      implicit none

      integer res
      double precision phi0,dphi0,phia,alpha,eps
      double precision d

      d=alpha*dphi0
      if (phia.gt.phi0+1D-4*d) then
         alpha=min(alpha,(-d*alpha)/(2D0*(phia-phi0-d)))
         if (alpha.lt.eps*1D-2) then
            res=-1
         else
            res=1
         end if
      else
         res=0
      end if

      end





