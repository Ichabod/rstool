      subroutine evcx(xsize,hsize,yb,sb,xi,rhs,ocx1,ocx2,err)

      implicit none

      double precision yb,sb,xi,rhs,ocx1,ocx2
      integer xsize,hsize,err

      integer i,j
      double precision t11,t12,t21,t22,d1
      double precision mb1,mb2,cx1,cx2

      dimension rhs(xsize+hsize,3),yb(xsize),sb(xsize)

      err=0
      mb1=0D0
      mb2=0D0
      do i=1,xsize
         mb1=mb1 + sb(i)*sb(i)*xi
         mb2=mb2 - sb(i)*yb(i)
      end do

      t11=1D0
      t21=0D0
      t12=0D0
      t22=1D0
      cx1=0D0
      cx2=0D0

      do i=1,xsize
         d1=sb(i)*xi/mb1
         t11=t11+d1*rhs(i,2)
         t12=t12+d1*rhs(i,3)
         cx1=cx1+d1*rhs(i,1)

         d1=yb(i)/mb2
         t21=t21+d1*rhs(i,2)
         t22=t22+d1*rhs(i,3)
         cx2=cx2+d1*rhs(i,1)
      end do

      if (abs(t11).gt.abs(t21)) then
         if (abs(t11)<1D-14) then
            err=100
            return
         end if
         d1=t21/t11
         t22=t22-t12*d1
         cx2=(cx2-cx1*d1)/t22
         cx1=(cx1-t12*cx2)/t11
         ocx1=cx1
         ocx2=cx2
      else
         if (abs(t21)<1D-14) then
            err=100
            return
         end if
         d1=t11/t21
         t12=t12-t22*d1
         cx1=(cx1-cx2*d1)/t12
         cx2=(cx2-t22*cx1)/t21
         ocx1=cx2
         ocx2=cx1
      end if


      end


      subroutine evpxsz(dim,g,
     / nlb,nub,plb,pub,nr,dr,
     / px,ps,pz,s,z,mu,cx1,cx2,ax1,ax2,hasB)

      implicit none

      integer dim,nlb,nub,nr,plb,pub,hasB
      integer i,j,idx

      double precision g,dr
      double precision px,ps,pz,s,z,mu,cx1,cx2,ax1,ax2
      double precision d1,d2

      dimension g(nlb+nub+nr),dr(nr,dim)
      dimension px(dim),ps(nlb+nub+nr),pz(nlb+nub+nr)
      dimension s(nlb+nub+nr),z(nlb+nub+nr)
      dimension ax1(dim),ax2(dim)
      dimension plb(nlb),pub(nub)

      idx=0
      do i=1,nub
         ps(i)=-px(pub(i)) - s(i) + g(i)
      end do
      idx=idx+nub
      do i=1,nlb
         ps(idx+i)=px(plb(i)) - s(idx+i) + g(idx+i)
      end do
      idx=idx+nlb
      do i=1,nr
         d1=0D0
         do j=1,dim
            d1=d1+dr(i,j)*px(j)
         end do
         ps(idx+i)=d1 - s(idx+i) + g(idx+i)
      end do
      do i=1,nub+nlb+nr
         pz(i)=(-z(i)*ps(i) + mu)/s(i) - z(i)
      end do

      if (hasB.ne.0) then
         idx=0
         do i=1,nub
            d1=-ax1(pub(i))*cx1-ax2(pub(i))*cx2
            ps(i)=ps(i)-d1
            pz(i)=pz(i)+z(i)/s(i)*d1
         end do
         idx=idx+nub
         do i=1,nlb
            d1=+ax1(plb(i))*cx1+ax2(plb(i))*cx2
            ps(idx+i)=ps(idx+i)-d1
            pz(idx+i)=pz(idx+i)+z(idx+i)/s(idx+i)*d1
         end do
         idx=idx+nlb
         do i=1,nr
            d1=0D0
            d2=0D0
            do j=1,dim
               d1=d1+dr(i,j)*ax1(j)
               d2=d2+dr(i,j)*ax2(j)
            end do
            d1=cx1*d1+cx2*d2
            ps(idx+i)=ps(idx+i) - d1
            pz(idx+i)=pz(idx+i) + z(idx+i)/s(idx+i)*d1
         end do
         do i=1,dim
            px(i)=px(i)-ax1(i)*cx1-ax2(i)*cx2
         end do
      end if

      end


      subroutine evpxszk(dim,px,cx1,cx2,ax1,ax2)

      implicit none

      integer dim
      integer i

      double precision px,cx1,cx2,ax1,ax2

      dimension px(dim)
      dimension ax1(dim),ax2(dim)

      do i=1,dim
         px(i)=px(i)-ax1(i)*cx1-ax2(i)*cx2
      end do

      end


      subroutine evpy(hsize,py,cx1,cx2,ax1,ax2,hasB)

      implicit none

      integer hsize,hasB
      integer i

      double precision py,cx1,cx2,ax1,ax2

      dimension py(hsize)
      dimension ax1(hsize),ax2(hsize)

      do i=1,hsize
         py(i)=-py(i)
      end do

      if (hasB.ne.0) then
         do i=1,hsize
            py(i)=py(i)+ax1(i)*cx1+ax2(i)*cx2
         end do
      end if

      end
