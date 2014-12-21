      subroutine evrhs(rhs,dimx,dimu,dimk,xsize,hsize,N,
     \ nub_u,nlb_u,nub_x,nlb_x,
     \ pub_u,plb_u,pub_x,plb_x,
     \ nr_u,nr_x,r_u,r_x,
     \ nq_u,nq_x,q_u,q_x,
     \ g,h,dfu,dfx,dfxA,rks,rkb,hstep,y,z,s,mu,
     \ yb,sb,hasB,xi)

      implicit none

      integer dimx,dimu,dimk,N,rks,xsize,hsize
      integer nq_u,nq_x
      integer nub_u,nlb_u,nub_x,nlb_x
      integer pub_u,plb_u,pub_x,plb_x
      integer nr_u,nr_x
      integer hasB
      double precision rhs,g,h,dfx,dfu,dfxA
      double precision rkb,hstep,y,z,s,mu
      double precision q_u,q_x,r_u,r_x
      double precision yb,sb,xi

      dimension rhs(*),g(*),h(*),rkb(*)
      dimension dfx(*),dfu(*),dfxA(*)
      dimension y(*),z(*),s(*)
      dimension nq_u(*),nq_x(*),q_u(*),q_x(*)
      dimension nub_u(*), nlb_u(*), nub_x(*), nlb_x(*)
      dimension pub_u(*), plb_u(*), pub_x(*), plb_x(*)
      dimension nr_u(*), nr_x(*),r_u(*),r_x(*)
      dimension yb(*),sb(*)

      integer irhsu,irhsx,irhsk,irhsh
      integer iy1,iy2
      
      irhsu=1
      irhsx=irhsu+N*dimu
      irhsk=irhsx+(N+1)*dimx
      irhsh=irhsk+N*dimk

      iy1=1
      iy2=iy1+N*(dimk+dimx)

      call evrhco(rhs(irhsu),rhs(irhsx),rhs(irhsk),rhs(irhsh),
     \ dimx,dimu,dimk,hsize,N,
     \ nub_u,nlb_u,nub_x,nlb_x,
     \ pub_u,plb_u,pub_x,plb_x,
     \ nr_u,nr_x,r_u,r_x,
     \ nq_u,nq_x,q_u,q_x,
     \ g,h,dfu,dfx,dfxA,rks,rkb,hstep,
     \ y(iy1),y(iy2),z,s,mu)

      if (hasB.eq.1) then
         call evbfco(rhs(xsize+hsize+1),xsize,hsize,yb,sb,xi)
      end if



      end

      subroutine evbfco(rhs,xsize,hsize,yb,sb,xi)

      implicit none

      double precision rhs,yb,sb,xi
      integer xsize,hsize,i

      dimension rhs(xsize+hsize,2),yb(xsize),sb(xsize)

      rhs=0D0

      do i=1,xsize
         rhs(i,1)=xi*sb(i)
         rhs(i,2)=yb(i)
      end do

      end


      subroutine evrhco(rhsu,rhsx,rhsk,rhsh,
     \ dimx,dimu,dimk,hsize,N,
     \ nub_u,nlb_u,nub_x,nlb_x,
     \ pub_u,plb_u,pub_x,plb_x,
     \ nr_u,nr_x,r_u,r_x,
     \ nq_u,nq_x,q_u,q_x,
     \ g,h,dfu,dfx,dfxA,rks,rkb,hstep,
     \ y1,y2,z,s,mu)

      implicit none

      integer dimx,dimu,dimk,N,rks,hsize
      integer nq_u,nq_x
      integer nub_u,nlb_u,nub_x,nlb_x
      integer pub_u,plb_u,pub_x,plb_x
      integer nr_u,nr_x

      double precision rhsu,rhsx,rhsk,rhsh,g,h,dfx,dfu,dfxA
      double precision rkb,hstep,y1,y2,z,s,mu
      double precision q_u,q_x,r_u,r_x
      double precision d1,d2

      integer i,j,m,l,iqu,iqx,iy,iub,ilb,ig,ir

      dimension rhsu(dimu,N),rhsx(dimx,N+1),rhsk(dimx,rks,N)
      dimension rhsh(hsize),g(*),h(hsize),rkb(*)
      dimension dfx(dimk,dimx,N),dfu(dimk,dimu,N)
      dimension dfxA(dimk,dimx,rks,N)
      dimension y1(dimk+dimx,N),y2(*),z(*),s(*)
      dimension nq_u(N),nq_x(N+1)
      dimension q_x(*),q_u(*)
      dimension nub_u(N), nlb_u(N), nub_x(N+1), nlb_x(N+1)
      dimension pub_u(*), plb_u(*), pub_x(*), plb_x(*)
      dimension nr_u(N), nr_x(N+1),r_u(*),r_x(*)

      do i=1,hsize
         rhsh(i)=-h(i)
      end do

      iqu=1
      iqx=1

      iy=0
      iub=1
      ilb=1
      ig=1
      ir=1
      do i=1,N
         do j=1,dimu
            d1 = 0D0
            do m=1,dimk
               d1=d1+dfu(m,j,i)*y1(m,i)
            end do
            if (nq_u(i).gt.0) then
               do m=1,nq_u(i)
                  d1=d1+q_u(iqu)*y2(iy+m)
                  iqu=iqu+1
               end do
            end if
            rhsu(j,i) = d1
         end do
         iy=iy+nq_u(i)
         if (nub_u(i).gt.0) then
            do j=1,nub_u(i)
               rhsu(pub_u(iub),i)=rhsu(pub_u(iub),i)-
     \           z(ig)+z(ig)/s(ig)*g(ig) - mu/s(ig)
               iub=iub+1
               ig=ig+1
            end do
         end if
         if (nlb_u(i).gt.0) then
            do j=1,nlb_u(i)
               rhsu(plb_u(ilb),i)=rhsu(plb_u(ilb),i)+
     \           z(ig)-z(ig)/s(ig)*g(ig) + mu/s(ig)
               ilb=ilb+1
               ig=ig+1
            end do
         end if
         if (nr_u(i).gt.0) then
            do m=1,nr_u(i)
              d1 = z(ig)-z(ig)/s(ig)*g(ig) + mu/s(ig)
              ig=ig+1
              do j=1,dimu
                 rhsu(j,i)=rhsu(j,i)+d1*r_u(ir)
                 ir=ir+1
              end do
            end do
         end if
      end do

      do j=1,dimx
         rhsx(j,1)=0D0
      end do

      ir=1
      iub=1
      ilb=1
      do i=1,N+1
         do j=1,dimx
            d1 = 0D0
            if (i.le.N) then            
               rhsx(j,i+1)=-y1(dimk+j,i)
               do m=1,dimk
                  d1=d1+dfx(m,j,i)*y1(m,i)
               end do
               d1=d1+y1(dimk+j,i)
            end if
            if (nq_x(i).gt.0) then
               do m=1,nq_x(i)
                  d1=d1+q_x(iqx)*y2(iy+m)
                  iqx=iqx+1
               end do
            end if
            rhsx(j,i) = rhsx(j,i) + d1
         end do
         iy=iy+nq_x(i)
         if (nub_x(i).gt.0) then
            do j=1,nub_x(i)
               rhsx(pub_x(iub),i)=rhsx(pub_x(iub),i)-
     \           z(ig)+z(ig)/s(ig)*g(ig) - mu/s(ig)
               iub=iub+1
               ig=ig+1
            end do
         end if
         if (nlb_x(i).gt.0) then
            do j=1,nlb_x(i)
               rhsx(plb_x(ilb),i)=rhsx(plb_x(ilb),i)+
     \           z(ig)-z(ig)/s(ig)*g(ig) + mu/s(ig)
               ilb=ilb+1
               ig=ig+1
            end do
         end if
         if (nr_x(i).gt.0) then
            do m=1,nr_x(i)
              d1 = z(ig)-z(ig)/s(ig)*g(ig) + mu/s(ig)
              ig=ig+1
              do j=1,dimx
                 rhsx(j,i)=rhsx(j,i)+d1*r_x(ir)
                 ir=ir+1
              end do
            end do
         end if
      end do

      do i=1,N
         do l=1,rks
            d2=rkb(l)*hstep
            do j=1,dimx
               d1 = d2*y1(dimk+j,i)
               do m=1,dimk
                  d1=d1+dfxA(m,j,l,i)*y1(m,i)
               end do
               rhsk(j,l,i) = d1
            end do
         end do
      end do

      end


