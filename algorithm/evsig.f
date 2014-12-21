      subroutine evsig(xi,vec_s,vec_z,N,dimx,dimu,gsize,
     \   nub_u,nlb_u,nub_x,nlb_x,
     \   pub_u,plb_u,pub_x,plb_x,
     \   nr_u,nr_x,
     \   iub_u,ilb_u,ir_u,iub_x,ilb_x,ir_x,
     \   isigu,isigx,sigu,sigx,
     \   dru,drx)

      implicit none
      integer i,N,dimx,dimu,gsize
      integer nub_u,nlb_u,nub_x,nlb_x
      integer pub_u,plb_u,pub_x,plb_x
      integer nr_u,nr_x
      integer iub_u,ilb_u,ir_u,iub_x,ilb_x,ir_x
      integer isigu,isigx,iub,ilb,ir
      double precision xi,vec_s,vec_z
      double precision dru,drx
      double precision sigu,sigx

      dimension nub_u(*),nlb_u(*),nub_x(*),nlb_x(*)
      dimension pub_u(*),plb_u(*),pub_x(*),plb_x(*)
      dimension iub_u(*),ilb_u(*),ir_u(*),iub_x(*),ilb_x(*),ir_x(*)
      dimension nr_u(*),nr_x(*)
      dimension dru(*),drx(*)
      dimension vec_s(*),vec_z(*)
      dimension isigu(*),isigx(*),sigu(*),sigx(*)

      iub=1
      ilb=1
      ir=1
      do i=1,N
         if (nr_u(i).eq.0) then
            call sigdia(sigu(isigu(i)),xi,dimu,
     \ nub_u(i), nlb_u(i), pub_u(iub), plb_u(ilb), 
     \ vec_s(iub_u(i)+1), vec_s(ilb_u(i)+1), 
     \ vec_z(iub_u(i)+1), vec_z(ilb_u(i)+1))
         else
            call sigsqr(sigu(isigu(i)),xi,dimu,
     \ nub_u(i), nlb_u(i), pub_u(iub), plb_u(ilb), 
     \ vec_s(iub_u(i)+1), vec_s(ilb_u(i)+1), 
     \ vec_z(iub_u(i)+1), vec_z(ilb_u(i)+1),
     \ nr_u(i),dru(ir), vec_s(ir_u(i)+1), vec_z(ir_u(i)+1));
            ir=ir+nr_u(i)*dimu
         end if
         iub=iub+nub_u(i)
         ilb=ilb+nlb_u(i)
      end do

      iub=1
      ilb=1
      ir=1
      do i=1,N+1
         if (nr_x(i).eq.0) then
            call sigdia(sigx(isigx(i)),xi,dimx,
     \ nub_x(i), nlb_x(i), pub_x(iub), plb_x(ilb), 
     \ vec_s(iub_x(i)+1), vec_s(ilb_x(i)+1), 
     \ vec_z(iub_x(i)+1), vec_z(ilb_x(i)+1))
         else
            call sigsqr(sigx(isigx(i)),xi,dimx,
     \ nub_x(i), nlb_x(i), pub_x(iub), plb_x(ilb), 
     \ vec_s(iub_x(i)+1), vec_s(ilb_x(i)+1), 
     \ vec_z(iub_x(i)+1), vec_z(ilb_x(i)+1),
     \ nr_x(i),drx(ir), vec_s(ir_x(i)+1), vec_z(ir_x(i)+1));
            ir=ir+nr_x(i)*dimx
         end if
         iub=iub+nub_x(i)
         ilb=ilb+nlb_x(i)
      end do

      end


      subroutine sigdia(sig,xi,dim,
     \ nub,nlb,pub,plb,s_ub,s_lb,z_ub,z_lb)

      implicit none

      integer dim,i
      integer nub,nlb,pub,plb
      double precision sig,xi,s_ub,s_lb,z_ub,z_lb

      dimension sig(dim),pub(nub),plb(nlb)
      dimension s_ub(nub),s_lb(nlb),z_ub(nub),z_lb(nlb)

      sig=xi

      do i=1,nub
         sig(pub(i))=sig(pub(i))+z_ub(i)/s_ub(i)
      end do
      do i=1,nlb
         sig(plb(i))=sig(plb(i))+z_lb(i)/s_lb(i)
      end do

      end


      subroutine sigsqr(sig,xi,dim,
     \ nub,nlb,pub,plb,s_ub,s_lb,z_ub,z_lb,
     \ nr,dr,s_r,z_r)

      implicit none

      integer dim,i,j,k
      integer nub,nlb,pub,plb,nr
      double precision sig,xi,s_ub,s_lb,z_ub,z_lb
      double precision dr,s_r,z_r
      double precision temp

      dimension sig(dim,dim),pub(nub),plb(nlb)
      dimension s_ub(nub),s_lb(nlb),z_ub(nub),z_lb(nlb)
      dimension dr(nr,dim),s_r(nr),z_r(nr)

      sig=0D0

      do i=1,dim
         sig(i,i)=xi
      end do
      do i=1,nub
         j=pub(i)
         sig(j,j)=sig(j,j)+z_ub(i)/s_ub(i)
      end do
      do i=1,nlb
         j=plb(i)
         sig(j,j)=sig(j,j)+z_lb(i)/s_lb(i)
      end do
      do i=1,dim
         do j=1,dim
            temp=0D0
            do k=1,nr
              temp=temp+dr(k,i)*z_r(k)/s_r(k)*dr(k,j)
            end do
            sig(i,j)=sig(i,j)+temp
         end do
      end do

      end
