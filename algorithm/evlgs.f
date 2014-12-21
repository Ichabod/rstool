      subroutine evbu(dimu,dimx,stages,sigu,isq,dfu,nrhs,
     \ mrhs,rhs,bidx,y13,err,pinfo)

      implicit none

      integer isq,nrhs,mrhs,bidx,dimu,dimx,stages
      integer i,j,ys,err,pinfo
      double precision sigu,dfu,rhs,y13

      dimension sigu(*),dfu(dimx*stages,dimu)
      dimension y13(dimu,dimx*stages+mrhs)
      dimension rhs(nrhs,mrhs),pinfo(*)

      ys=dimx*stages
      err=0

      if (isq.eq.0) then
         do i=1,dimu
            do j=1,ys
               y13(i,j)=dfu(j,i)/sigu(i)
            end do
            do j=1,mrhs
               y13(i,ys+j)=rhs(bidx+i,j)/sigu(i)
            end do
         end do
      else
         do i=1,dimu
            do j=1,ys
               y13(i,j)=dfu(j,i)
            end do
            do j=1,mrhs
               y13(i,ys+j)=rhs(bidx+i,j)
            end do
         end do
C         call choslv(dimu,ys+mrhs,dimu,0,sigu,y13,err)
         call chlsco(dimu,ys+mrhs,dimu,0,sigu,y13,err,pinfo)
      end if

      end



      subroutine evbx(dimx,stages,sigx,isq,dfx,nrhs,
     \ mrhs,rhs,bidx,y579,err,pinfo)

      implicit none

      integer isq,nrhs,mrhs,bidx,dimx,stages
      integer i,j,ys1,ys2,err,pinfo
      double precision sigx,dfx,rhs,y579

      dimension sigx(*),dfx(dimx*stages,dimx)
      dimension y579(dimx,dimx*stages + dimx +mrhs)
      dimension rhs(nrhs,mrhs),pinfo(*)

      ys1=dimx*stages
      ys2=ys1+dimx
      y579=0
      err=0

      if (isq.eq.0) then
         do i=1,dimx
            do j=1,ys1
               y579(i,j)=dfx(j,i)/sigx(i)
            end do
            y579(i,ys1+i)=1D0/sigx(i)
            do j=1,mrhs
               y579(i,ys2+j)=rhs(bidx+i,j)/sigx(i)
            end do
         end do
      else
         do i=1,dimx
            do j=1,ys1
               y579(i,j)=dfx(j,i)
            end do
            y579(i,ys1+i)=1D0
            do j=1,mrhs
               y579(i,ys2+j)=rhs(bidx+i,j)
            end do
         end do
C         call choslv(dimx,ys2+mrhs,dimx,0,sigx,y579,err)
         call chlsco(dimx,ys2+mrhs,dimx,0,sigx,y579,err,pinfo)
      end if

      end

      subroutine evbxl(dimx,sigx,isq,nrhs,
     \ mrhs,rhs,bidx,Q,R,DN,cNidx,err,pinfo)

      implicit none

      integer isq,nrhs,mrhs,bidx,dimx
      integer i,j,err,cNidx,pinfo
      double precision sigx,rhs,Q,R,DN

      dimension sigx(*)
      dimension Q(dimx,dimx),R(dimx,mrhs)
      dimension rhs(nrhs,mrhs),DN(dimx,dimx),pinfo(*)

      err=0
      Q=0

      if (isq.eq.0) then
         do i=1,dimx
            Q(i,i)=-1D0/sigx(i)
            do j=1,mrhs
               R(i,j)=rhs(bidx+i,j)/sigx(i)
            end do
         end do
      else
         do i=1,dimx
            Q(i,i)=-1D0
            do j=1,mrhs
               R(i,j)=rhs(bidx+i,j)
            end do
         end do
C         call choslv(dimx,dimx,dimx,0,sigx,Q,err)
         call chlsco(dimx,dimx,dimx,0,sigx,Q,err,pinfo)
         if (err.eq.0) then
C            call rrslv(dimx,mrhs,dimx,0,sigx,R)
            call rrslvo(dimx,mrhs,dimx,0,sigx,R,pinfo)
         end if
      end if

      do i=1,dimx
         do j=1,dimx
            DN(i,j)=DN(i,j)-Q(i,j)
         end do
         do j=1,mrhs
            rhs(cNidx+i,j)=rhs(cNidx+i,j)-R(i,j)
         end do
      end do

      end



      subroutine evbur(dimu,dimx,stages,sigu,isq,dfu,nrhs,
     \ mrhs,rhs,dqu,nqu,bidx,y13,y24,D,temp,err,pinfo)

      implicit none

      integer isq,nrhs,mrhs,bidx,dimu,dimx,stages,nqu
      integer i,j,k,ys,err,pinfo
      double precision sigu,dfu,rhs,y13,y24,dqu,D,temp,val

      dimension sigu(*),dfu(dimx*stages,dimu),dqu(nqu,dimu)
      dimension y13(dimu,dimx*stages+mrhs)
      dimension y24(nqu,dimx*stages+mrhs)
      dimension rhs(nrhs,mrhs),D(nqu,nqu),temp(dimu,nqu)
      dimension pinfo(*)

      ys=dimx*stages
      err=0

      if (isq.eq.0) then
         do i=1,dimu
            do j=1,ys
               y13(i,j)=dfu(j,i)/sigu(i)
            end do
            do j=1,mrhs
               y13(i,ys+j)=rhs(bidx+i,j)/sigu(i)
            end do
            do j=1,nqu
               temp(i,j)=dqu(j,i)/sigu(i)
            end do
         end do
      else
         do i=1,dimu
            do j=1,ys
               y13(i,j)=dfu(j,i)
            end do
            do j=1,mrhs
               y13(i,ys+j)=rhs(bidx+i,j)
            end do
            do j=1,nqu
               temp(i,j)=dqu(j,i)
            end do
         end do
C         call choslv(dimu,ys+mrhs,dimu,0,sigu,y13,err)
         call chlsco(dimu,ys+mrhs,dimu,0,sigu,y13,err,pinfo)
         if (err.ne.0) then
            return
         end if
C         call rrslv(dimu,nqu,dimu,0,sigu,temp)
         call rrslvo(dimu,nqu,dimu,0,sigu,temp,pinfo)
      end if

      do i=1,nqu
         do j=1,ys+mrhs
            val=0
            do k=1,dimu
               val=val+dqu(i,k)*y13(k,j)
            end do
            y24(i,j)=val
         end do
         do j=1,mrhs
            y24(i,ys+j)=y24(i,ys+j)-rhs(bidx+dimu+i,j)
         end do
      end do
      do i=1,nqu
         do j=1,nqu
            val=0
            do k=1,dimu
               val=val+dqu(i,k)*temp(k,j)
            end do
            D(i,j)=val
         end do
      end do

      call choslv(nqu,ys+mrhs,nqu,0,D,y24,err)
      if (err.ne.0) then
         return
      end if

      do i=1,dimu
         do j=1,ys+mrhs
            val=0
            do k=1,nqu
               val=val+temp(i,k)*y24(k,j)
            end do
            y13(i,j)=y13(i,j)-val;
         end do
      end do

      end


      subroutine evbxr(dimx,stages,sigx,isq,dfx,nrhs,
     \ mrhs,rhs,dqx,nqx,bidx,y579,y6810,D,temp,err,pinfo)

      implicit none

      integer isq,nrhs,mrhs,bidx,dimx,stages,nqx
      integer i,j,ys1,ys2,err,k,pinfo
      double precision sigx,dfx,rhs,y579,y6810,dqx,D,temp,val

      dimension sigx(*),dfx(dimx*stages,dimx),dqx(nqx,dimx)
      dimension y579(dimx,dimx*stages + dimx +mrhs)
      dimension y6810(nqx,dimx*stages + dimx +mrhs)
      dimension rhs(nrhs,mrhs),D(nqx,nqx),temp(dimx,nqx)
      dimension pinfo(*)

      ys1=dimx*stages
      ys2=ys1+dimx
      y579=0
      err=0

      if (isq.eq.0) then
         do i=1,dimx
            do j=1,ys1
               y579(i,j)=dfx(j,i)/sigx(i)
            end do
            y579(i,ys1+i)=1D0/sigx(i)
            do j=1,mrhs
               y579(i,ys2+j)=rhs(bidx+i,j)/sigx(i)
            end do
            do j=1,nqx
               temp(i,j)=dqx(j,i)/sigx(i)
            end do
         end do
      else
         do i=1,dimx
            do j=1,ys1
               y579(i,j)=dfx(j,i)
            end do
            y579(i,ys1+i)=1D0
            do j=1,mrhs
               y579(i,ys2+j)=rhs(bidx+i,j)
            end do
            do j=1,nqx
               temp(i,j)=dqx(j,i)
            end do
         end do

C         call choslv(dimx,ys2+mrhs,dimx,0,sigx,y579,err)
         call chlsco(dimx,ys2+mrhs,dimx,0,sigx,y579,err,pinfo)
         if (err.ne.0) then
            return
         end if
C         call rrslv(dimx,nqx,dimx,0,sigx,temp)
         call rrslvo(dimx,nqx,dimx,0,sigx,temp,pinfo)
      end if

      do i=1,nqx
         do j=1,ys1
            val=0
            do k=1,dimx
               val=val+dqx(i,k)*y579(k,j)
            end do
            y6810(i,j)=val
         end do
         do j=1,dimx
            y6810(i,ys1+j)=temp(j,i)
         end do
         do j=1,mrhs
            val=0
            do k=1,dimx
               val=val+dqx(i,k)*y579(k,ys2+j)
            end do
            y6810(i,ys2+j)=val-rhs(bidx+dimx+i,j)
         end do
      end do

      do i=1,nqx
         do j=1,nqx
            val=0
            do k=1,dimx
               val=val+dqx(i,k)*temp(k,j)
            end do
            D(i,j)=val
         end do
      end do

      call choslv(nqx,ys2+mrhs,nqx,0,D,y6810,err)
      if (err.ne.0) then
         return
      end if

      do i=1,dimx
         do j=1,ys2+mrhs
            val=0
            do k=1,nqx
               val=val+temp(i,k)*y6810(k,j)
            end do
            y579(i,j)=y579(i,j)-val;
         end do
      end do

      end

      subroutine evbxrl(dimx,sigx,isq,nrhs,
     \ mrhs,rhs,dqx,nqx,bidx,Q,R,temp,D,DN,cNidx,err,pinfo)

      implicit none

      integer isq,nrhs,mrhs,bidx,dimx,nqx
      integer i,j,k,err,cNidx,pinfo
      double precision sigx,rhs,Q,R,dqx,temp,D,d1,DN

      dimension sigx(*)
      dimension Q(dimx+nqx,dimx),R(dimx+nqx,mrhs)
      dimension rhs(nrhs,mrhs)
      dimension temp(dimx,nqx),dqx(nqx,dimx)
      dimension D(nqx,nqx),DN(dimx,dimx),pinfo(*)

      err=0
      Q=0

      if (isq.eq.0) then
         do i=1,dimx
            do j=1,nqx
               temp(i,j)=dqx(j,i)/sigx(i)
            end do
            Q(i,i)=-1D0/sigx(i)
            do j=1,mrhs
               R(i,j)=rhs(bidx+i,j)/sigx(i)
            end do
         end do
      else
         do i=1,dimx
            do j=1,nqx
               temp(i,j)=dqx(j,i)
            end do
            Q(i,i)=-1D0
            do j=1,mrhs
               R(i,j)=rhs(bidx+i,j)
            end do
         end do
         call chlsco(dimx,dimx,dimx+nqx,0,sigx,Q,err,pinfo)
C         call choslv(dimx,dimx,dimx+nqx,0,sigx,Q,err)
         if (err.ne.0) then
            return
         end if
C         call rrslv(dimx,nqx,dimx,0,sigx,temp)
C         call rrslv(dimx,mrhs,dimx+nqx,0,sigx,R)
         call rrslvo(dimx,nqx,dimx,0,sigx,temp,pinfo)
         call rrslvo(dimx,mrhs,dimx+nqx,0,sigx,R,pinfo)
      end if

      do i=1,nqx
         do j=1,nqx
            d1=0D0
            do k=1,dimx
               d1=d1+dqx(i,k)*temp(k,j)
            end do
            D(i,j)=d1
         end do
         do j=1,dimx
            Q(dimx+i,j)=-temp(j,i)
         end do

         do j=1,mrhs
            d1=-rhs(bidx+dimx+i,j)
            do k=1,dimx
               d1=d1+dqx(i,k)*R(k,j)
            end do
            R(dimx+i,j)=d1
         end do
      end do

      call choslv(nqx,dimx,dimx+nqx,dimx,D,Q,err)
      if (err.ne.0) then
         return
      end if
      call rrslv(nqx,mrhs,dimx+nqx,dimx,D,R)

      do i=1,dimx
         do j=1,dimx
            d1=0D0
            do k=1,nqx
               d1=d1+temp(i,k)*Q(dimx+k,j)
            end do
            Q(i,j)=Q(i,j)-d1
            DN(i,j)=DN(i,j)-Q(i,j)
         end do
         do j=1,mrhs
            d1=0D0
            do k=1,nqx
               d1=d1+temp(i,k)*R(dimx+k,j)
            end do
            R(i,j)=R(i,j)-d1
            rhs(cNidx+i,j)=rhs(cNidx+i,j)-R(i,j)
         end do
      end do

      end


      subroutine evmrpr(dimu,dimx,stages,dfu,dfx,dfxA,xi,
     \ y1,y3,y5,y7,y9,bidx,bsz,rhs,nrhs,mrhs,
     \ nqu,nqx,rkb,hstep,
     \ Q1,Q2,R,err,D,iQ2)

      implicit none

      integer dimu,dimx,stages,bidx,bsz,nrhs,mrhs,nqu,nqx,err,iQ2
      integer i,j,k,dimk,tidx5,tidx4,tidx1,tidx2,tidx3
      integer ridx5

      double precision dfu,dfx,dfxA,xi
      double precision y1,y3,y5,y7,y9
      double precision rhs,rkb,hstep
      double precision Q1,Q2,R,D
      double precision d1

      dimension dfx(dimx*stages,dimx)
      dimension dfu(dimx*stages,dimu)
      dimension dfxA(dimx*stages,dimx*stages)
      dimension y1(dimu,dimx*stages)
      dimension y3(dimu,mrhs)
      dimension y5(dimx,dimx*stages)
      dimension y7(dimx,dimx)
      dimension y9(dimx,mrhs)
      dimension rhs(nrhs,mrhs),D(stages*dimx,stages*dimx)
      dimension Q1(bsz,dimx),Q2(bsz,dimx),R(bsz,3)
      dimension rkb(stages)

      dimk=dimx*stages
      tidx1=bidx+dimu
      tidx2=tidx1+nqu
      tidx3=tidx2+dimx
      tidx4=tidx3+nqx
      tidx5=tidx4+dimk

      ridx5=tidx5-bidx

      err=0

      do i=1,dimk
         do j=1,dimk
            d1=0D0
            do k=1,dimu
               d1=d1+dfu(i,k)*y1(k,j)
            end do
            do k=1,dimx
               d1=d1+dfx(i,k)*y5(k,j)
            end do
            do k=1,dimk
               d1=d1+dfxA(i,k)*dfxA(j,k)/xi
            end do
            D(i,j)=d1
         end do
      end do

      do i=1,dimk
         do j=1,dimx
            d1=0D0
            do k=1,dimx
               d1=d1+dfx(i,k)*y7(k,j)
            end do
            if (iQ2.eq.1) then
               Q2(ridx5+i,j)=-d1
            end if
            do k=1,stages
               d1=d1+dfxA(i,(k-1)*dimx+j)*rkb(k)*hstep/xi
            end do
            Q1(ridx5+i,j)=d1
         end do
         do j=1,mrhs
            d1=-rhs(tidx5+i,j)
            do k=1,dimu
               d1=d1+dfu(i,k)*y3(k,j)
            end do
            do k=1,dimx
               d1=d1+dfx(i,k)*y9(k,j)
            end do
            do k=1,dimk
               d1=d1+dfxA(i,k)*rhs(tidx4+k,j)/xi
            end do
            R(ridx5+i,j)=d1
         end do
      end do

      end

      subroutine evmrge(dimu,dimx,stages,dfxA,xi,
     \ y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,bidx,bsz,rhs,nrhs,mrhs,
     \ nqu,nqx,rkb,hstep,
     \ Q1,Q2,R,err,D,iQ2,pinfo)

      implicit none

      integer dimu,dimx,stages,bidx,bsz,nrhs,mrhs,nqu,nqx,err,iQ2
      integer i,j,k,dimk,tidx5,tidx4,tidx1,tidx2,tidx3,pinfo
      integer ridx1,ridx2,ridx3,ridx4,ridx5,ridx

      double precision dfxA,xi
      double precision y1,y2,y3,y4,y5,y6,y7,y8,y9,y10
      double precision rhs,rkb,hstep
      double precision Q1,Q2,R,D
      double precision d1

      dimension dfxA(dimx*stages,dimx*stages)
      dimension y1(dimu,dimx*stages)
      dimension y2(nqu,dimx*stages)
      dimension y3(dimu,mrhs)
      dimension y4(nqu,mrhs)
      dimension y5(dimx,dimx*stages)
      dimension y6(nqx,dimx*stages)
      dimension y7(dimx,dimx)
      dimension y8(nqx,dimx)
      dimension y9(dimx,mrhs)
      dimension y10(nqx,mrhs)
      dimension rhs(nrhs,mrhs),D(stages*dimx,stages*dimx)
      dimension Q1(bsz,dimx),Q2(bsz,dimx),R(bsz,3)
      dimension rkb(stages)
      dimension pinfo(*)

      dimk=dimx*stages
      tidx1=bidx+dimu
      tidx2=tidx1+nqu
      tidx3=tidx2+dimx
      tidx4=tidx3+nqx
      tidx5=tidx4+dimk

      ridx1=tidx1-bidx
      ridx2=tidx2-bidx
      ridx3=tidx3-bidx
      ridx4=tidx4-bidx
      ridx5=tidx5-bidx

      ridx=dimx+mrhs

      err=0

C      call choslv(dimk,dimx,bsz,ridx5,D,Q1,err)
      call chlsco(dimk,dimx,bsz,ridx5,D,Q1,err,pinfo)
      if (err.ne.0) then
         return
      end if

C      call rrslv(dimk,mrhs,bsz,ridx5,D,R)
      call rrslvo(dimk,mrhs,bsz,ridx5,D,R,pinfo)
      if (iQ2.eq.1) then
C         call rrslv(dimk,dimx,bsz,ridx5,D,Q2)
         call rrslvo(dimk,dimx,bsz,ridx5,D,Q2,pinfo)
      end if

      do j=1,dimx
         do i=1,dimu
            d1=0d0
            do k=1,dimk
               d1=d1-y1(i,k)*Q1(ridx5+k,j)
            end do
            Q1(i,j)=d1
         end do
         do i=1,nqu
            d1=0d0
            do k=1,dimk
               d1=d1-y2(i,k)*Q1(ridx5+k,j)
            end do
            Q1(ridx1+i,j)=d1
         end do
         do i=1,dimx
            d1=y7(i,j)
            do k=1,dimk
               d1=d1-y5(i,k)*Q1(ridx5+k,j)
            end do
            Q1(ridx2+i,j)=d1
         end do
         do i=1,nqx
            d1=y8(i,j)
            do k=1,dimk
               d1=d1-y6(i,k)*Q1(ridx5+k,j)
            end do
            Q1(ridx3+i,j)=d1
         end do
         do i=1,dimk
            d1=0d0
            do k=1,dimk
               d1=d1-dfxA(k,i)*Q1(ridx5+k,j)
            end do
            Q1(ridx4+i,j)=d1/xi
         end do
      end do

      do i=1,stages
         d1=hstep*rkb(i)/xi
         do j=1,dimx
            Q1(ridx4+(i-1)*dimx+j,j)=Q1(ridx4+(i-1)*dimx+j,j)+d1
         end do
      end do

      if (iQ2.eq.1) then
         do j=1,dimx
            do i=1,dimu
               d1=0d0
               do k=1,dimk
                  d1=d1-y1(i,k)*Q2(ridx5+k,j)
               end do
               Q2(i,j)=d1
            end do
            do i=1,nqu
               d1=0d0
               do k=1,dimk
                  d1=d1-y2(i,k)*Q2(ridx5+k,j)
               end do
               Q2(ridx1+i,j)=d1
            end do
            do i=1,dimx
               d1=-y7(i,j)
               do k=1,dimk
                  d1=d1-y5(i,k)*Q2(ridx5+k,j)
               end do
               Q2(ridx2+i,j)=d1
            end do
            do i=1,nqx
               d1=-y8(i,j)
               do k=1,dimk
                  d1=d1-y6(i,k)*Q2(ridx5+k,j)
               end do
               Q2(ridx3+i,j)=d1
            end do
            do i=1,dimk
               d1=0d0
               do k=1,dimk
                  d1=d1-dfxA(k,i)*Q2(ridx5+k,j)
               end do
               Q2(ridx4+i,j)=d1/xi
            end do
         end do
      end if

      do j=1,mrhs
         do i=1,dimu
            d1=y3(i,j)
            do k=1,dimk
               d1=d1-y1(i,k)*R(ridx5+k,j)
            end do
            R(i,j)=d1
         end do
         do i=1,nqu
            d1=y4(i,j)
            do k=1,dimk
               d1=d1-y2(i,k)*R(ridx5+k,j)
            end do
            R(ridx1+i,j)=d1
         end do
         do i=1,dimx
            d1=y9(i,j)
            do k=1,dimk
               d1=d1-y5(i,k)*R(ridx5+k,j)
            end do
            R(ridx2+i,j)=d1
         end do
         do i=1,nqx
            d1=y10(i,j)
            do k=1,dimk
               d1=d1-y6(i,k)*R(ridx5+k,j)
            end do
            R(ridx3+i,j)=d1
         end do
         do i=1,dimk
            d1=rhs(tidx4+i,j)
            do k=1,dimk
               d1=d1-dfxA(k,i)*R(ridx5+k,j)
            end do
            R(ridx4+i,j)=d1/xi
         end do
      end do

      end

      subroutine evdn(dimu,dimx,stages,
     \ bidx,bsz,rhs,nrhs,mrhs,
     \ nqu,nqx,rkb,hstep,
     \ Q1,Q2,R,iQ2,DN1,DN2,DN3,lbsz,cNidx)

      implicit none

      integer dimu,dimx,stages,bidx,bsz,nrhs,mrhs,nqu,nqx,iQ2
      integer i,j,k,dimk,tidx4,tidx1,tidx2,tidx3
      integer ridx2,ridx4,lbsz,cNidx

      double precision rhs,rkb,hstep
      double precision DN1,DN2,DN3
      double precision Q1,Q2,R
      double precision d1

      dimension rhs(nrhs,mrhs)
      dimension Q1(bsz,dimx),Q2(bsz,dimx),R(bsz,3)
      dimension rkb(stages)
      dimension DN1(dimx,dimx),DN2(dimx,dimx),DN3(dimx,dimx)

      dimk=dimx*stages
      tidx1=bidx+dimu
      tidx2=tidx1+nqu
      tidx3=tidx2+dimx
      tidx4=tidx3+nqx

      ridx2=tidx2-bidx
      ridx4=tidx4-bidx

      do i=1,dimx
         do j=1,dimx
            d1=Q1(ridx2+i,j)
            do k=1,stages
               d1=d1 + rkb(k)*hstep*Q1(ridx4+(k-1)*dimx+i,j)
            end do
            DN1(i,j)=d1
         end do
         do j=1,mrhs
            d1=-rhs(lbsz+cnIdx+i,j)+R(ridx2+i,j)
            do k=1,stages
               d1=d1 + rkb(k)*hstep*R(ridx4+(k-1)*dimx+i,j)
            end do
            rhs(lbsz+cnIdx+i,j)=d1
         end do
      end do

      if (iQ2.eq.1) then
         do i=1,dimx
            do j=1,dimx
               DN2(i,j)=DN2(i,j)-Q2(ridx2+i,j)
               DN3(i,j)=-Q1(ridx2+j,i)
            end do
            do j=1,mrhs
               rhs(lbsz+cnIdx-dimx+i,j)=rhs(lbsz+cnIdx-dimx+i,j)-
     \  R(ridx2+i,j)
            end do
         end do
      end if


      end

      subroutine evsol(N,dimx,bidx,ynidx,bsz,rhs,nrhs,mrhs,
     \ Q1,Q2,R,iQ2,P,sol)

      implicit none

      integer N,dimx,bidx,bsz,nrhs,mrhs,iQ2,ynidx,P
      integer i,j,k

      double precision rhs,Q1,Q2,R,sol
      double precision d1

      dimension rhs(nrhs,mrhs)
      dimension Q1(bsz,dimx),Q2(bsz,dimx),R(bsz,3)
      dimension P(*),sol(nrhs,mrhs)

      do i=1,bsz
         do j=1,mrhs
            d1=R(i,j)
            do k=1,dimx
               d1=d1-Q1(i,k)*rhs(ynidx+k,j)
            end do
            if (iQ2.eq.1) then
               do k=1,dimx
                  d1=d1-Q2(i,k)*rhs(ynidx-dimx+k,j)
               end do
            end if
            sol(P(bidx+i),j)=d1
         end do
      end do

      end