*     Temporärer Speicher: dimx + dimx*dimx + dimx*dimu

      subroutine evrest(x,p,t0,N,hstep,dimx,dimu,xsize,gsize,hsize,
     \   s,A,b,c,
     \   rhs,rhsx,rhsu,r_u,dr_u,r_x,dr_x,q_u,dq_u,q_x,dq_x,
     \   nub_u,nlb_u,nub_x,nlb_x,
     \   iub_u,ilb_u,iub_x,ilb_x,
     \   pub_u,plb_u,pub_x,plb_x,
     \   ub_u,lb_u,ub_x,lb_x,
     \   nr_u,nr_x,nq_u,nq_x,
     \   idr_u,idr_x,idq_u,idq_x,
     \   igub_u,iglb_u,igr_u,igub_x,iglb_x,igr_x,
     \   ihq_u,ihq_x,
     \   ihk,ihx,
     \   g,h,dfu,dfx,dfxA,dqu,dqx,dru,drx,work,lwork,dograd,err)

      implicit none
      integer N,dimx,dimu,s,dograd,xsize,gsize,hsize
      integer nub_u,nlb_u,nub_x,nlb_x
      integer iub_u,ilb_u,iub_x,ilb_x
      integer pub_u,plb_u,pub_x,plb_x
      integer nr_u,nr_x,nq_u,nq_x
      integer idr_u,idr_x,idq_u,idq_x
      integer igub_u,iglb_u,igr_u,igub_x,iglb_x,igr_x
      integer ihq_u,ihq_x
      integer ihk,ihx
      integer lwork,err
      double precision x,p,t0,hstep,A,b,c,g,h,work
      double precision ub_u,lb_u,ub_x,lb_x
      double precision dfu,dfx,dfxA,dqu,dqx,dru,drx

      dimension nub_u(*),nlb_u(*),nub_x(*),nlb_x(*)
      dimension iub_u(*),ilb_u(*),iub_x(*),ilb_x(*)
      dimension pub_u(*),plb_u(*),pub_x(*),plb_x(*)
      dimension ub_u(*),lb_u(*),ub_x(*),lb_x(*)
      dimension igub_u(*),iglb_u(*),igr_u(*)
      dimension igub_x(*),iglb_x(*),igr_x(*)
      dimension ihq_u(*),ihq_x(*)
      dimension ihk(*),ihx(*)
      dimension nr_u(*),nq_u(*),nr_x(*),nq_x(*)
      dimension idr_u(*),idq_u(*),idr_x(*),idq_x(*)
      dimension x(*),p(*),A(*),b(*),c(*),g(*),h(*),work(*)
      dimension dfu(*),dfx(*),dfxA(*),dqu(*),dqx(*),dru(*),drx(*)

      external rhs,rhsx,rhsu,r_u,dr_u,r_x,dr_x,q_u,dq_u,q_x,dq_x

      integer iu,ix,ik
      integer itx,itfx,itfu

      iu=1
      ix=iu+N*dimu
      ik=ix+(N+1)*dimx
      itx=1
      itfx=itx+dimx
      itfu=itfx+dimx*dimx
      if (itfu+dimx*dimu .gt. lwork+1) then
         err = 1
         return
      end if


      call evreco(x(iu),x(ix),x(ik),p,t0,N,hstep,dimx,dimu,dimx*s,
     \ xsize,gsize,hsize,s,A,b,c,
     \ rhs,rhsx,rhsu,r_u,dr_u,r_x,dr_x,q_u,dq_u,q_x,dq_x,
     \ nub_u,nlb_u,nub_x,nlb_x,
     \ iub_u,ilb_u,iub_x,ilb_x,
     \ pub_u,plb_u,pub_x,plb_x,
     \ ub_u,lb_u,ub_x,lb_x,
     \ nr_u,nr_x,nq_u,nq_x,
     \ idr_u,idr_x,idq_u,idq_x,
     \ igub_u,iglb_u,igr_u,igub_x,iglb_x,igr_x,
     \ ihq_u,ihq_x,
     \ ihk,ihx,
     \ g,h,dfu,dfx,dfxA,dqu,dqx,dru,drx,
     \ work(itx),work(itfx),work(itfu),dograd)


      end


      subroutine evreco(u,x,k,p,t0,N,hstep,dimx,dimu,dimk,
     \ xsize,gsize,hsize,s,A,b,c,
     \ rhs,rhsx,rhsu,r_u,dr_u,r_x,dr_x,q_u,dq_u,q_x,dq_x,
     \ nub_u,nlb_u,nub_x,nlb_x,
     \ iub_u,ilb_u,iub_x,ilb_x,
     \ pub_u,plb_u,pub_x,plb_x,
     \ ub_u,lb_u,ub_x,lb_x,
     \ nr_u,nr_x,nq_u,nq_x,
     \ idr_u,idr_x,idq_u,idq_x,
     \ igub_u,iglb_u,igr_u,igub_x,iglb_x,igr_x,
     \ ihq_u,ihq_x,
     \ ihk,ihx,
     \ g,h,dfu,dfx,dfxA,dqu,dqx,dru,drx,
     \ tx,tfx,tfu,dograd)

      implicit none
      integer N,dimx,dimu,dimk,s,dograd,xsize,gsize,hsize
      integer nub_u,nlb_u,nub_x,nlb_x
      integer iub_u,ilb_u,iub_x,ilb_x
      integer pub_u,plb_u,pub_x,plb_x
      integer nr_u,nr_x,nq_u,nq_x
      integer idr_u,idr_x,idq_u,idq_x
      integer igub_u,iglb_u,igr_u,igub_x,iglb_x,igr_x
      integer ihq_u,ihq_x
      integer ihk,ihx
      integer i,j,l,m,o
      integer ofs,ofs2
      double precision u,x,k,p,t0,hstep,A,b,c,g,h
      double precision ub_u,lb_u,ub_x,lb_x
      double precision t,tx,tfx,tfu,tt
      double precision dfu,dfx,dfxA,dqu,dqx,dru,drx

      dimension nub_u(*),nlb_u(*),nub_x(*),nlb_x(*)
      dimension iub_u(*),ilb_u(*),iub_x(*),ilb_x(*)
      dimension pub_u(*),plb_u(*),pub_x(*),plb_x(*)
      dimension ub_u(*),lb_u(*),ub_x(*),lb_x(*)
      dimension igub_u(*),iglb_u(*),igr_u(*)
      dimension igub_x(*),iglb_x(*),igr_x(*)
      dimension ihq_u(*),ihq_x(*)
      dimension ihk(s,*),ihx(*)
      dimension nr_u(*),nq_u(*),nr_x(*),nq_x(*)
      dimension idr_u(*),idq_u(*),idr_x(*),idq_x(*)
      dimension u(dimu,N),x(dimx,N+1),k(dimx,s,N),p(*)
      dimension A(s,s),b(s),c(s)
      dimension g(*),h(*)
      dimension tx(dimx),tfx(dimx,dimx),tfu(dimx,dimu)
      dimension dfu(dimk,dimu,N),dfx(dimk,dimx,N)
      dimension dfxA(dimk,dimk,N),dqu(*),dqx(*),dru(*),drx(*)

      external rhs,rhsx,rhsu,r_u,dr_u,r_x,dr_x,q_u,dq_u,q_x,dq_x

      t=t0

      do i=1,N

*        RK-Stufen
         ofs=0
         do j=1,s
            call evrc1(j,ofs,u(1,i),x(1,i),k(1,1,i),p,t,
     \ hstep,dimx,dimu,dimk,s,A,b,c,rhs,rhsx,rhsu,
     \ h(ihk(j,i)+1),dfu(1,1,i),dfx(1,1,i),dfxA(1,1,i),
     \ tx,tfx,tfu,dograd)

            ofs=ofs+dimx
         end do

         call evrc3(dimx,s,x(1,i),x(1,i+1),b,hstep,
     \ k(1,1,i),tx, h(ihx(i)+1))

         call evrc2(i,u(1,i),p,N,r_u,dr_u,q_u,dq_u,
     \ nub_u(i),nlb_u(i),pub_u(iub_u(i)),plb_u(ilb_u(i)),
     \ ub_u(iub_u(i)),lb_u(ilb_u(i)),nr_u(i),nq_u(i),
     \ g(igub_u(i)+1),g(iglb_u(i)+1),g(igr_u(i)+1),h(ihq_u(i)+1),
     \ dqu(idq_u(i)+1),dru(idr_u(i)+1),dograd)

         t=t+hstep
      end do

      do i=1,N+1

         call evrc2(i,x(1,i),p,N,r_x,dr_x,q_x,dq_x,
     \ nub_x(i),nlb_x(i),pub_x(iub_x(i)),plb_x(ilb_x(i)),
     \ ub_x(iub_x(i)),lb_x(ilb_x(i)),nr_x(i),nq_x(i),
     \ g(igub_x(i)+1),g(iglb_x(i)+1),g(igr_x(i)+1),h(ihq_x(i)+1),
     \ dqx(idq_x(i)+1),drx(idr_x(i)+1),dograd)      
         
      end do


      end


      subroutine evrc1(j,ofs,u,x,k,p,t,hstep,dimx,dimu,dimk,
     \ s,A,b,c,
     \ rhs,rhsx,rhsu,
     \ h,dfu,dfx,dfxA,
     \ tx,tfx,tfu,dograd)

      implicit none
      integer dimx,dimu,dimk,s,dograd
      integer j,l,m,o
      integer ofs,ofs2
      double precision u,x,k,p,t,hstep,A,b,c,h
      double precision tx,tfx,tfu,tt
      double precision dfu,dfx,dfxA

      dimension u(dimu),x(dimx),k(dimx,s),p(*)
      dimension A(s,s),b(s),c(s)
      dimension h(*)
      dimension tx(dimx),tfx(dimx,dimx),tfu(dimx,dimu)
      dimension dfu(dimk,dimu),dfx(dimk,dimx)
      dimension dfxA(dimk,dimk)

      external rhs,rhsx,rhsu


      tx=0D0
      do l=1,s
         do m=1,dimx
            tx(m)=tx(m)+A(j,l)*k(m,l)
         end do
      end do
      do m=1,dimx
         tx(m)=hstep*tx(m)+x(m)
      end do
      tt=t+hstep*c(j)

      if (dograd.eq.1) then

*        Funktionen auswerten
         call rhsx(tt,tx,u,p,tfx)
         call rhsu(tt,tx,u,p,tfu)

*        Werte in dfx und dfu kopieren
         do l=1,dimx
            do m=1,dimx
               dfx(ofs+l,m)=tfx(l,m)
            end do
            do m=1,dimu
               dfu(ofs+l,m)=tfu(l,m)
            end do
         end do

*        dfxA berechnen...
         ofs2=0
         do l=1,s
            do m=1,dimx
               do o=1,dimx
                  dfxA(ofs+m,ofs2+o)=tfx(m,o)*hstep*A(j,l)
               end do
            end do
            if (l.eq.j) then
               do  m=1,dimx
                  dfxA(ofs+m,ofs2+m)=dfxA(ofs+m,ofs2+m)-1D0
               end do
            end if
            ofs2=ofs2+dimx
         end do
      end if


*     Gleichungsrestriktionen speichern
      call rhs(tt,tx,u,p,h)
      do m=1,dimx
         h(m)=h(m)-k(m,j)
      end do


      end


      subroutine evrc2(i,x,p,N,
     \ fr,fdr,fq,fdq,
     \ nub,nlb,
     \ pub,plb,
     \ ub,lb,
     \ nr,nq,
     \ gub,glb,gr,h,dq,dr,dograd)

      implicit none
      integer N,dograd
      integer nub,nlb
      integer pub,plb
      integer nr,nq
      integer i,j
      double precision x,p
      double precision gub,glb,gr,h
      double precision ub,lb
      double precision dq,dr

      dimension pub(*),plb(*)
      dimension ub(*),lb(*)
      dimension x(*),p(*)
      dimension glb(*),gub(*),gr(*),h(*)
      dimension dq(*),dr(*)

      external fr,fdr,fq,fdq


*     Boxed Constraints berechnen
      if (nub.gt.0) then
         do j=1,nub
            gub(j)=ub(j)-x(pub(j))
         end do
      end if
      if (nlb.gt.0) then
         do j=1,nlb
            glb(j)=x(plb(j))-lb(j)
         end do
      end if

*     Nichtlineare Ungleichungen
      if (nr.gt.0) then
         call fr(i,N,x,p,gr)
         if (dograd.eq.1) then
            call fdr(i,N,x,p,dr)
         end if
      end if

*     Nichtlineare Gleichungen
      if (nq.gt.0) then
         call fq(i,N,x,p,h)
         if (dograd.eq.1) then
            call fdq(i,N,x,p,dq)
         end if
      end if

      end


      subroutine evrc3(dimx,s,x1,x2,b,hstep,k,tx,h)

      implicit none

      integer dimx,s,l,m
      double precision x1,x2,b,hstep,k,tx,h
      
      dimension x1(dimx),x2(dimx)
      dimension b(s),k(dimx,s),tx(dimx),h(dimx)

      tx=0D0
      do l=1,s
         do m=1,dimx
            tx(m)=tx(m)+b(l)*k(m,l)
         end do
      end do
      do m=1,dimx
         h(m)=x1(m)+hstep*tx(m) - x2(m)
      end do

      end