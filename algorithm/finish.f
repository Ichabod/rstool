      subroutine ifin(gsize,hsize,g,h,ieps,eps,cinf,infstp,
     \  stp,maxstp,mu,err,iloop)

      implicit none

      integer gsize,hsize,cinf,err,infstp,iloop
      integer stp,maxstp
      integer i

      double precision g,h,ieps,eps,mu
      double precision d

      dimension g(gsize),h(hsize)

      d=ieps

      err=0
      ieps=0D0
      iloop=0
      do i=1,gsize
         if (g(i).lt.0D0) then
            ieps=max(ieps,abs(g(i)))
         end if
      end do
      do i=1,hsize
         ieps=max(ieps,abs(h(i)))
      end do
      if (ieps/d.gt.0.95D0) then
         cinf=cinf+1
         if (cinf.ge.infstp) then
            err=1
            return
         end if
      else
         cinf=0
      end if

      stp=stp+1

      if (ieps.gt.mu.and.ieps.gt.eps.and.stp.lt.maxstp) then
         iloop=1
      end if


      end


      subroutine ofin(gsize,s,z,mu,mtau,oeps,eps,
     \  stp,maxstp,ostp,oloop)

      implicit none

      integer gsize,oloop
      integer stp,maxstp,ostp
      integer i

      double precision s,z,oeps,eps,mu,mtau
      double precision d,minpr,sz,xi,sig

      dimension s(gsize),z(gsize)

      minpr=huge(0D0)
      sz=0D0
      do i=1,gsize
         d=s(i)*z(i)
         minpr=min(minpr,abs(d))
         sz=sz+d
      end do
      sz=sz/gsize
      xi=minpr/sz
      sig=5D-2*(1D0-xi)/xi
      if (sig.gt.2D0) then
         sig=2D0
      end if
      sig=sig*1D-1
      mu=sig*sz
      mtau=sig*1E-1

      ostp=ostp+1

      if (oeps.gt.eps.and.stp.lt.maxstp) then
         oloop=1
      else 
         oloop=0;
      end if

      end