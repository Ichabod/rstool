      subroutine initip(xsize, gsize, hsize, istart,
     \      mu, mtau, ceps, steps, osteps, cinf, 
     \      s, z, y)

      implicit none

      integer xsize,gsize,hsize,istart,steps,osteps,cinf
      double precision mu,mtau,ceps,s,y,z

      dimension s(gsize),z(gsize),y(hsize)

      mu=0.99D0
      mtau=0.005D0
      ceps=1D19
      steps=0
      osteps=0
      cinf=0

      if (istart.eq.1) then
         s=0.2D0
         z=0.01D0
         y=0.01D0
      end if

      end

