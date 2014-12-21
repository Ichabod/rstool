      subroutine rhs (t, x, u, p, cgret)
        double precision t
        double precision x(*)
        double precision u(*)
        double precision p(*)
        double precision cgret(13)
        double precision t1
        double precision t11
        double precision t12
        double precision t13
        double precision t14
        double precision t15
        double precision t16
        double precision t17
        double precision t18
        double precision t2
        double precision t20
        double precision t22
        double precision t23
        double precision t26
        double precision t28
        double precision t29
        double precision t32
        double precision t36
        double precision t4
        double precision t46
        double precision t47
        double precision t67
        double precision t68
        double precision t7
        double precision t70
        double precision t71
        double precision t76
        double precision t79
        t1 = x(4)
        t2 = x(5)
        t4 = p(1)
        t7 = t4 ** 2
        t11 = x(10)
        t12 = t11 ** 2
        t13 = x(11)
        t14 = t13 ** 2
        t15 = x(12)
        t16 = t15 ** 2
        t17 = x(13)
        t18 = t17 ** 2
        t20 = u(1)
        t22 = t11 * t13
        t23 = t15 * t17
        t26 = u(2)
        t28 = t11 * t15
        t29 = t13 * t17
        t32 = u(3)
        t36 = 0.1D1 / p(3)
        t46 = t11 * t17
        t47 = t13 * t15
        t67 = x(8)
        t68 = x(9)
        t70 = p(5)
        t71 = p(6)
        t76 = p(4)
        t79 = x(7)
        cgret(1) = t1
        cgret(2) = t2
        cgret(3) = x(6)
        cgret(4) = 0.2D1 * t4 * t2 + 0.3D1 * t7 * x(1) + ((t12 - t14 - t
     #16 + t18) * t20 + (0.2D1 * t22 - 0.2D1 * t23) * t26 + (0.2D1 * t28
     # + 0.2D1 * t29) * t32) * t36
        cgret(5) = -0.2D1 * t4 * t1 + ((0.2D1 * t22 + 0.2D1 * t23) * t20
     # + (-t12 + t14 - t16 + t18) * t26 + (-0.2D1 * t46 + 0.2D1 * t47) *
     # t32) * t36
        cgret(6) = -t7 * x(3) + ((0.2D1 * t28 - 0.2D1 * t29) * t20 + (0.
     #2D1 * t46 + 0.2D1 * t47) * t26 + (-t12 - t14 + t16 + t18) * t32) *
     # t36
        cgret(7) = (t67 * t68 * (t70 - t71) + u(4)) / t76
        cgret(8) = (t79 * t68 * (t71 - t76) + u(5)) / t70
        cgret(9) = (t79 * t67 * (t76 - t70) + u(6)) / t71
        cgret(10) = t13 * t68 / 0.2D1 - t15 * t67 / 0.2D1 + t17 * t79 / 
     #0.2D1
        cgret(11) = -t11 * t68 / 0.2D1 + t15 * t79 / 0.2D1 + t17 * t67 /
     # 0.2D1
        cgret(12) = t11 * t67 / 0.2D1 - t13 * t79 / 0.2D1 + t17 * t68 / 
     #0.2D1
        cgret(13) = -t11 * t79 / 0.2D1 - t13 * t67 / 0.2D1 - t15 * t68 /
     # 0.2D1
      end

      subroutine rhsx (t, x, u, p, cgret)
        double precision t
        double precision x(*)
        double precision u(*)
        double precision p(*)
        double precision cgret(13,13)
        double precision t11
        double precision t12
        double precision t14
        double precision t15
        double precision t20
        double precision t21
        double precision t24
        double precision t27
        double precision t28
        double precision t33
        double precision t39
        double precision t4
        double precision t40
        double precision t45
        double precision t5
        double precision t50
        double precision t51
        double precision t52
        double precision t53
        double precision t55
        double precision t56
        double precision t58
        double precision t62
        double precision t64
        double precision t66
        double precision t7
        double precision t70
        double precision t72
        double precision t77
        double precision t78
        double precision t79
        double precision t8
        double precision t80
        double precision t81
        double precision t82
        double precision t84
        double precision t9
        t4 = p(1)
        t5 = t4 ** 2
        t7 = 0.2D1 * t4
        t8 = x(10)
        t9 = u(1)
        t11 = x(11)
        t12 = u(2)
        t14 = x(12)
        t15 = u(3)
        t20 = 0.1D1 / p(3)
        t21 = (0.2D1 * t11 * t12 + 0.2D1 * t14 * t15 + 0.2D1 * t8 * t9) 
     #* t20
        t24 = x(13)
        t27 = -0.2D1 * t11 * t9 + 0.2D1 * t12 * t8 + 0.2D1 * t15 * t24
        t28 = t27 * t20
        t33 = -0.2D1 * t12 * t24 - 0.2D1 * t14 * t9 + 0.2D1 * t15 * t8
        t39 = 0.2D1 * t11 * t15 - 0.2D1 * t12 * t14 + 0.2D1 * t24 * t9
        t40 = t39 * t20
        t45 = -t33 * t20
        t50 = x(9)
        t51 = p(5)
        t52 = p(6)
        t53 = t51 - t52
        t55 = p(4)
        t56 = 0.1D1 / t55
        t58 = x(8)
        t62 = t52 - t55
        t64 = 0.1D1 / t51
        t66 = x(7)
        t70 = t55 - t51
        t72 = 0.1D1 / t52
        t77 = t24 / 0.2D1
        t78 = t14 / 0.2D1
        t79 = t11 / 0.2D1
        t80 = t50 / 0.2D1
        t81 = t58 / 0.2D1
        t82 = t66 / 0.2D1
        t84 = t8 / 0.2D1
        cgret(1,1) = 0D0
        cgret(1,2) = 0D0
        cgret(1,3) = 0D0
        cgret(1,4) = 1D0
        cgret(1,5) = 0D0
        cgret(1,6) = 0D0
        cgret(1,7) = 0D0
        cgret(1,8) = 0D0
        cgret(1,9) = 0D0
        cgret(1,10) = 0D0
        cgret(1,11) = 0D0
        cgret(1,12) = 0D0
        cgret(1,13) = 0D0
        cgret(2,1) = 0D0
        cgret(2,2) = 0D0
        cgret(2,3) = 0D0
        cgret(2,4) = 0D0
        cgret(2,5) = 1D0
        cgret(2,6) = 0D0
        cgret(2,7) = 0D0
        cgret(2,8) = 0D0
        cgret(2,9) = 0D0
        cgret(2,10) = 0D0
        cgret(2,11) = 0D0
        cgret(2,12) = 0D0
        cgret(2,13) = 0D0
        cgret(3,1) = 0D0
        cgret(3,2) = 0D0
        cgret(3,3) = 0D0
        cgret(3,4) = 0D0
        cgret(3,5) = 0D0
        cgret(3,6) = 1D0
        cgret(3,7) = 0D0
        cgret(3,8) = 0D0
        cgret(3,9) = 0D0
        cgret(3,10) = 0D0
        cgret(3,11) = 0D0
        cgret(3,12) = 0D0
        cgret(3,13) = 0D0
        cgret(4,1) = 0.3D1 * t5
        cgret(4,2) = 0D0
        cgret(4,3) = 0D0
        cgret(4,4) = 0D0
        cgret(4,5) = t7
        cgret(4,6) = 0D0
        cgret(4,7) = 0D0
        cgret(4,8) = 0D0
        cgret(4,9) = 0D0
        cgret(4,10) = t21
        cgret(4,11) = t28
        cgret(4,12) = t33 * t20
        cgret(4,13) = t40
        cgret(5,1) = 0D0
        cgret(5,2) = 0D0
        cgret(5,3) = 0D0
        cgret(5,4) = -t7
        cgret(5,5) = 0D0
        cgret(5,6) = 0D0
        cgret(5,7) = 0D0
        cgret(5,8) = 0D0
        cgret(5,9) = 0D0
        cgret(5,10) = -t27 * t20
        cgret(5,11) = t21
        cgret(5,12) = t40
        cgret(5,13) = t45
        cgret(6,1) = 0D0
        cgret(6,2) = 0D0
        cgret(6,3) = -t5
        cgret(6,4) = 0D0
        cgret(6,5) = 0D0
        cgret(6,6) = 0D0
        cgret(6,7) = 0D0
        cgret(6,8) = 0D0
        cgret(6,9) = 0D0
        cgret(6,10) = t45
        cgret(6,11) = -t39 * t20
        cgret(6,12) = t21
        cgret(6,13) = t28
        cgret(7,1) = 0D0
        cgret(7,2) = 0D0
        cgret(7,3) = 0D0
        cgret(7,4) = 0D0
        cgret(7,5) = 0D0
        cgret(7,6) = 0D0
        cgret(7,7) = 0D0
        cgret(7,8) = t50 * t53 * t56
        cgret(7,9) = t58 * t53 * t56
        cgret(7,10) = 0D0
        cgret(7,11) = 0D0
        cgret(7,12) = 0D0
        cgret(7,13) = 0D0
        cgret(8,1) = 0D0
        cgret(8,2) = 0D0
        cgret(8,3) = 0D0
        cgret(8,4) = 0D0
        cgret(8,5) = 0D0
        cgret(8,6) = 0D0
        cgret(8,7) = t50 * t62 * t64
        cgret(8,8) = 0D0
        cgret(8,9) = t66 * t62 * t64
        cgret(8,10) = 0D0
        cgret(8,11) = 0D0
        cgret(8,12) = 0D0
        cgret(8,13) = 0D0
        cgret(9,1) = 0D0
        cgret(9,2) = 0D0
        cgret(9,3) = 0D0
        cgret(9,4) = 0D0
        cgret(9,5) = 0D0
        cgret(9,6) = 0D0
        cgret(9,7) = t58 * t70 * t72
        cgret(9,8) = t66 * t70 * t72
        cgret(9,9) = 0D0
        cgret(9,10) = 0D0
        cgret(9,11) = 0D0
        cgret(9,12) = 0D0
        cgret(9,13) = 0D0
        cgret(10,1) = 0D0
        cgret(10,2) = 0D0
        cgret(10,3) = 0D0
        cgret(10,4) = 0D0
        cgret(10,5) = 0D0
        cgret(10,6) = 0D0
        cgret(10,7) = t77
        cgret(10,8) = -t78
        cgret(10,9) = t79
        cgret(10,10) = 0D0
        cgret(10,11) = t80
        cgret(10,12) = -t81
        cgret(10,13) = t82
        cgret(11,1) = 0D0
        cgret(11,2) = 0D0
        cgret(11,3) = 0D0
        cgret(11,4) = 0D0
        cgret(11,5) = 0D0
        cgret(11,6) = 0D0
        cgret(11,7) = t78
        cgret(11,8) = t77
        cgret(11,9) = -t84
        cgret(11,10) = -t80
        cgret(11,11) = 0D0
        cgret(11,12) = t82
        cgret(11,13) = t81
        cgret(12,1) = 0D0
        cgret(12,2) = 0D0
        cgret(12,3) = 0D0
        cgret(12,4) = 0D0
        cgret(12,5) = 0D0
        cgret(12,6) = 0D0
        cgret(12,7) = -t79
        cgret(12,8) = t84
        cgret(12,9) = t77
        cgret(12,10) = t81
        cgret(12,11) = -t82
        cgret(12,12) = 0D0
        cgret(12,13) = t80
        cgret(13,1) = 0D0
        cgret(13,2) = 0D0
        cgret(13,3) = 0D0
        cgret(13,4) = 0D0
        cgret(13,5) = 0D0
        cgret(13,6) = 0D0
        cgret(13,7) = -t84
        cgret(13,8) = -t79
        cgret(13,9) = -t78
        cgret(13,10) = -t82
        cgret(13,11) = -t81
        cgret(13,12) = -t80
        cgret(13,13) = 0D0
      end


      subroutine rhsu (t, x, u, p, cgret)
        double precision t
        double precision x(*)
        double precision u
        double precision p(*)
        double precision cgret(13,6)

        double precision t10,t11,t12,t13
        double precision tp3
        double precision t10s,t11s,t12s,t13s
        double precision t1011,t1213,t1012,t1113,t1013,t1112

        t10=x(10)
        t11=x(11)
        t12=x(12)
        t13=x(13)

        tp3=p(3)
        t10s = t10**2
        t11s = t11**2
        t12s = t12**2
        t13s = t13**2
        t1011 = 0.2D1 * t10*t11
        t1213 = 0.2D1 *t12*t13
        t1012 = 0.2D1 * t10*t12
        t1113 = 0.2D1 * t11*t13
        t1013 = 0.2D1 * t10*t13
        t1112 = 0.2D1 * t11*t12



        cgret(1,1) = 0D0
        cgret(1,2) = 0D0
        cgret(1,3) = 0D0
        cgret(1,4) = 0D0
        cgret(1,5) = 0D0
        cgret(1,6) = 0D0
        cgret(2,1) = 0D0
        cgret(2,2) = 0D0
        cgret(2,3) = 0D0
        cgret(2,4) = 0D0
        cgret(2,5) = 0D0
        cgret(2,6) = 0D0
        cgret(3,1) = 0D0
        cgret(3,2) = 0D0
        cgret(3,3) = 0D0
        cgret(3,4) = 0D0
        cgret(3,5) = 0D0
        cgret(3,6) = 0D0
        cgret(4,1) = (t10s - t11s -t12s + t13s)/ tp3
        cgret(4,2) = (t1011 - t1213) / tp3
        cgret(4,3) = (t1012 + t1113) / tp3
        cgret(4,4) = 0D0
        cgret(4,5) = 0D0
        cgret(4,6) = 0D0
        cgret(5,1) = (t1011 + t1213) / tp3
        cgret(5,2) = (-t10s + t11s -t12s + t13s) / tp3
        cgret(5,3) = (-t1013 + t1112) / tp3
        cgret(5,4) = 0D0
        cgret(5,5) = 0D0
        cgret(5,6) = 0D0
        cgret(6,1) = (t1012 - t1113) / tp3
        cgret(6,2) = (t1013 + t1112) / tp3
        cgret(6,3) = (-t10s - t11s +t12s + t13s) / tp3
        cgret(6,4) = 0D0
        cgret(6,5) = 0D0
        cgret(6,6) = 0D0
        cgret(7,1) = 0D0
        cgret(7,2) = 0D0
        cgret(7,3) = 0D0
        cgret(7,4) = 0.1D1 / p(4)
        cgret(7,5) = 0D0
        cgret(7,6) = 0D0
        cgret(8,1) = 0D0
        cgret(8,2) = 0D0
        cgret(8,3) = 0D0
        cgret(8,4) = 0D0
        cgret(8,5) = 0.1D1 / p(5)
        cgret(8,6) = 0D0
        cgret(9,1) = 0D0
        cgret(9,2) = 0D0
        cgret(9,3) = 0D0
        cgret(9,4) = 0D0
        cgret(9,5) = 0D0
        cgret(9,6) = 0.1D1 / p(6)
        cgret(10,1) = 0D0
        cgret(10,2) = 0D0
        cgret(10,3) = 0D0
        cgret(10,4) = 0D0
        cgret(10,5) = 0D0
        cgret(10,6) = 0D0
        cgret(11,1) = 0D0
        cgret(11,2) = 0D0
        cgret(11,3) = 0D0
        cgret(11,4) = 0D0
        cgret(11,5) = 0D0
        cgret(11,6) = 0D0
        cgret(12,1) = 0D0
        cgret(12,2) = 0D0
        cgret(12,3) = 0D0
        cgret(12,4) = 0D0
        cgret(12,5) = 0D0
        cgret(12,6) = 0D0
        cgret(13,1) = 0D0
        cgret(13,2) = 0D0
        cgret(13,3) = 0D0
        cgret(13,4) = 0D0
        cgret(13,5) = 0D0
        cgret(13,6) = 0D0
      end

      subroutine q_u(i,N,u,p,res)
         integer i,N
         double precision u,p,res

         dimension u(6)
         dimension p(*)
         dimension res(*)

      end

      subroutine dq_u(i,N,u,p,res)
         integer i,N
         double precision u,p,res

         dimension u(6)
         dimension p(*)
         dimension res(*)
      end

      subroutine nq_u(N,p,res)
         integer N,res
         double precision p
         dimension res(N)
         res=0
      end

      subroutine r_u(i,N,u,p,res)
         integer i,N
         double precision u,p,res

         dimension u(6)
         dimension p(*)
         dimension res(*)
      end

      subroutine dr_u(i,N,u,p,res)
         integer i,N
         double precision u,p,res

         dimension u(6)
         dimension p(*)
         dimension res(*)

      end

      subroutine nr_u(N,p,res)
         integer N,res
         double precision p
         dimension res(N)
         res=0
      end

      subroutine q_x(i,N,x,p,res)
         integer i,N,tidx,rtidx,drtidx,sz
         double precision x,p,rot,rsd,rsw,res

         dimension x(13)
         dimension p(*)
         dimension rot(3,3)
         dimension rsd(3),rsw(3)
         dimension res(*)

         if (i.eq.N+1) then
            call rotate(x,rot)
            call musmm(3,1,3,3,0,0,3,0,0,3,0,0,rsd,rot,p(7),1D0)
            call musmm(3,1,3,3,0,0,3,0,0,3,0,0,rsw,rot,x(7),1D0)
            tidx=13
            sz=3*(N+1)
            rtidx=tidx+sz
            drtidx=rtidx+sz
            call q_xsub(i,N,x,res,rsd,rsw,p(tidx),p(rtidx),p(drtidx))
         else if (i.eq.1) then
            res(1) = x(4)
            res(2) = x(5)
            res(3) = x(6)
            res(4) = x(7)
            res(5) = x(8)
            res(6) = x(9)
            res(7) = 30D0*x(10)
            res(8) = 30D0*x(11)
            res(9) = 30D0*x(12)
            res(10) = 30D0*(x(13)-1D0)
         end if
      end

      subroutine q_xsub(i,N,x,res,rsd,rsw,target,Rtd,dRtd)
         integer i,N
         double precision x,rsd,rsw,target,Rtd,dRtd,res

         dimension x(13)
         dimension rsd(3),rsw(3)
         dimension target(3,N+1)
         dimension Rtd(3,N+1)
         dimension dRtd(3,N+1)
         dimension res(*)

         res(1) = x(1) + rsd(1)-target(1,i)-Rtd(1,i)
         res(2) = x(2) + rsd(2)-target(2,i)-Rtd(2,i)
         res(3) = x(3) + rsd(3)-target(3,i)-Rtd(3,i)
         res(4) = x(4) + (rsw(2)*rsd(3) - rsw(3)*rsd(2))-dRtd(1,i)
         res(5) = x(5) + (rsw(3)*rsd(1) - rsw(1)*rsd(3))-dRtd(2,i)
         res(6) = x(6) + (rsw(1)*rsd(2) - rsw(2)*rsd(1))-dRtd(3,i)
      end

      subroutine dq_x(i,N,x,p,cgret)
         integer i,N
         double precision x(*)
         double precision p(*)
         double precision cgret(*)

         if (i.eq.N+1) then
            call dq_x1(i,N,x,p,cgret)  
         else if (i.eq.1) then
            call dq_x2(i,N,x,p,cgret)  
         end if
      end

      subroutine dq_x1(i,N,x,p,cgret)
         integer i,N
         double precision x(*)
         double precision p(*)
         double precision cgret(6,13)
         double precision t1
         double precision t101
         double precision t102
         double precision t108
         double precision t11
         double precision t110
         double precision t112
         double precision t113
         double precision t114
         double precision t115
         double precision t116
         double precision t119
         double precision t120
         double precision t121
         double precision t125
         double precision t128
         double precision t131
         double precision t133
         double precision t14
         double precision t143
         double precision t147
         double precision t148
         double precision t155
         double precision t157
         double precision t159
         double precision t17
         double precision t2
         double precision t22
         double precision t27
         double precision t29
         double precision t30
         double precision t32
         double precision t34
         double precision t35
         double precision t37
         double precision t38
         double precision t39
         double precision t4
         double precision t41
         double precision t43
         double precision t44
         double precision t46
         double precision t48
         double precision t49
         double precision t5
         double precision t50
         double precision t51
         double precision t52
         double precision t54
         double precision t57
         double precision t60
         double precision t62
         double precision t7
         double precision t71
         double precision t73
         double precision t75
         double precision t78
         double precision t8
         double precision t83
         double precision t89
         double precision t94

         t1 = x(10)
         t2 = p(7)
         t4 = x(11)
         t5 = p(8)
         t7 = x(12)
         t8 = p(9)
         t11 = 0.2D1 * t1 * t2 + 0.2D1 * t4 * t5 + 0.2D1 * t7 * t8
         t14 = x(13)
         t17 = -0.2D1 * t4 * t2 + 0.2D1 * t1 * t5 + 0.2D1 * t14 * t8
         t22 = -0.2D1 * t7 * t2 - 0.2D1 * t14 * t5 + 0.2D1 * t1 * t8
         t27 = 0.2D1 * t14 * t2 - 0.2D1 * t7 * t5 + 0.2D1 * t4 * t8
         t29 = -t17
         t30 = -t22
         t32 = -t27
         t34 = t1 * t4
         t35 = t7 * t14
         t37 = 0.2D1 * t34 + 0.2D1 * t35
         t38 = t1 * t7
         t39 = t4 * t14
         t41 = 0.2D1 * t38 - 0.2D1 * t39
         t43 = t4 * t7
         t44 = t1 * t14
         t46 = 0.2D1 * t43 + 0.2D1 * t44
         t48 = t1 ** 2
         t49 = t4 ** 2
         t50 = t7 ** 2
         t51 = t14 ** 2
         t52 = -t48 - t49 + t50 + t51
         t54 = t41 * t2 + t46 * t5 + t52 * t8
         t57 = -t48 + t49 - t50 + t51
         t60 = 0.2D1 * t43 - 0.2D1 * t44
         t62 = t37 * t2 + t57 * t5 + t60 * t8
         t71 = x(7)
         t73 = x(8)
         t75 = x(9)
         t78 = 0.2D1 * t4 * t71 - 0.2D1 * t1 * t73 - 0.2D1 * t14 * t75
         t83 = t37 * t71 + t57 * t73 + t60 * t75
         t89 = 0.2D1 * t7 * t71 + 0.2D1 * t14 * t73 - 0.2D1 * t1 * t75
         t94 = t41 * t71 + t46 * t73 + t52 * t75
         t101 = 0.2D1 * t1 * t71 + 0.2D1 * t4 * t73 + 0.2D1 * t7 * t75
         t102 = t101 * t54
         t108 = -0.2D1 * t14 * t71 + 0.2D1 * t7 * t73 - 0.2D1 * t4 * t75
         t110 = t94 * t11
         t112 = -t108
         t113 = t112 * t54
         t114 = t83 * t11
         t115 = t101 * t62
         t116 = t94 * t27
         t119 = t83 * t17
         t120 = -t78
         t121 = t120 * t62
         t125 = t48 - t49 - t50 + t51
         t128 = 0.2D1 * t34 - 0.2D1 * t35
         t131 = 0.2D1 * t38 + 0.2D1 * t39
         t133 = t125 * t2 + t128 * t5 + t131 * t8
         t143 = t89 * t133
         t147 = t125 * t71 + t128 * t73 + t131 * t75
         t148 = t147 * t30
         t155 = t101 * t133
         t157 = -t89
         t159 = t147 * t11
         cgret(1,1) = 1D0
         cgret(1,2) = 0D0
         cgret(1,3) = 0D0
         cgret(1,4) = 0D0
         cgret(1,5) = 0D0
         cgret(1,6) = 0D0
         cgret(1,7) = 0D0
         cgret(1,8) = 0D0
         cgret(1,9) = 0D0
         cgret(1,10) = t11
         cgret(1,11) = t17
         cgret(1,12) = t22
         cgret(1,13) = t27
         cgret(2,1) = 0D0
         cgret(2,2) = 1D0
         cgret(2,3) = 0D0
         cgret(2,4) = 0D0
         cgret(2,5) = 0D0
         cgret(2,6) = 0D0
         cgret(2,7) = 0D0
         cgret(2,8) = 0D0
         cgret(2,9) = 0D0
         cgret(2,10) = t29
         cgret(2,11) = t11
         cgret(2,12) = t27
         cgret(2,13) = t30
         cgret(3,1) = 0D0
         cgret(3,2) = 0D0
         cgret(3,3) = 1D0
         cgret(3,4) = 0D0
         cgret(3,5) = 0D0
         cgret(3,6) = 0D0
         cgret(3,7) = 0D0
         cgret(3,8) = 0D0
         cgret(3,9) = 0D0
         cgret(3,10) = t30
         cgret(3,11) = t32
         cgret(3,12) = t11
         cgret(3,13) = t17
         cgret(4,1) = 0D0
         cgret(4,2) = 0D0
         cgret(4,3) = 0D0
         cgret(4,4) = 1D0
         cgret(4,5) = 0D0
         cgret(4,6) = 0D0
         cgret(4,7) = t37 * t54 - t41 * t62
         cgret(4,8) = t57 * t54 - t46 * t62
         cgret(4,9) = t60 * t54 - t52 * t62
         cgret(4,10) = t78 * t54 + t83 * t30 - t89 * t62 - t94 * t29
         cgret(4,11) = t102 + t83 * t32 - t108 * t62 - t110
         cgret(4,12) = t113 + t114 - t115 - t116
         cgret(4,13) = t89 * t54 + t119 - t121 - t94 * t30
         cgret(5,1) = 0D0
         cgret(5,2) = 0D0
         cgret(5,3) = 0D0
         cgret(5,4) = 0D0
         cgret(5,5) = 1D0
         cgret(5,6) = 0D0
         cgret(5,7) = t41 * t133 - t125 * t54
         cgret(5,8) = t46 * t133 - t128 * t54
         cgret(5,9) = t52 * t133 - t131 * t54
         cgret(5,10) = t143 + t110 - t102 - t148
         cgret(5,11) = t108 * t133 + t94 * t17 - t120 * t54 - t147 * t32
         cgret(5,12) = t155 + t94 * t22 - t157 * t54 - t159
         cgret(5,13) = t120 * t133 + t116 - t113 - t147 * t17
         cgret(6,1) = 0D0
         cgret(6,2) = 0D0
         cgret(6,3) = 0D0
         cgret(6,4) = 0D0
         cgret(6,5) = 0D0
         cgret(6,6) = 1D0
         cgret(6,7) = t125 * t62 - t37 * t133
         cgret(6,8) = t128 * t62 - t57 * t133
         cgret(6,9) = t131 * t62 - t60 * t133
         cgret(6,10) = t115 + t147 * t29 - t78 * t133 - t114
         cgret(6,11) = t121 + t159 - t155 - t119
         cgret(6,12) = t157 * t62 + t147 * t27 - t112 * t133 - t83 * t22
         cgret(6,13) = t112 * t62 + t148 - t143 - t83 * t27    
      end

      subroutine dq_x2(i,N,x,p,cgret)
         integer i,N
         double precision x(*)
         double precision p(*)
         double precision cgret(10,13)

         cgret(1,1) = 0D0
         cgret(1,2) = 0D0
         cgret(1,3) = 0D0
         cgret(1,4) = 1D0
         cgret(1,5) = 0D0
         cgret(1,6) = 0D0
         cgret(1,7) = 0D0
         cgret(1,8) = 0D0
         cgret(1,9) = 0D0
         cgret(1,10) = 0D0
         cgret(1,11) = 0D0
         cgret(1,12) = 0D0
         cgret(1,13) = 0D0
         cgret(2,1) = 0D0
         cgret(2,2) = 0D0
         cgret(2,3) = 0D0
         cgret(2,4) = 0D0
         cgret(2,5) = 1D0
         cgret(2,6) = 0D0
         cgret(2,7) = 0D0
         cgret(2,8) = 0D0
         cgret(2,9) = 0D0
         cgret(2,10) = 0D0
         cgret(2,11) = 0D0
         cgret(2,12) = 0D0
         cgret(2,13) = 0D0
         cgret(3,1) = 0D0
         cgret(3,2) = 0D0
         cgret(3,3) = 0D0
         cgret(3,4) = 0D0
         cgret(3,5) = 0D0
         cgret(3,6) = 1D0
         cgret(3,7) = 0D0
         cgret(3,8) = 0D0
         cgret(3,9) = 0D0
         cgret(3,10) = 0D0
         cgret(3,11) = 0D0
         cgret(3,12) = 0D0
         cgret(3,13) = 0D0
         cgret(4,1) = 0D0
         cgret(4,2) = 0D0
         cgret(4,3) = 0D0
         cgret(4,4) = 0D0
         cgret(4,5) = 0D0
         cgret(4,6) = 0D0
         cgret(4,7) = 1D0
         cgret(4,8) = 0D0
         cgret(4,9) = 0D0
         cgret(4,10) = 0D0
         cgret(4,11) = 0D0
         cgret(4,12) = 0D0
         cgret(4,13) = 0D0
         cgret(5,1) = 0D0
         cgret(5,2) = 0D0
         cgret(5,3) = 0D0
         cgret(5,4) = 0D0
         cgret(5,5) = 0D0
         cgret(5,6) = 0D0
         cgret(5,7) = 0D0
         cgret(5,8) = 1D0
         cgret(5,9) = 0D0
         cgret(5,10) = 0D0
         cgret(5,11) = 0D0
         cgret(5,12) = 0D0
         cgret(5,13) = 0D0
         cgret(6,1) = 0D0
         cgret(6,2) = 0D0
         cgret(6,3) = 0D0
         cgret(6,4) = 0D0
         cgret(6,5) = 0D0
         cgret(6,6) = 0D0
         cgret(6,7) = 0D0
         cgret(6,8) = 0D0
         cgret(6,9) = 1D0
         cgret(6,10) = 0D0
         cgret(6,11) = 0D0
         cgret(6,12) = 0D0
         cgret(6,13) = 0D0
         cgret(7,1) = 0D0
         cgret(7,2) = 0D0
         cgret(7,3) = 0D0
         cgret(7,4) = 0D0
         cgret(7,5) = 0D0
         cgret(7,6) = 0D0
         cgret(7,7) = 0D0
         cgret(7,8) = 0D0
         cgret(7,9) = 0D0
         cgret(7,10) = 30D0
         cgret(7,11) = 0D0
         cgret(7,12) = 0D0
         cgret(7,13) = 0D0
         cgret(8,1) = 0D0
         cgret(8,2) = 0D0
         cgret(8,3) = 0D0
         cgret(8,4) = 0D0
         cgret(8,5) = 0D0
         cgret(8,6) = 0D0
         cgret(8,7) = 0D0
         cgret(8,8) = 0D0
         cgret(8,9) = 0D0
         cgret(8,10) = 0D0
         cgret(8,11) = 30D0
         cgret(8,12) = 0D0
         cgret(8,13) = 0D0
         cgret(9,1) = 0D0
         cgret(9,2) = 0D0
         cgret(9,3) = 0D0
         cgret(9,4) = 0D0
         cgret(9,5) = 0D0
         cgret(9,6) = 0D0
         cgret(9,7) = 0D0
         cgret(9,8) = 0D0
         cgret(9,9) = 0D0
         cgret(9,10) = 0D0
         cgret(9,11) = 0D0
         cgret(9,12) = 30D0
         cgret(9,13) = 0D0
         cgret(10,1) = 0D0
         cgret(10,2) = 0D0
         cgret(10,3) = 0D0
         cgret(10,4) = 0D0
         cgret(10,5) = 0D0
         cgret(10,6) = 0D0
         cgret(10,7) = 0D0
         cgret(10,8) = 0D0
         cgret(10,9) = 0D0
         cgret(10,10) = 0D0
         cgret(10,11) = 0D0
         cgret(10,12) = 0D0
         cgret(10,13) = 30D0
      end

      subroutine nq_x(N,p,res)
         integer N,res
         double precision p
         dimension res(N+1)
         res=0
         res(1)=10
         res(N+1)=6
      end

      subroutine r_x(i,N,x,p,res)
         integer i,N
         double precision x,p,res

         dimension x(13)
         dimension p(*)
         dimension res(*)

         call r_xsub(i,N,x,p,p(13),res)

      end

      subroutine r_xsub(i,N,x,p,target,res)
         integer i,N
         double precision x,p,res,target

         dimension x(13)
         dimension target(3,N+1)
         dimension p(*)
         dimension res(*)

         res(1)= (x(1)-target(1,i))**2 + (x(2)-target(2,i))**2 + 
     \     (x(3)-target(3,i))**2 - p(10)**2;
      end

      subroutine dr_x(i,N,x,p,res)
         integer i,N
         double precision x,p,res

         dimension x(13)
         dimension p(*)
         dimension res(*)

         call dr_xsub(i,N,x,p(13),res)

      end

      subroutine dr_xsub(i,N,x,target,res)
         integer i,N
         double precision x,res,target

         dimension x(13)
         dimension target(3,N+1)
         dimension res(13)

         res(1)= 2D0*(x(1)-target(1,i))
         res(2)= 2D0*(x(2)-target(2,i))
         res(3)= 2D0*(x(3)-target(3,i))
         res(4)=0D0
         res(5)=0D0
         res(6)=0D0
         res(7)=0D0
         res(8)=0D0
         res(9)=0D0
         res(10)=0D0
         res(11)=0D0
         res(12)=0D0
         res(13)=0D0
      end

      subroutine nr_x(N,p,res)
         integer N,res
         double precision p
         dimension p(*)
         dimension res(N+1)
         res=1
      end

      subroutine nbox_u(N,p,lres,ures)
         integer N,lres,ures
         double precision p
         dimension lres(N),ures(N)
         dimension p(*)
         lres=6
         ures=6
      end

      subroutine Pbox_u(i,N,p,lres,ures)
         integer i,N,lres,ures
         double precision p
         dimension lres(6),ures(6)
         dimension p(*)
         lres(1)=1
         lres(2)=2
         lres(3)=3
         lres(4)=4
         lres(5)=5
         lres(6)=6
         ures(1)=1
         ures(2)=2
         ures(3)=3
         ures(4)=4
         ures(5)=5
         ures(6)=6
      end

      subroutine box_u(i,N,p,lb,ub)
         integer i,N
         double precision p,lb,ub
         dimension lb(6),ub(6)
         dimension p(*)
         lb(1)=-p(11)
         lb(2)=-p(11)
         lb(3)=-p(11)
         lb(4)=-p(12)
         lb(5)=-p(12)
         lb(6)=-p(12)
         ub(1)=p(11)
         ub(2)=p(11)
         ub(3)=p(11)
         ub(4)=p(12)
         ub(5)=p(12)
         ub(6)=p(12)
      end

      subroutine nbox_x(N,p,lres,ures)
         integer N,lres,ures
         double precision p
         dimension lres(N+1),ures(N+1)
         dimension p(*)
         lres=0
         ures=0
         lres(1)=3
         ures(1)=3
      end

      subroutine Pbox_x(i,N,p,lres,ures)
         integer i,N,lres,ures
         double precision p
         dimension lres(*),ures(*)
         dimension p(*)

         if (i.eq.1) then
             lres(1)=1
             lres(2)=2
             lres(3)=3
             ures(1)=1
             ures(2)=2
             ures(3)=3
        end if
      end

      subroutine box_x(i,N,p,lb,ub)
         integer i,N
         double precision p,lb,ub
         dimension lb(*),ub(*)
         dimension p(*)

         if (i.eq.1) then
             lb(1)=-1D2
             lb(2)=-1D2
             lb(3)=-1D2
             ub(1)=1D2
             ub(2)=1D2
             ub(3)=1D2
         end if
      end

      subroutine readnx (cgret)
        integer cgret(1)
        cgret(1) = 13
      end

      subroutine readnu (cgret)
        integer cgret(1)
        cgret(1) = 6
      end

      subroutine rdrks(s)
      integer s
      s=3
      end

      subroutine readrk(rkA,rkb,rkc)
      double precision rkA,rkb,rkc
      dimension rkA(3,3),rkb(3),rkc(3)
      rkA=0D0
      rkb=0D0
      rkc=0D0
      rkA(1,1) = 1D0/9D0
      rkA(1,2) = (-1D0-sqrt(6D0))/18D0
      rkA(1,3) = (-1D0+sqrt(6D0))/18D0
      rkA(2,1) = 1D0/9D0
      rkA(2,2) = (88D0+7D0*sqrt(6D0))/360D0
      rkA(2,3) = (88D0-43D0*sqrt(6D0))/360D0
      rkA(3,1) = 1D0/9D0
      rkA(3,2) = (88D0+43D0*sqrt(6D0))/360D0
      rkA(3,3) = (88D0-7D0*sqrt(6D0))/360D0
      rkb(1) = 1D0/9D0
      rkb(2) = (16D0+sqrt(6D0))/36D0
      rkb(3) = (16D0-sqrt(6D0))/36D0
      rkc(2) = (6D0-sqrt(6D0))/10D0
      rkc(3) = (6D0+sqrt(6D0))/10D0
      end

      subroutine rotate(x,res)
      double precision x
      double precision res
      dimension x(13)
      dimension res(3,3)

      double precision t1
      double precision t2
      double precision t3
      double precision t4
      double precision t5
      double precision t6
      double precision t7
      double precision t8
      double precision t9
      double precision t10
      double precision t11
      double precision t12
      double precision t13
      double precision t14

      t1 = x(10)
      t2 = x(11)
      t3 = x(12)
      t4 = x(13)
      t5 = t1**2
      t6 = t2**2
      t7 = t3**2
      t8 = t4**2
      t9 = t1*t2
      t10 = t3*t4
      t11 = t1*t3
      t12 = t2*t4
      t13 = t2*t3
      t14 = t1*t4


      res(1,1) = t5-t6-t7+t8
      res(2,1) = 2D0*(t9 + t10)
      res(3,1) = 2D0*(t11 - t12)
      res(1,2) = 2D0*(t9 - t10)
      res(2,2) = -t5+t6-t7+t8
      res(3,2) = 2D0*(t13 + t14)
      res(1,3) = 2D0*(t11 + t12)
      res(2,3) = 2D0*(t13 - t14)
      res(3,3) =  -t5-t6+t7+t8

      end