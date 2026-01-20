       program Labinit
       implicit none
       integer N, Lo, Lp, Lm, j
       integer idum, iy, ir(97)
       double precision nu, kT, ks
       double precision g, alpha, Tdif
       double precision H, lambda, Uo, f
       double precision C2, ZC2, chi
       double precision dt, tmax
       double precision bT, Pr
       double precision Up, fp
       double precision z, u, v, w, T
       double precision small, random
 10    format(3i12)
 20    format(3d15.4)
 30    format(5d16.8)
       small = 2.d-10
       nu = 1.41d-5
       kT = 1.96d-5
       g = 9.81d0
       alpha = 3.5d-3
       Tdif = 12.d0
       H = 1.d0
       Uo = 0.d0
       f = 0.d0
       C2 = 1.5d3
       ZC2 = 1.d5
       chi = 0.d0
       N = 6000
       Lo = 25
       Lp = 100
       Lm = (N/3)
       dt = 1.d3/(1.d0*N*N)
c       dt = 0.00015
       tmax = 1.d-3
       bT = (8.d0*g*alpha*Tdif*C2*H*H*H)/(27.d0*nu*nu)
       Pr = nu/kT
       Up = dsqrt(C2)*Uo*H/nu
       fp = f*H*H/nu
       open(1,file="LabExppar.dat", status="unknown")
       write(1,10) N
       write(1,20) bT, Pr
       write(1,20) Up, fp
       write(1,20) dt, tmax
       write(1,20) ZC2, chi
       write(1,10) Lo, Lp, Lm
       write(1,10) -1234
       close(1)
       idum = -5555
       open(1, file="U.dat", status="unknown")
       open(2, file="V.dat", status="unknown")
       open(3, file="W.dat", status="unknown")
       open(4, file="T.dat", status="unknown")
       do j=1, N
        z = (1.d0*j)/(1.d0*N)
        if (j .eq. N) small = 0.d0
        u = small*(random(idum,iy,ir) - 0.5d0) + (z*Up)
        v = small*(random(idum,iy,ir) - 0.5d0)
        w = small*(random(idum,iy,ir) - 0.5d0)
        T = z
        write(1, 30) z, u
        write(2, 30) z, v
        write(3, 30) z, w
        write(4, 30) z, T
       enddo
       close(1)
       close(2)
       close(3)
       close(4)
       stop
       end
       
       
       function random(idum,iy,ir)
       implicit none
       integer idum, iy, m, ia, ic, j
       double precision random, rm
       integer ir(97)
       m = 714025
       ia = 1366
       ic = 150889
       rm = 1.4005112d-6
       if (idum .lt. 0) then
        idum = mod(ic-idum,m)
        do j=1,97
         idum = mod(ia*idum+ic,m)
         ir(j) = idum
        enddo
        idum = mod(ia*idum+ic,m)
        iy = idum
       endif
       j = 1 + (97*iy)/m
       if (j .gt. 97 .or. j .lt. 1) then
        j = mod(j,97)
        if (j .ge. 0) j=-j
        j = j + 1
       endif
       iy = ir(j)
       idum = mod(ia*idum+ic,m)
       ir(j) = idum
       random = iy*rm
       return
       end 