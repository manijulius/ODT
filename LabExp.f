       program LabExp
       implicit none
       integer M, L, idum, iy, Length
       integer N, Lo, Lp, Lm
       integer NL(100000), ir(97), num_args, realz_int, j
       character(len=32), dimension(:), allocatable :: args
       character(len=:), allocatable :: input_dir,cmd_str,exp_name
       character(len=:), allocatable :: realization
       double precision Nt, Np, Na, Nd, Co, Cm, uK, vK, wK, PE
       double precision pmax, p, prob, pa
       double precision dt, time, td, tmax, te, ts, random
       double precision u(100000), ua(100000), ur(100000), eu(100000)
       double precision v(100000), va(100000), vr(100000), ev(100000)
       double precision w(100000), wa(100000), wr(100000), ew(100000)
       double precision T(100000), Ta(100000), Tr(100000), eT(100000)
       double precision PL(100000), a(14), pdf(5,1000)

c      Command line arguments      
      num_args = command_argument_count()
      allocate(args(num_args))
      do j = 1, num_args
        call get_command_argument(j,args(j))
      end do

c      get input_dir, experiment name and realization value     
      if (num_args   < 3)then
          write(*,*) "Incorrect number of arguments"
          return
      else
          input_dir      = trim(args(1))
          exp_name       = trim(args(2)) 
          realization    = trim(args(3)) 
      endif     
      read(realization,*) realz_int
      write(*,'(6A)') "input:", input_dir," exp:",exp_name,
     +     " realization:", realization

c      change to input directory 
      cmd_str           = "input/"//input_dir
      call chdir(cmd_str)
      call system("pwd")

10     format(5g16.6)
20     format(g20.1,3g12.3,2i8)
       pmax = 0.1d0
       Nd   = 0.d0

       call readpar(idum,N,Lo,Lp,Lm,dt,td,tmax,a)
       idum   = idum + realz_int;
       write(*,"(A,I8)") "Random seed Idum:", idum

       call init(N,u,v,w,T)
       call zeroparam(NL,Nt,Np,Na,pdf,time,te,ts)
       call zerovars(N,ua,ur,eu,va,vr,ev,wa,wr,ew,Ta,Tr,eT)
       call LenProb(Lo,Lm,PL,a(14),Co,Cm)

c      change to parent directory       
       call chdir("../../") 
       write(*,*) "now at parent dir:"
       call system("pwd");
      
       cmd_str   = "mkdir -p output/"//exp_name//"/"//realization
       write(*,*) "cmd: ",cmd_str
       call system(cmd_str)
       cmd_str   = "output/"//exp_name//"/"//realization
c      excecuting cd did not work. so using chdir       
       call chdir(cmd_str)
       call system("pwd")

       open(101, file="Warnings.dat", status="unknown")
       open(102, file="EddyDiagram.dat",status="unknown")

       do while (time .le. tmax)
        time = time + dt
        Nt = Nt + 1.d0
c        write(*,*) "Nt:", Nt, " time:", time
c       track time and time step
        if (mod(Nt,1000000.0) == 0) then
           write(*,'(A24, F12.0, G12.2)') "Nt, time:",        Nt, time
           write(*,'(A24, 4F12.0)') "Nt, Np, Nd, Na:",  Nt, Np, Nd, Na
           write(*,'(A24, 2G12.2)') "tmax, dt:",        tmax, dt
           write(*,*)
           
c         exit
        endif

        if ((time-te) .ge. td) then
         call vis(N,u,v,w,time-te,a(9),a(10))
         call difT(N,T,time-te,a(4))
         call statsvel(N,u,ua,ur,eu,v,va,vr,ev,w,wa,wr,ew,time-te)
         call statsdens(N,T,Ta,Tr,eT,time-te)
         call comppdf(N,t,pdf,time-te)
         Nd = Nd + 1.d0
         te = time
        endif

        L = 3*Length(PL,a(14),Co,Cm,idum,iy,ir)
        M = 1 + int(random(idum,iy,ir)*(N-L))
        p = dt*prob(N,M,L,u,v,w,T,a,uK,vK,wK,PE)
        if (p .gt. 0.d0) then
         if (p .gt. pmax) then
          write(101,20) "Time step warning: ", time, dt, p, M, L
          dt = dt*pmax/p
          p = pmax
         endif
         pa = pa + p
         Np = Np + 1.d0
        endif

        if (random(idum,iy,ir) .lt. p) then 
c         call energy(N,M,M+L,a(2),u,v,w,T)
         call eddy(N,M,L,u,v,w,T,uK,vK,wK,PE,a(12))
c         call energy(N,M,M+L,a(2),u,v,w,T)
         call vis(N,u,v,w,time-te,a(9),a(10))
         call difT(N,T,time-te,a(4))
         call statsvel(N,u,ua,ur,eu,v,va,vr,ev,w,wa,wr,ew,time-te)
         call statsdens(N,T,Ta,Tr,eT,time-te)
         call comppdf(N,t,pdf,time-te)
         write(102,10) time, 1.d0*(M+(L/2))/(1.d0*N), (1.d0*L)/(2.d0*N)
         Na = Na + 1.d0
         NL(L/3) = NL(L/3) + 1
         te = time
        endif

        if (Np > 1.d6) then 
          call raisedt(Np,dt,pa)
        endif

       enddo

       call vis(N,u,v,w,time-te,a(9),a(10))
       call difT(N,T,time-te,a(4))
       call statsvel(N,u,ua,ur,eu,v,va,vr,ev,w,wa,wr,ew,time-te)
       call statsdens(N,T,Ta,Tr,eT,time-te)
       call comppdf(N,t,pdf,time-te)
       close(101)
       close(102)

       call outputstats(N,Na,Nt,Lo,Lm,NL,PL,pdf,time,dt)
       call saveveldata(N,u,ua,ur,eu,v,va,vr,ev,w,wa,wr,ew,time)
       call savedensdata(N,T,Ta,Tr,eT,time)
       
       stop
       end
    
       function prob(N,M,L,u,v,w,T,a,uK,vK,wK,PE)
       implicit none
       integer N, M, L
       double precision prob, u(100000), v(100000), w(100000)
       double precision T(100000), a(14), uK, vK, wK, PE, KE
       double precision TK, p, x, psiK
       x = (1.d0*L)/(1.d0*N)
       uK = psiK(N,M,L,u)
       vK = psiK(N,M,L,v)
       wK = psiK(N,M,L,w)
       TK = -psiK(N,M,L,T)
       PE = a(2)*TK*x
       KE = ((1.d0 - a(12))*wK*wK) + (a(12)*((uK*uK) + (vK*vK))/2.d0)
       p = ((KE + PE)*x*x) - a(11)

       if (p .gt. 0.d0) then 
        prob = a(1)*(1.d0-x)*dsqrt(p)*dexp(3.d0*a(14)/(x*N))/(x*x)
       else 
        prob = 0.d0
       endif

       return
       end

       function psiK(N,M,L,psi)
       implicit none
       integer N, M, L, j
       double precision psiK, sum, z
       double precision psi(100000)
       sum = psi(M)*L/2.d0
       do j=1, L-1
        z = L -( 2*j)
        sum = sum + (psi(j+M)*z)
       enddo
       sum = sum - (psi(M+L)*L/2.d0)
       psiK = 4.d0*sum/(9.d0*L*L)
       return
       end
       
       subroutine eddy(N,M,L,u,v,w,T,uK,vK,wK,PE,c)
       implicit none
       integer N, M, L
       double precision u(100000), v(100000), w(100000)
       double precision T(100000), uK, vK, wK, PE, c
       double precision qu, qv, qw, cu, cv, cw
       qu = dsqrt((c*((vK*vK) + (wK*wK))/2.d0) + ((1.d0-c)*uK*uK))
       if (uK .gt. 0.d0) then
        cu = 6.75d0*(qu - uK)/(1.d0*L)
       else
        cu = -6.75d0*(qu + uK)/(1.d0*L)
       endif
       qv = dsqrt((c*((uK*uK) + (wK*wK))/2.d0) + ((1.d0-c)*vK*vK))
       if (vK .gt. 0.d0) then
        cv = 6.75d0*(qv - vK)/(1.d0*L)
       else
        cv = -6.75d0*(qv + vK)/(1.d0*L)
       endif
       qw = dsqrt((c*((uK*uK) + (vK*vK))/2.d0) + ((1.d0-c)*wK*wK) + PE)
       if (wK .gt. 0.d0) then
        cw = 6.75d0*(qw - wK)/(1.d0*L)
       else
        cw = -6.75d0*(qw + wK)/(1.d0*L)
       endif
       call triplet(N,M,L,u)
       call triplet(N,M,L,v)
       call triplet(N,M,L,w)
       call triplet(N,M,L,T)
       call addK(N,M,L,u,cu)
       call addK(N,M,L,v,cv)
       call addK(N,M,L,w,cw)
       return
       end

       subroutine triplet(N,M,L,psi)
       implicit none
       integer N, M, L, Lo, j, k
       double precision psi(100000), x(100000)
       Lo = L/3
       do j = 1, Lo
        k = M + 3*(j-1)
        x(j) = psi(k)
       enddo
       do j=1, Lo
        k = M + L + 1 - (3*j)
        x(j+Lo) = psi(k)
       enddo
       do j = 1, Lo
        k = M + (3*j) - 1
        x(j+Lo+Lo) = psi(k)
       enddo
       do j=1, L
        k = M+j-1
        psi(k) = x(j)
       enddo
       return
       end
       
       subroutine addK(N,M,L,u,c)
       implicit none
       integer N, M, L, Lo, j, j1, j2, j3
       double precision c, y1, y2, y3
       double precision u(100000)
       Lo = L/3
       do j = 1, Lo
        y1 = -(2.d0*j)
        y2 = (4.d0*(j+Lo)) - (2.d0*L)
        y3 = (2.d0*L) - (2.d0*(j+Lo+Lo))
        j1 = M + j
        j2 = M + j + Lo
        j3 = M + j + Lo + Lo
        u(j1) = u(j1) + (c*y1)
        u(j2) = u(j2) + (c*y2)
        u(j3) = u(j3) + (c*y3)
       enddo
       return
       end
       
       subroutine difT(N,T,dt,Pr)
       implicit none
       integer N, j
       double precision T(100000), x(100000)
       double precision l(100000), d(100000), r(100000)
       double precision dt, Pr, De
       De = (dt*N*N)/(2.d0*Pr)
       l(1) = 0.d0
       d(1) = 1.d0 + (2.d0*De)
       r(1) = -De
       do j = 2, N-1
         l(j) = -De
         d(j) = 1.d0 + (2.d0*De)
         r(j) = -De
       enddo
       d(N) = 1.d0
       l(N) = 0.d0
       r(N) = 0.d0
       x(1) = ((1.d0 - (2.d0*De))*T(1)) + (De*T(2))
       do j = 2, N-1
         x(j) = ((1.d0-(2.d0*De))*T(j)) + (De*(T(j+1) + T(j-1)))
       enddo
       x(N) = T(N)
       call tridiagonal(N,l,d,r,x,T)
       return
       end
       
       subroutine vis(N,u,v,w,dt,Uo,f)
       implicit none
       integer N, j
       double precision u(100000), v(100000), w(100000), dt, Uo, f, De
       double precision l(100000), d(100000), r(100000)
       double precision xu(100000), xv(100000), xw(100000)
       double precision duf(100000), dvf(100000), cft, sft
       De = (dt*N*N)/(2.d0)
       l(1) = 0.d0
       d(1) = 1.d0 + (2.d0*De)
       r(1) = -De
       do j = 2, N-1
         l(j) = -De
         d(j) = 1.d0 + (2.d0*De)
         r(j) = -De
       enddo
       d(N) = 1.d0
       l(N) = 0.d0
       r(N) = 0.d0
       cft = dcos(f*dt) - 1.d0
       sft = dsin(f*dt)
       do j = 1, N-1
         duf(j) = ((u(j) - Uo)*cft) + (v(j)*sft)
         dvf(j) = (v(j)*cft) - ((u(j) - Uo)*sft)
       enddo
       xu(1) = ((1.d0-(2.d0*De))*u(1)) + (De*u(2))
       xv(1) = ((1.d0-(2.d0*De))*v(1)) + (De*v(2))
       xw(1) = ((1.d0-(2.d0*De))*w(1)) + (De*w(2))
       do j = 2, N-1
         xu(j) = ((1.d0-(2.d0*De))*u(j)) + (De*(u(j+1)+u(j-1)))
         xv(j) = ((1.d0-(2.d0*De))*v(j)) + (De*(v(j+1)+v(j-1)))
         xw(j) = ((1.d0-(2.d0*De))*w(j)) + (De*(w(j+1)+w(j-1)))
       enddo
       xu(N) = u(N)
       xv(N) = v(N)
       xw(N) = w(N)
       call tridiagonal(N,l,d,r,xu,u)
       call tridiagonal(N,l,d,r,xv,v)
       call tridiagonal(N,l,d,r,xw,w)
       do j = 1, N-1
         u(j) = u(j) + duf(j)
         v(j) = v(j) + dvf(j)
       enddo
       return
       end
       
       subroutine tridiagonal(N,l,d,r,x,y)  
       implicit none
       integer N, j
       double precision l(100000), d(100000), r(100000)
       double precision x(100000), y(100000), b, g(100000)
       b = d(1)
       y(1) = x(1)/b
       do j = 2, N
         g(j) = r(j-1)/b
         b = d(j) - (l(j)*g(j))
         if (b .eq. 0.d0) pause 'Tridiagonal:  failure'
         y(j) = (x(j) - (l(j)*y(j-1)))/b
       enddo
       do j = N-1, 1, -1
         y(j) = y(j) - (g(j+1)*y(j+1))
       enddo
       return 
       end
       
       subroutine raisedt(Np,dt,p)
       implicit none
       double precision Np, dt, p, pmin
       pmin = 1.d-3
       p = p/Np
       if (p .lt. (pmin/2.d0)) then 
         dt = dt*2.d0;
       else
         dt = dt*pmin/p
       endif
       p = 0.d0
       Np = 0.d0
       return
       end
   
       function Length(PL,xp,Co,Cm,idum,iy,ir)
       implicit none
       integer Length, idum, iy, n
       integer ir(97)
       double precision xp, Co, Cm, r, random, x
       double precision PL(100000)
       r = random(idum,iy,ir)
       x = -xp/dlog((Co*r)+(Cm*(1.d0-r)))
       n = int(x)-1
       if (r .gt. PL(n)) then
 10     n = n + 1
        if (r .gt. PL(n)) goto 10
       endif
       if (r .lt. PL(n-1)) then
 20     n = n - 1
        if (r .lt. PL(n-1)) goto 20
       endif
       Length = n
       return
       end

       subroutine LenProb(Lo,Lm,P,xp,Co,Cm)
       implicit none
       integer Lo, Lm, L
       double precision P(100000), xp, Co, Cm, C, z
       Co = dexp(-xp/(1.d0*Lo))
       Cm = dexp(-xp/(1.d0*Lm))
       C = 0.d0
       do L = Lo, Lm
        z = dexp(-xp/(1.d0*L))*(dexp(xp/(L*(L+1.d0)))-1.d0)
        C = C + z
       enddo
       C = 1.d0/C
       do L = 1, Lo-1
        P(L) = 0.d0
       enddo
       do L = Lo, Lm
        z = dexp(-xp/(1.d0*L))*(dexp(xp/(L*(L+1.d0)))-1.d0)
        P(L) = P(L-1) + (C*z)
       enddo
       do L = Lm+1, 100000
        P(L) = 0.d0
       enddo
       return
       end
       
       subroutine comppdf(N,T,pdf,dt)
       implicit none
       integer N, j, k(5), i, w
       double precision dt, T(100000), pdf(5,1000), X
c removed '-' in X
       X = 1.d3                                                     
       w = 5
       k(1) = N/2
       k(2) = N/4
       k(3) = N/8
       k(4) = 5*w
       k(5) = 2*w
       do j = -w, w
        i = int(X*(T(j+k(1))))
        pdf(1,i) = pdf(1,i) + dt
        i = int(X*(T(j+k(2))))
        pdf(2,i) = pdf(2,i) + dt
        i = int(X*(T(j+k(3))))
        pdf(3,i) = pdf(3,i) + dt
        i = int(X*(T(j+k(4))))
        pdf(4,i) = pdf(4,i) + dt
        i = int(X*(T(j+k(5))))
        pdf(5,i) = pdf(5,i) + dt
       enddo
       return
       end
   
       subroutine statsvel(N,u,ua,ur,eu,v,va,vr,ev,w,wa,wr,ew,dt)
       implicit none
       integer N, j
       double precision u(100000), ua(100000), ur(100000), eu(100000)
       double precision v(100000), va(100000), vr(100000), ev(100000)
       double precision w(100000), wa(100000), wr(100000), ew(100000)
       double precision dt, dudz, dvdz, dwdz
       do j = 1, N
        ua(j) = ua(j)+(u(j)*dt)
        va(j) = va(j)+(v(j)*dt)
        wa(j) = wa(j)+(w(j)*dt)
        ur(j) = ur(j)+(u(j)*u(j)*dt)
        vr(j) = vr(j)+(v(j)*v(j)*dt)
        wr(j) = wr(j)+(w(j)*w(j)*dt)
       enddo
       dudz = (u(2) - u(1))*N
       dvdz = (v(2) - v(1))*N
       dwdz = (w(2) - w(1))*N
       eu(1) = eu(1)+(dudz*dudz*dt)
       ev(1) = ev(1)+(dvdz*dvdz*dt)
       ew(1) = ew(1)+(dwdz*dwdz*dt)
       do j = 2, N-1
         dudz = (u(j+1) - u(j-1))*N/2.d0
         dvdz = (v(j+1) - v(j-1))*N/2.d0
         dwdz = (w(j+1) - w(j-1))*N/2.d0
         eu(j) = eu(j)+(dudz*dudz*dt)
         ev(j) = ev(j)+(dvdz*dvdz*dt)
         ew(j) = ew(j)+(dwdz*dwdz*dt)
       enddo
       dudz = (u(N) - u(N-1))*N
       dvdz = (v(N) - v(N-1))*N
       dwdz = (w(N) - w(N-1))*N
       eu(N) = eu(N)+(dudz*dudz*dt)
       ev(N) = ev(N)+(dvdz*dvdz*dt)
       ew(N) = ew(N)+(dwdz*dwdz*dt)
       return
       end
       
       subroutine statsdens(N,T,Ta,Tr,eT,dt)
       implicit none
       integer N, j
       double precision T(100000), Ta(100000), Tr(100000), eT(100000)
       double precision dt, delT
       do j = 1, N
        Ta(j) = Ta(j)+(T(j)*dt)
        Tr(j) = Tr(j)+(T(j)*T(j)*dt)
       enddo
       delT = T(2) - T(1)
       eT(1) = eT(1)+(delT*delT*N*N*dt)
       do j = 2, N-1
         delT = T(j+1) - T(j-1)
         eT(j) = eT(j)+(delT*delT*N*N*dt/4.d0)
       enddo
       delT = T(N) - T(N-1)
       eT(N) = eT(N)+(delT*delT*N*N*dt)
       return
       end
       
       subroutine saveveldata(N,u,ua,ur,eu,v,va,vr,ev,w,wa,wr,ew,tf)
       implicit none
       integer N, j, k
       double precision u(100000), ua(100000), ur(100000), eu(100000)
       double precision v(100000), va(100000), vr(100000), ev(100000)
       double precision w(100000), wa(100000), wr(100000), ew(100000)
       double precision tf, z, small, dis, Re, etot
10     format(5g16.6)
20     format(g30.1,g15.3)
       small = 1.d-10
       dis = 0.d0
       do j = 1, N
        ua(j) = ua(j)/tf
        va(j) = va(j)/tf
        wa(j) = wa(j)/tf
        ur(j) = (ur(j)/tf)-(ua(j)*ua(j))
        vr(j) = (vr(j)/tf)-(va(j)*va(j))
        wr(j) = (wr(j)/tf)-(wa(j)*wa(j))
        eu(j) = eu(j)/tf
        ev(j) = ev(j)/tf
        ew(j) = ew(j)/tf
        dis = dis + eu(j) + ev(j) + ew(j)
       enddo
       dis = dis/(1.d0*N)
       Re = 0.d0
       k = 3*N/8
       do j = 1, N/4
        Re = Re + ur(j+k) + vr(j+k) + wr(j+k)
       enddo
       Re = dsqrt(4.d0*Re/(1.d0*N))
       write(6,20) "Core Reynolds Number, ReC: ", Re
       write(6,20) "Total KE dissipation:      ", dis
       write(6,20) "Wall stress, u:            ", ua(1)*N
       write(6,20) "                           ", (ua(N)-ua(N-1))*N
       write(6,20) "Wall stress, v:            ", va(1)*N
       write(6,20) "                           ", (va(N)-va(N-1))*N
       write(6,20) "Wall stress, w:            ", wa(1)*N
       write(6,20) "                           ", (wa(N)-wa(N-1))*N
       write(6,20)
       open(100, file="U.dat")
       open(200, file="V.dat")
       open(300, file="W.dat")
       open(400, file="eu.dat")
       do j = 1, N
        z = (1.d0*j)/(1.d0*N)
        etot = eu(j) + ev(j) + ew(j)
        if (abs(ua(j)) .lt. small) ua(j) = 0.d0
        if (abs(va(j)) .lt. small) va(j) = 0.d0
        if (abs(wa(j)) .lt. small) wa(j) = 0.d0
        if (ur(j) .lt. 0.d0) ur(j) = 0.d0
        if (vr(j) .lt. 0.d0) vr(j) = 0.d0
        if (wr(j) .lt. 0.d0) wr(j) = 0.d0
        write(100,10) z, u(j), ua(j), dsqrt(ur(j))
        write(200,10) z, v(j), va(j), dsqrt(vr(j))
        write(300,10) z, w(j), wa(j), dsqrt(wr(j))
        write(400,10) z, eu(j), ev(j), ew(j), etot
       enddo
       close(100)
       close(200)
       close(300)
       close(400)
       return
       end
       
       subroutine savedensdata(N,T,Ta,Tr,eT,tf)
       implicit none
       integer N, j, k
       double precision T(100000), Ta(100000), Tr(100000), eT(100000)
       double precision tf, z, small, dis, x
10     format(5g16.6)
20     format(g30.1,g15.4)
       small = 1.d-10
       dis = 0.d0
       do j = 1, N
        Ta(j) = Ta(j)/tf
        Tr(j) = (Tr(j)/tf)-(Ta(j)*Ta(j))
        eT(j) = eT(j)/tf
        dis = dis + eT(j)
       enddo
       dis = dis/(1.d0*N)
       x = 0.d0
       k = 3*(N/8)
       do j = 1, N/4
        x = x + Tr(j+k)
       enddo
       x = dsqrt(4.d0*x/(1.d0*N))
       write(6,20) "Core Temperature Fluct.:   ", x
       write(6,20) "Total Temperature dissip.: ", dis
       write(6,20) "Nusselt Number:            ", Ta(1)*N
       write(6,20) "                           ", (Ta(N)-Ta(N-1))*N
       open(100, file="T.dat", status="unknown")
       do j = 1, N
        z = (1.d0*j)/(1.d0*N)
        if (abs(Ta(j)) .lt. small) Ta(j) = 0.d0
        if (Tr(j) .lt. 0.d0) Tr(j) = 0.d0
        write(100,10) z, T(j), Ta(j), dsqrt(Tr(j))
       enddo
       close(100)
       return
       end

       subroutine outputstats(N,Na,Nt,Lo,Lm,NL,PL,pdf,tf,dt)
       implicit none
       integer N, NL(100000), Lo, Lm, j
       double precision Na, Nt, tf, dt, x, xo
       double precision PL(100000), pdf(5,1000)
 10    format(6g16.8)
 20    format(g30.1,g15.4) 
       open(1, file="PL.dat", status="unknown")
       do j = Lo, Lm
        if (NL(j) .gt. 0) then
         x = (1.d0*NL(j))/Na
         xo = PL(j) - PL(j-1)
         write(1,10) dlog((3.d0*j)/(1.d0*N)), dlog(x), dlog(xo)
        endif
       enddo
       close(1)
       write(6,20) "Final Time Step:           ", dt
       write(6,20) "Average Eddy Time Step:    ", tf/Nt
       write(6,20) "Time Between Eddies:       ", tf/Na
       write(6,20) "Total Number of Eddies:    ", Na
       write(6,20) "Eddie Acceptance Rate:     ", Na/Nt
       write(6,20)
       open(1, file="PDF.dat", status="unknown")
       do j = 1 , 1000
        x = -j/1.d3
        pdf(1,j) = pdf(1,j)/tf
        pdf(2,j) = pdf(2,j)/tf
        pdf(3,j) = pdf(3,j)/tf
        pdf(4,j) = pdf(4,j)/tf
        pdf(5,j) = pdf(5,j)/tf
        write(1,10) x, pdf(1,j), pdf(2,j), pdf(3,j), pdf(4,j), pdf(5,j)
       enddo
       close(1)
       return
       end

       subroutine init(N,u,v,w,T)
       implicit none
       integer N, j
       double precision u(100000), v(100000), w(100000)
       double precision T(100000), z
 10    format(5g16.8)
       open(100, file="U.dat")
       open(200, file="V.dat")
       open(300, file="W.dat")
       open(400, file="T.dat")
       do j = 1, N
        read(100,10) z, u(j)
        read(200,10) z, v(j)
        read(300,10) z, w(j)
        read(400,10) z, T(j)
       enddo
       close(100)
       close(200)
       close(300)
       close(400)
       return
       end

       subroutine zeroparam(NL,Nt,Np,Na,pdf,time,te,ts)
       implicit none
       integer NL(100000), j
       double precision Nt, Np, Na, pdf(5,1000), time, te, ts
       Nt = 0.d0
       Np = 0.d0
       Na = 0.d0
       time = 0.d0
       te = 0.d0
       ts = 0.d0
       do j=1, 1000
        pdf(1,j) = 0.d0
        pdf(2,j) = 0.d0
        pdf(3,j) = 0.d0
        pdf(4,j) = 0.d0
        pdf(5,j) = 0.d0
       enddo
       do j=1, 100000
        NL(j) = 0
       enddo
       return
       end
       
       subroutine zerovars(N,ua,ur,eu,va,vr,ev,wa,wr,ew,Ta,Tr,eT)
       implicit none
       integer N, j
       double precision ua(100000), ur(100000), eu(100000)
       double precision va(100000), vr(100000), ev(100000)
       double precision wa(100000), wr(100000), ew(100000)
       double precision Ta(100000), Tr(100000), eT(100000)
       do j=1, N
        ua(j) = 0.d0
        ur(j) = 0.d0
        eu(j) = 0.d0
        va(j) = 0.d0
        vr(j) = 0.d0
        ev(j) = 0.d0
        wa(j) = 0.d0
        wr(j) = 0.d0
        ew(j) = 0.d0
        Ta(j) = 0.d0
        Tr(j) = 0.d0
        eT(j) = 0.d0
       enddo
       return
       end
       
       subroutine readpar(idum,N,Lo,Lp,Lm,dt,td,tmax,a)
       implicit none
       integer idum, N, Lo, Lp, Lm
       double precision dt, tmax, td, a(14), ignore
 10    format(3i12)
 20    format(3d15.4)
 30    format(g30.1,g15.4) 
 40    format(g30.1,i15) 
       open(1,file="LabExppar.dat", status="old")
       read(1,10) N
       read(1,20) a(2), a(4), ignore
       read(1,20) a(9), a(10)
       read(1,20) dt, tmax
       read(1,20) a(11), a(12)
       read(1,10) Lo, Lp, Lm
       read(1,10) idum
       close(1)
       a(14) = 2.d0*Lp
       a(1) = dexp(-a(14)/(1.d0*Lm))-dexp(-a(14)/(1.d0*Lo))
       a(1) = a(1)*N/(3.d0*a(14))
       td = 1.d1/(1.d0*N*N)
       if (a(10)*td .gt. 0.1d0) td = 0.1d0/a(10)
       write(6,30) "Dimensionless Parameters:  "
       write(6,30) "   Buoyancy, Temperature:  ", a(2)
       write(6,30) "   Prandtl Number:         ", a(4)
       write(6,30) "   Wind Speed:             ", a(9)
       write(6,30) "   Coriolis Parameter:     ", a(10)
       write(6,30) 
       write(6,30) "Model Parameters:          "
       write(6,30) "   ZC2 =                   ", a(11)
       write(6,30) "   KE mixing ratio =       ", a(12)
       write(6,30) "   a1 =                    ", a(1)
       write(6,30) "   Smallest eddy:          ", (3.d0*Lo)/(1.d0*N)
       write(6,30) "   Most likely eddy:       ", (3.d0*Lp)/(1.d0*N)
       write(6,30) "   Largest eddy:           ", (3.d0*Lm)/(1.d0*N)
       write(6,30) 
       write(6,40) "Grid Points:               ", N
       write(6,30) "Total Simulation Time:     ", tmax
       write(6,30) "Diffusive Time Step:       ", td
       write(6,30) "Initial Eddy Time Step:    ", dt
       return
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
       
       subroutine energy(N,M1,M2,bT,u,v,w,T)
       implicit none
       integer N, M1,M2, j
       double precision bT, Eu, Ev, Ew, ET, EK
       double precision u(100000), v(100000), w(100000), T(050000)
 10    format(i12,3g16.6)
       Eu = 0.d0
       Ev = 0.d0
       Ew = 0.d0
       ET = 0.d0
       do j = M1, M2
        Eu = Eu + (u(j)*u(j))
        Ev = Ev + (v(j)*v(j))
        Ew = Ew + (w(j)*w(j))
        ET = ET + (T(j)*j)
       enddo
       Eu = Eu/(2.d0*N)
       Ev = Ev/(2.d0*N)
       Ew = Ew/(2.d0*N)
       EK = Eu + Ev + Ew
       ET = -ET*27.d0*bT/(8.d0*N*N)
       write(6, 10) M2-M1, EK+ET, EK, ET
       return
       end
