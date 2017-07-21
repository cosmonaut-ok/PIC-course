!Particle in cell code with Yee mesh for fields, Boris particle pusher and periodic BC
!not finished yet, thing to do is generate initial particles/fields
!not relativistic stream test
program pic_2d
    
    implicit none
    
    ! data structure for particles
    type part
        double precision r(3), v(3), u(3), m, q ! u is for real velocity
    end type part
    
    ! parameters
    integer, parameter :: n=100, nt=1000, ppc=10, np=n*n*ppc !grid size, number of timesteps, particles per cell and number of particles
    double precision, parameter :: L=5d0, sigma=1d0, qm=1d0 !L=Debay radius/dx, sigma is magnetisation parameter (not used yet), qm is charge-mass ratio
    double precision, parameter :: dt=1d0, dx=1d0, dy=1d0, c=.4d0
    double precision, parameter :: vth=.2d0*c, vst=.1d0*c !thermal and stream velocities
    double precision qe, me
    
    ! variables
    integer t, i, j, k, i1, j1, i2, j2
    double precision pi, wx, wy, gf, x1, y1, x2, y2, xr, yr, wx1, wx2, wy1, wy2, fx1, fx2, fy1, fy2, rn(4)
    double precision, dimension(3) :: E,B, vm, vp, vs, tt, s
    double precision, dimension(0:n+1,0:n+1) :: Ex,Ey,Bz, Jx, Jy
    type(part) p(2*np)
    
    !Important note
    !Ex_{i+1/2,j} for Ex(i,j)
    !Ey_{i,j+1/2} for Ey(i,j)
    !Bz_{i+1/2,j+1/2} for Bz(i,j), initially retarded by 1/2 in time
    
    pi=acos(-1d0)
    qe=c*c/(L*L*4*pi*ppc)
    me=qe/qm
    call random_seed()
    ! generate particles
    ! positrons
    do k = 1, np/2
        p(k)%q=qe
        p(k)%m=me
        ! generate positions and velocities
        call random_number(rn)
        p(k)%r(1:2)=1+n*rn(1:2)
        p(k)%r(3)=0
        p(k)%v(1)=vst+vth*sqrt(-log(rn(3)))*cos(2*pi*rn(4))
        p(k)%v(2)=vth*sqrt(-log(rn(3)))*sin(2*pi*rn(4))
        p(k)%v(3)=0
        p(k)%u=p(k)%v/sqrt(1+scalp(p(k)%v,p(k)%v)/c/c)
    end do
    ! electrons
    do k = np/2+1, np
        p(k)%q=-qe
        p(k)%m=me
        ! generate positions and velocities
        call random_number(rn)
        p(k)%r(1:2)=1+n*rn(1:2)
        p(k)%r(3)=0
        p(k)%v(1)=-vst+vth*sqrt(-log(rn(3)))*cos(2*pi*rn(4))
        p(k)%v(2)=vth*sqrt(-log(rn(3)))*sin(2*pi*rn(4))
        p(k)%v(3)=0
    end do
    ! generate fields
    !now just zeros
    Ex=0
    Ey=0
    Bz=0
    
    ! main loop
    do t = 1, nt
        !half-move B
        Bz(1:n,1:n)=Bz(1:n,1:n)-c*dt/2*(Ex(1:n,1:n)*dx+Ey(2:n+1,1:n)*dy-Ex(1:n,2:n+1)*dx-Ey(1:n,1:n)*dy)/dx/dy
        !set boundaries
        Bz(0,:)=Bz(n,:)
        Bz(n+1,:)=Bz(1,:)
        Bz(:,0)=Bz(:,n)
        Bz(:,n+1)=Bz(:,1)
        !clear currents
        Jx=0
        Jy=0
        do k = 1, np
            !interpolate fields to particle
            E(3)=0
            B(1)=0
            B(2)=0
            !Ex
            i=int(p(k)%r(1)/dx-.5d0)
            j=int(p(k)%r(2)/dy)
            wx=p(k)%r(1)/dx-.5d0-i
            wy=p(k)%r(2)/dy-j
            E(1)=Ex(i,j)*(1-wx)*(1-wy)+Ex(i,j+1)*(1-wx)*wy+Ex(i+1,j)*wx*(1-wy)+Ex(i+1,j+1)*wx*wy
            !Bz, changed only j and wy
            j=int(p(k)%r(2)/dy-.5d0)
            wy=p(k)%r(2)/dy-.5d0-j
            B(3)=Bz(i,j)*(1-wx)*(1-wy)+Bz(i,j+1)*(1-wx)*wy+Bz(i+1,j)*wx*(1-wy)+Bz(i+1,j+1)*wx*wy
            !Ey, changed only i and wx
            i=int(p(k)%r(1)/dx)
            wx=p(k)%r(1)/dx-i
            E(2)=Ey(i,j)*(1-wx)*(1-wy)+Ey(i,j+1)*(1-wx)*wy+Ey(i+1,j)*wx*(1-wy)+Ey(i+1,j+1)*wx*wy
            !particle push
            vm=p(k)%v+(p(k)%q/p(k)%m*dt/2)*E
            gf=sqrt(1+scalp(vm,vm)/c/c)
            tt=(p(k)%q*dt/(2*p(k)%m*c*gf))*B
            vs=vm+vecp(vm, tt)
            s=2/(1+scalp(tt,tt))*tt
            vp=vm+vecp(vs, s)
            p(k)%v=vp+(p(k)%q/p(k)%m*dt/2)*E
            p(k)%u=p(k)%v/(1+scalp(p(k)%v,p(k)%v)/c/c)
            x1=p(k)%r(1)
            y1=p(k)%r(2)
            p(k)%r=p(k)%r+p(k)%u*dt
            x2=p(k)%r(1)
            y2=p(k)%r(2)
            !interpolate J
            !calculate relay point
            i1=int(x1/dx)
            i2=int(x2/dx)
            j1=int(y1/dy)
            j2=int(y2/dy)
            xr=min(min(i1,i2)*dx+dx,max(max(i1,i2)*dx,(x1+x2)/2))
            yr=min(min(j1,j2)*dy+dy,max(max(j1,j2)*dy,(y1+y2)/2))
            !calculate different things
            Fx1=p(k)%q/dt*(xr-x1)
            Fy1=p(k)%q/dt*(yr-y1)
            Fx2=p(k)%q/dt*(x2-xr)
            Fy2=p(k)%q/dt*(y2-yr)
            Wx1=(x1+xr)/2/dx-i1
            Wy1=(y1+yr)/2/dy-j1
            Wx2=(x2+xr)/2/dx-i2
            Wy2=(y2+yr)/2/dy-j2
            !finally interpolate current
            Jx(i1,j1)=Jx(i1,j1)+Fx1*(1-Wy1)/dx/dy
            Jx(i2,j2)=Jx(i2,j2)+Fx2*(1-Wy2)/dx/dy
            Jy(i1,j1)=Jy(i1,j1)+Fy1*(1-Wx1)/dx/dy
            Jy(i2,j2)=Jy(i2,j2)+Fy2*(1-Wx2)/dx/dy
            Jx(i1,j1+1)=Jx(i1,j1+1)+Fx1*Wy1/dx/dy
            Jy(i1+1,j1)=Jy(i1+1,j1)+Fy1*Wx1/dx/dy
            Jx(i2,j2+1)=Jx(i2,j2+1)+Fx2*Wy2/dx/dy
            Jy(i2+1,j2)=Jy(i2+1,j2)+Fy2*Wx2/dx/dy
            !move particles crossing outer boundary
            p(k)%r=modulo(p(k)%r-1,float(n))+1
        end do
        !add boundaries
        Jx(n,:)=Jx(n,:)+Jx(0,:)
        Jx(1,:)=Jx(1,:)+Jx(n+1,:)
        Jx(:,n)=Jx(:,n)+Jx(:,0)
        Jx(:,1)=Jx(:,1)+Jx(:,n+1)
        Jy(n,:)=Jy(n,:)+Jy(0,:)
        Jy(1,:)=Jy(1,:)+Jy(n+1,:)
        Jy(:,n)=Jy(:,n)+Jy(:,0)
        Jy(:,1)=Jy(:,1)+Jy(:,n+1)
        !one more half-move B
        Bz(1:n,1:n)=Bz(1:n,1:n)-c*dt/2*(Ex(1:n,1:n)*dx+Ey(2:n+1,1:n)*dy-Ex(1:n,2:n+1)*dx-Ey(1:n,1:n)*dy)/dx/dy
        !set boundaries
        Bz(0,:)=Bz(n,:)
        Bz(n+1,:)=Bz(1,:)
        Bz(:,0)=Bz(:,n)
        Bz(:,n+1)=Bz(:,1)
        !move E
        Ex(1:n,1:n)=Ex(1:n,1:n)+c*dt*(Bz(1:n,1:n)-Bz(1:n,0:n-1))-4*pi*Jx(1:n,1:n)*dt
        Ey(1:n,1:n)=Ey(1:n,1:n)-c*dt*(Bz(1:n,1:n)-Bz(0:n-1,1:n))-4*pi*Jy(1:n,1:n)*dt
        !set boundaries
        Ex(0,:)=Ex(n,:)
        Ex(n+1,:)=Ex(1,:)
        Ex(:,0)=Ex(:,n)
        Ex(:,n+1)=Ex(:,1)
        Ey(0,:)=Ey(n,:)
        Ey(n+1,:)=Ey(1,:)
        Ey(:,0)=Ey(:,n)
        Ey(:,n+1)=Ey(:,1)
    end do
    
    !some outputs
    !positrons
    open(2,file="stream-positrons.dat")
    do k = 1, np/2
        write (2,*) p(k)%r(1:2), p(k)%u(1:2)
    end do
    !electrons
    open(3,file="stream-electrons.dat")
    do k = np/2+1, np
        write (3,*) p(k)%r(1:2), p(k)%u(1:2)
    end do
    
    contains
    
    function vecp(a,b)
        double precision, dimension(3) :: vecp
        double precision, dimension(3) :: a,b
        vecp(1)=a(2)*b(3)-a(3)*b(2)
        vecp(2)=-a(1)*b(3)+a(3)*b(1)
        vecp(3)=a(1)*b(2)-a(2)*b(1)
    end function vecp
    
    function scalp(a,b)
        double precision scalp
        double precision a(3),b(3)
        scalp=sum(a*b)
    end function scalp
    
!    function sf(p1,p2) ! shape function
!        double precision p1,p2
!        double precision sf
!        sf=1d0-min(abs(p1-p2),1d0)
!    end function sf
    
end program
