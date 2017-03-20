!Simple particle moving program using external E and B distribution
!implies relativistci Boris pusher

program pmove

    implicit none
    !parameters
    integer, parameter :: n=20000
    double precision, parameter :: dt=1e-3, q=1, m=1, c=1
    double precision, parameter :: r0(3)=(/1.,1.,0./), u0(3)=(/0.05,0.,0./) !u will be real velocity
    !variables
    double precision, dimension(n,3) :: r,v,u
    double precision, dimension(3) :: vs,t,gf,v2,tau,Es,Bs
    double precision :: vss,sigma,s
    integer i
    
    !set initial conditions
    r(1,:)=r0 !r(i) is going to be r(i+1/2)
    u(1,:)=u0 !v in Vay's terms
    v(1,:)=u0/sqrt(1-scalp(u0,u0)) !u in Vay's terms
    !open file
    open(2,file="NS-sim.dat")
    !main loop
    do i = 1,n-1
        Es=E(r(i,:))
        Bs=B(r(i,:))
        v2=v(i,:)+q*dt/(2*m)*(Es+vecp(u(i,:),Bs)) !u_{i+1/2}
        vs=v2+q*dt/(2*m)*Es !u'
        tau=q*dt/(2*m)*Bs
        vss=scalp(vs,tau) !u*
        sigma=1+scalp(vs,vs)-scalp(tau,tau)
        gf=sqrt((sigma+sqrt(sigma*sigma+4*(scalp(tau,tau)+vss*vss)))/2) !gamma_{i+1}
        t=tau/gf
        s=1/(1+scalp(t,t))
        v(i+1,:)=s*(vs+scalp(vs,t)*t+vecp(vs,t)) !u_{i+1}
        u(i+1,:)=v(i+1,:)/gf !v_{i+1}
        r(i+1,:)=r(i,:)+dt*u(i+1,:)
        write (2,*) r(i+1,:)
        write (*,*) r(i+1,:)
    end do
    
    contains
    
    function E(r)
        double precision E(3)
        double precision r(3)
        double precision, parameter :: B0=1e-2,R0=1,Omega=1
        double precision rr,theta,phi,etheta(3)
        rr=norm(r)
        theta=acos(r(3)/rr)
        phi=acos(r(1)/sqrt(sum(r(:2)*r(:2))))
        if (r(1)<0) phi=-phi
        etheta(1)=-cos(theta)*cos(phi)
        etheta(2)=-cos(theta)*sin(phi)
        etheta(3)=sin(theta)
        E=-(B0*R0**2*Omega/rr*sin(theta))*etheta
    end function E
    
    function B(r)
        double precision B(3)
        double precision r(3)
        double precision, parameter :: B0=1e-2,R0=1,Omega=1,sigma=1e5
        double precision rr,theta,phi,er(3),ephi(3)
        rr=norm(r)
        theta=acos(r(3)/rr)
        phi=acos(r(1)/sqrt(sum(r(:2)*r(:2))))
        if (r(1)<0) phi=-phi
        er=r/rr
        ephi(1)=-sin(phi)
        ephi(2)=cos(phi)
        ephi(3)=0
        B=(B0*(R0/rr)**2)*er-(B0*R0**2*Omega/rr*sin(theta)*(1+1/sqrt(sigma)))*ephi
    end function B
    
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
    
    function norm(a)
        double precision norm
        double precision a(3)
        norm=sqrt(sum(a*a))
    end function norm
    
end program pmove
