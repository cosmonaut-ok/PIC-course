!Simple particle moving program using external E and B distribution
!implies relativistic Boris pusher

program pmove

    implicit none
    !parameters
    integer, parameter :: n=20000
    double precision, parameter :: dt=1e-3, q=1, m=1, c=1
    double precision, parameter :: r0(3)=(/0.,0.,0./), u0(3)=(/0.5,0.,0./) !u will be real velocity
    !variables
    double precision, dimension(n,3) :: r,v,u
    double precision, dimension(3) :: vm,vp,vs,s,t,gf
    integer i
    
    !set initial conditions
    r(1,:)=r0
    u(1,:)=u0
    v(1,:)=u0/sqrt(1-scalp(u0,u0)) !v(i) is going to be v_{n-1/2}
    !open file
    open(2,file="pmove.dat")
    !main loop
    do i = 1,n-1
        vm=v(i,:)+(q/m*dt/2)*E(r(i,:))
        gf=sqrt(1+scalp(vm,vm))
        t=(q*dt/(2*m*c*gf))*B(r(i,:))
        vs=vm+vecp(vm, t)
        s=2/(1+scalp(t,t))*t
        vp=vm+vecp(vs, s)
        v(i+1,:)=vp+(q/m*dt/2)*E(r(i,:))
        u(i+1,:)=v(i+1,:)/(1+scalp(v(i+1,:),v(i+1,:)))
        r(i+1,:)=r(i,:)+dt*u(i+1,:)
        write (2,*) r(i+1,:)
        write (*,*) r(i+1,:)
    end do
    
    contains
    
    function E(r)
        double precision E(3)
        double precision r(3)
        E=(/0.,0.5,0./)
    end function E
    
    function B(r)
        double precision B(3)
        double precision r(3)
        B=(/0,0,1/)
        !B(3)=1+r(1)
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
    
end program pmove
