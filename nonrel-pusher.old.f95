!Simple particle moving program using external E and B distribution
!implies non-relativistic Boris pusher

program pmove
type vec
    sequence
    double precision x,y,z
end type vec

    !parameters
    integer, parameter :: n=20000
    double precision, parameter :: dt=1e-3, q=1, m=1, c=1
    type(vec), parameter :: r0=vec(0,0,0), v0=vec(1,0,0)
    !functions
    type(vec) vecp,muln,addv,B,E
    double precision scalp
    !variables
    type(vec) r(n),v(n)
    type(vec) vm,vp,vs,s,t
    
    !set initial conditions
    r(1)=r0
    v(1)=v0 !v(i) is going to be v_{n-1/2}
    !open file
    open(2,file="pmove.dat")
    !main loop
    do i = 1,n-1
        vm=addv(v(i), muln((q/m*dt/2), E(r(i))))
        t=muln((q*dt/(2*m*c)),B(r(i)))
        vs=addv(vm, vecp(vm, t))
        s=muln(2/(1+scalp(t,t)), t)
        vp=addv(vm, vecp(vs, s))
        v(i+1)=addv(vp, muln((q/m*dt/2), E(r(i))))
        r(i+1)=addv(r(i), muln(dt, v(i+1)))
        write (2,*) r(i+1)
        write (*,*) r(i+1), sqrt(scalp(v(i+1),v(i+1)))
    end do
    
end program pmove

function E(r)
type vec
    sequence
    double precision x,y,z
end type vec

    type(vec) E
    type(vec) r
    E=vec(0,0,0)

end function E

function B(r)
type vec
    sequence
    double precision x,y,z
end type vec

    type(vec) B
    type(vec) r
    B=vec(0,0,1+r%x)

end function B

function addv(a,b)
type vec
    sequence
    double precision x,y,z
end type vec

    type(vec) addv
    type(vec) a,b
    addv=vec(a%x+b%x,a%y+b%y,a%z+b%z)

end function addv

function muln(n,a)
type vec
    sequence
    double precision x,y,z
end type vec

    type(vec) muln
    double precision n
    type(vec) a
    muln=vec(n*a%x,n*a%y,n*a%z)

end function muln

function vecp(a,b)
type vec
    sequence
    double precision x,y,z
end type vec

    type(vec) vecp
    type(vec) a,b
    vecp=vec(a%y*b%z-a%z*b%y, -a%x*b%z+a%z*b%x, a%x*b%y-a%y*b%x)
    
end function vecp

function scalp(a,b)
type vec
    sequence
    double precision x,y,z
end type vec

    double precision scalp
    type(vec) a,b
    scalp=a%x*b%x+a%y*b%y+a%z*b%z

end function scalp
