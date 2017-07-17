!2d electromagnetic field evolution, no сharge and current
program field_evol_2d
    
    implicit none
    !parameters
    integer, parameter :: n=200,nt=200
    double precision, parameter :: dt=1d0,dx=1d0,dy=1d0,c=.4d0,k1x=5d-3,k1y=1d-2
    !variables
    double precision pi,omega,k,kx,ky
    double precision, dimension(0:n+1,0:n+1) :: Ex,Ey,Bz!,Exo,Eyo,Bzo
    integer i,j,t
    
    pi=acos(-1d0)
    kx=2*pi*k1x
    ky=2*pi*k1y
    k=sqrt(kx*kx+ky*ky)
    omega=c*k
    
    !initial field setting
    t=0
    do i = 0,n+1
        do j = 0,n+1
            Ex(i,j)=-ky/k*cos(kx*i+ky*j-omega*t) !Ex^t_{i+1/2,j}
            Ey(i,j)=kx/k*cos(kx*i+ky*j-omega*t) !Ey^t_{i,j+1/2}
            Bz(i,j)=cos(kx*(i+.5d0)+ky*j-omega*(t+.5d0)) !Bz^{t+1/2}_{i+1/2,j+1/2}
        end do
    end do
    
    !main loop
    do t = 1,nt
        !time t
        write (*,*) t
        Ex(1:n,1:n)=Ex(1:n,1:n)+c*(Bz(1:n,1:n)-Bz(1:n,0:n-1))
        Ey(1:n,1:n)=Ey(1:n,1:n)-c*(Bz(1:n,1:n)-Bz(0:n-1,1:n))
        !set boundaries
        Ex(0,:)=Ex(n,:)
        Ex(n+1,:)=Ex(1,:)
        Ex(:,0)=Ex(:,n)
        Ex(:,n+1)=Ex(:,1)
        Ey(0,:)=Ey(n,:)
        Ey(n+1,:)=Ey(1,:)
        Ey(:,0)=Ey(:,n)
        Ey(:,n+1)=Ey(:,1)
        !time t+1/2
        write (*,*) t,"+1/2"
        Bz(1:n,1:n)=Bz(1:n,1:n)-c*(Ex(1:n,1:n)+Ey(2:n+1,1:n)-Ex(1:n,2:n+1)-Ey(1:n,1:n))
        !set boundaries
        Bz(0,:)=Bz(n,:)
        Bz(n+1,:)=Bz(1,:)
        Bz(:,0)=Bz(:,n)
        Bz(:,n+1)=Bz(:,1)
    end do
    
    !output
    open(2,file="field-evol-Ex.dat")
    open(3,file="field-evol-Ey.dat")
    open(4,file="field-evol-Bz.dat")
    do i=1,n
        write (2,*) Ex(i,1:n)
        write (3,*) Ey(i,1:n)
        write (4,*) Bz(i,1:n)
    end do
    
end program
