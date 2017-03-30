!2d electromagnetic field evolution, no сharge and current
program field_evol
    
    implicit none
    !parameters
    integer, parameter :: n=200,nt=200
    double precision, parameter :: dt=1,dx=1,dy=1,c=.4,k1x=0.005,k1y=0
    !variables
    double precision pi,omega,kx,ky
    double precision, dimension(0:n+1,0:n+1) :: Ex,Ey,Bz!,Exo,Eyo,Bzo
    integer i,j,t
    
    pi=acos(-1.)
    kx=2*pi*k1x
    ky=2*pi*k1y
    omega=c*sqrt(kx*kx+ky*ky)
    
    !initial field setting
    t=0
    do i = 0,n+1
        do j = 0,n+1
            Ex(i,j)=0 !Ex^t_{i+1/2,j}
            Ey(i,j)=cos(kx*i+ky*j-omega*t) !Ey^t_{i,j+1/2}
            Bz(i,j)=cos(kx*i+ky*j-omega*(t+.5)) !Bz^{t+1/2}_{i+1/2,j+1/2}
        end do
    end do
    
    !main loop
    do t = 1,nt
        !time t
        write (*,*) t
        do i = 1,n
            do j = 1,n
                Ex(i,j)=Ex(i,j)+c*(Bz(i,j)-Bz(i,j-1))
                Ey(i,j)=Ey(i,j)-c*(Bz(i,j)-Bz(i-1,j))
            end do
        end do
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
        do i = 1,n
            do j = 1,n
                Bz(i,j)=Bz(i,j)-c*(Ex(i,j)+Ey(i+1,j)-Ex(i,j+1)-Ey(i,j))
            end do
        end do
        !set boundaries
        Bz(0,:)=Bz(n,:)
        Bz(n+1,:)=Bz(1,:)
        Bz(:,0)=Bz(:,n)
        Bz(:,n+1)=Bz(:,1)
    end do
    
    !output
    open(2,file="field_evol_Ex.dat")
    open(3,file="field_evol_Ey.dat")
    open(4,file="field_evol_Bz.dat")
    do i=1,n
        write (2,*) Ex(i,1:n)
        write (3,*) Ey(i,1:n)
        write (4,*) Bz(i,1:n)
    end do
    
end program field_evol
