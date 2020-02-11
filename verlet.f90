program verlet
implicit none
integer,parameter::n=9
real,dimension(n)::x,y,xp,yp,fx,fy,vx,vy,m
real::dt,t,D,g,ms, rcmx, rcmy
integer::Nend,i, j

Nend=10000
dt=0.001
t=0
rcmx=0.0
rcmy=0.0

D=1.0

g=1.15e-4

m(9)=333054.25  !msol (le pongo como un vector para el centro de masa)
m(1)=0.055
m(2)=0.815
m(3)=1.0
m(4)=0.107
m(5)=317.8
m(6)=95.16
m(7)=14.54
m(8)=17.15

y(1)=0.39*D
y(2)=0.72*D
y(3)=D
y(4)=1.52*D
y(5)=5.2*D
y(6)=9.54*D
y(7)=19.19*D
y(8)=30.06*D
y(9)=0.0

x(1)=0.0
x(2)=0.0
x(3)=0.0
x(4)=0.0
x(5)=0.0
x(6)=0.0
x(7)=0.0
x(8)=0.0
x(9)=0.0


vx(1)=10.09
vx(2)=-7.3827
vx(3)=6.0682
vx(4)=5.087
vx(5)=2.7552
vx(6)=2.039
vx(7)=-1.4409
vx(8)=1.1547
vx(9)=0.0

vy(1)=0.0
vy(2)=0.0
vy(3)=0.0
vy(4)=0.0
vy(5)=0.0
vy(6)=0.0
vy(7)=0.0
vy(8)=0.0
vy(9)=0.0

fx=0.0
fy=0.0

open(file='datos.dat',unit=15,status='unknown')
open(file='datoscm.dat',unit=16,status='unknown')

  write(15,*) x(1),y(1),x(2),y(2),x(3),y(3),x(4),y(4),x(5),y(5),x(6),y(6),x(7),y(7),x(8),y(8),x(9),y(9)

do i=1,Nend
  t=i*dt
  x(:)=x(:)+vx(:)*dt+(0.5*fx(:)*dt*dt)
  y(:)=y(:)+vy(:)*dt+(0.5*fy(:)*dt*dt)
  write(15,*) x(1),y(1),x(2),y(2),x(3),y(3),x(4),y(4),x(5),y(5),x(6),y(6),x(7),y(7),x(8),y(8),x(9),y(9)
  vx(:)=vx(:)+(0.5*fx(:)*dt)
  vy(:)=vy(:)+(0.5*fy(:)*dt)
  call crear_fuerza(n,m,x,y,fx,fy)
  vx(:)=vx(:)+(0.5*fx(:)*dt)
  vy(:)=vy(:)+(0.5*fy(:)*dt)
  call cm(n, m, x, y, rcmx, rcmy)
  do j=1,n
    xp(j)=x(j)-rcmx
    yp(j)=y(j)-rcmy
  enddo
  write(16,*) xp(1),yp(1),xp(2),yp(2),xp(3),yp(3),xp(4),yp(4),xp(5),yp(5),xp(6),yp(6),xp(7),yp(7),xp(8),yp(8),xp(9),yp(9)
enddo

close(15)
close(16)
contains
subroutine crear_fuerza(n,m,x,y,fx,fy)
implicit none
!recibe el vector de posiciones de las particulas y entrega el vector de fuerzas
integer::n, i, j
real::m(n),x(n),y(n),fx(n),fy(n)
!auxiliares
!real::dis(n-1)

fx=0.0
fy=0.0
               
do i=1, n
  !dis(i)=sqrt(x(i)*x(i)+y(i)*y(i))
  !fx(i)=-g*m(9)*x(i)/(dis(i)**(3))
  !fy(i)=-g*m(9)*y(i)/(dis(i)**(3))

  do j=1, n
  
  if(i .ne. j) then

  fx(i)=fx(i)-(g*m(j)*(x(i)-x(j))/((x(i)-x(j))*(x(i)-x(j))+(y(i)-y(j))*(y(i)-y(j)))**(1.5))
  fy(i)=fy(i)-(g*m(j)*(y(i)-y(j))/((x(i)-x(j))*(x(i)-x(j))+(y(i)-y(j))*(y(i)-y(j)))**(1.5))

  endif

  enddo
enddo

end subroutine crear_fuerza


subroutine cm(n, m, x, y, rcmx, rcmy)
implicit none
integer::i, n
real::m(n), x(n), y(n)
real::mcm, sumax, sumay, rcmx, rcmy 

mcm=0.0
sumax=0.0
sumay=0.0

do i=1, n
 mcm=mcm+m(i)
enddo


do i=1, n
  sumax=sumax+m(i)*x(i)
  sumay=sumay+m(i)*y(i)
enddo

rcmx=sumax/mcm
rcmy=sumay/mcm


end subroutine cm


end program
