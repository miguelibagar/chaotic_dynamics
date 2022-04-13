program ode

double precision, parameter :: pi=3.14159265
double precision alpha,eps,omega,r,a,ar

integer*8, parameter :: Niter = 1e4
double precision t,tf,h
double precision y0,v0,yold,vold
double precision f,phi,integr1,integr2
double precision k11,k12,k21,k22,k31,k32,k41,k42

double precision :: y(Niter+1),v(Niter+1)
  

external f,phi
  
! Parametros del sistema
alpha = 13.0d0
eps = 30.0d0
omega = 1.d0
r = 0.007d0
a = 0.02d0

ar = a/r

! Parametros del integrador
t = 0.d0
tf = 20.d0

h = (tf-t)/dble(Niter)

! Condiciones iniciales
y0 = (2.d0 + a/r)/2.01d0 ! y0\in [1,a/r]
v0 = 0

! Inicializacion
y(1) = y0
v(1) = v0

yold = y0
vold = v0

i = 2

open(1,file = 'ode_solved_2.txt')
write(1,*) "t   y   v"
write(1,*) t,y0,v0

do while (t .le. tf)
   integr1 = phi(yold,f)
   integr2 = phi(yold-ar,f)
   
   k11 = h*vold
   k12 = -alpha*h*vold - 2.d0*pi*alpha*h*yold + &
        2.d0*alpha*eps*h*cos(omega*t)*(integr1 + integr2)

   integr1 = phi(yold + k11/2.d0,f)
   integr2 = phi(yold-ar + k11/2.d0,f)
   
   k21 = h*vold + h*k12/2.d0
   k22 = -alpha*h*(vold+k12/2.d0) - &
        2.d0*pi*h*alpha*(yold+k11/2.d0) + &
        2.d0*alpha*h*eps*cos(omega*(t+h/2.d0))*(integr1 + integr2)

   integr1 = phi(yold + k21/2.d0,f)
   integr2 = phi(yold-ar + k21/2.d0,f)
   
   k31 = h*vold + h*k22/2.d0
   k32 = -alpha*h*(vold+k22/2.d0) - &
        2.d0*pi*h*alpha*(yold+k21/2.d0) + &
        2.d0*alpha*h*eps*cos(omega*(t+h/2.d0))*(integr1 + integr2)

   integr1 = phi(yold + k31,f)
   integr2 = phi(yold-ar + k31,f)
   
   k41 = h*vold + h*k32
   k42 = -alpha*h*(vold+k32) - &
        2.d0*pi*alpha*h*(yold+k31) + &
        2.d0*alpha*h*eps*cos(omega*(t+h))*(integr1 + integr2)

   yold = yold + k11/6.d0 + k21/3.d0 + k31/3.d0 + k41/6.d0
   vold = vold + k12/6.d0 + k22/3.d0 + k32/3.d0 + k42/6.d0

   if (yold.le.1.or.yold.ge.(1+ar)) exit
   
   y(i) = yold
   v(i) = vold

   write(*,*) "Iteracion ", i
   
   t = t + h
   i = i+1

   write(1,*) t,yold,vold
enddo

end program ode


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


double precision function phi(x,f)

double precision, parameter :: pi=3.14159265
double precision h,x,t0,tf,I0
integer*8 N
double precision f

external f


! Parametros del integrador
N = 1e4 !Discretizacion del dominio de integracion
h = pi/(N*2.d0)

! Calculo de la integral
I0 = 0.d0

t0 = 0.d0
tf = pi/2
   
t0 = t0 + h
tf = tf - h
I0 = I0 + dble(55.d0/24.d0)*(f(x,t0) + f(x,tf))

t0 = t0 + h
tf = tf - h
I0 = I0 - dble(1.d0/6.d0)*(f(x,t0) + f(x,tf))

t0 = t0 + h
tf = tf - h
I0 = I0 + dble(11.d0/8.d0)*(f(x,t0) + f(x,tf))

do j=1,(N-7)
   t0 = t0 + h
   I0 = I0 + f(x,t0)
enddo

phi = h*I0
return
end function phi


!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function f(x,t)
double precision x,t
double precision num,den,ct

ct = cos(t)
num = x-ct
den = (1+x*x - 2.d0*x*ct)**(1.5d0)

f = num/den
return
end function f
