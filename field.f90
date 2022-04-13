program field

double precision, parameter :: pi=3.14159265
double precision a,b,Dx,h,x,t0,tf,I0
integer*8 Npts,N
double precision f

double precision, dimension(:), allocatable :: integr

external f

! Integral parametrica
a = 1.05d0
b = 5.d0
Npts = 1e2
Dx = dble((b-a)/Npts)

allocate(integr(Npts))

! Parametros del integrador
N = 1e4 !Discretizacion del dominio de integracion
h = pi/(N*2.d0)

! Inicializacion
x = a

! Calculo de la integral
open(unit=1,file='data.dat')

do i=1,(Npts+1)
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

   I0 = h*I0
   integr(i) = I0

   write(1,*) x,I0
   
   x = x+Dx
enddo

close(1)


end program field


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
