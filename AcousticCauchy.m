clear all, clc
close all
%Constants
r0=1;
c0=1;

%Eigenvectors and eigenvalues
P=(1/(sqrt(1+(r0*c0)^2)))*[1 1; -r0*c0 r0*c0];
L=[-c0 0; 0 c0];

%Problem dimension
a=-3; b=3; t0=0; tf=1;
dimx=60;
deltax=(b-a)/(dimx+1);
cfl=0.4;
dimt=max(max(abs(L)))*(tf-t0)/(cfl*deltax)-1;
deltat=(tf-t0)/(dimt+1);

%Initial conditions
uu=zeros(1,dimx+2);
pp=zeros(1,dimx+2);
x=a:deltax:b;
for i=1:dimx+2
uu(i)=u0fun(x(i));
pp(i)=p0fun(x(i));
end
%
%Initialization
u=uu;
p=pp;
w1=uu;
w2=pp;
%
%Plot the initial condition   
movegui(figure(1),[120 120])
hold off
plot(x,uu,'-black')
xlabel('x')
title('Initial condition u0(x)');
ylabel('u_0(x)')
xlim([a b]);
%
movegui(figure(2),[684 120])
plot(x,pp,'-black')
xlabel('x')
title('Initial condition p0(x)');
ylabel('p_0(x)')
xlim([a b]);
pause
%Solver
for n=1:dimt+1
%Exact solution w1,w2
  for i=1:dimx+2
    w1(i) = (1/(2*r0*c0))*(r0*c0*(u0fun(x(i)+c0*n*deltat) +u0fun(x(i)-c0*n*deltat)) + p0fun(x(i)-c0*n*deltat) -p0fun(x(i)+c0*n*deltat) );   
    w2(i) = (1/2)*(r0*c0*(u0fun(x(i)-c0*n*deltat) -u0fun(x(i)+c0*n*deltat)) + p0fun(x(i)-c0*n*deltat) +p0fun(x(i)+c0*n*deltat) );
  end 
%Godunov solver
  for i=2:dimx+1
    u(i)=uu(i)-(deltat/(2*deltax))*((1/r0)*(pp(i+1)-pp(i-1)) - abs(c0)*(uu(i+1)-2*uu(i)+uu(i-1)));
    p(i)=pp(i)-(deltat/(2*deltax))*((r0*c0^2)*(uu(i+1)-uu(i-1)) - abs(c0)*(pp(i+1)-2*pp(i)+pp(i-1)));
  end
%Boundary conditions
  u(1)=w1(1);
  u(dimx+2)=w1(dimx+2);
  p(1)=w2(1);
  p(dimx+2)=w2(dimx+2);
%Redefinition        
  uu=u; pp=p;
%Plot the solution   
  figure(3)
  set(figure(3),'units','normalized','outerposition',[0 0 1 1])
  subplot(2,2,1)
  plot(x,w1,'-black',x,uu,'.-blue')
  title(['t = ',num2str(n*deltat)]);
  xlabel('x')
  ylabel('u(x,t)')
  xlim([a b]);  
  ylim([min(min(u),min(w1)) max(max(u),max(w1))]);
%   
  subplot(2,2,3)
  plot(x,w2,'-black',x,pp,'.-red')
  title(['t = ',num2str(n*deltat)]);
  xlabel('x')
  ylabel('p(x,t)')
  xlim([a b]);
  ylim([min(min(p),min(w2)) max(max(p),max(w2))]);
%Plot error
  subplot(2,2,2)
  plot(x,abs(w1-u),'.-blue')
  xlabel('x'); ylabel('Absolute error u(x,t)');
  xlim([a b])
  title(['Absolute error u(x,t), t = ',num2str(n*deltat)]);
%
  subplot(2,2,4)
  plot(x,abs(w2-p),'.-red')
  xlabel('x'); ylabel('Absolute error p(x,t)');
  xlim([a b])
  title(['Absolute error p(x,t), t = ',num2str(n*deltat)]);
drawnow
end

