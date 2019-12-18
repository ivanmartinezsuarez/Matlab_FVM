clear all, clc
close all
%Constants
r0=1;
c0=1;

%Initial conditions
ul=1; pl=5;
ur=0; pr=10;
wl=[ul; pl];
wr=[ur; pr];

%Eigenvectors and eigenvalues
P=(1/(sqrt(1+(r0*c0)^2)))*[1 1; -r0*c0 r0*c0];
L=[-c0 0; 0 c0];
alphal=P\wl;
alphar=P\wr;
wm=wl+(alphar(1)-alphal(1))*P(:,1);

%Problem dimension
a=-3; b=3; t0=0; tf=2;
dimx=200;
%dimt=60;
deltax=(b-a)/(dimx+1);
%deltat=(tf-t0)/(dimt+1);
%cfl=max(max(abs(L)))*deltat/deltax;
cfl=0.8;
dimt=max(max(abs(L)))*(tf-t0)/(cfl*deltax)-1;
deltat=(tf-t0)/(dimt+1);

%Initial conditions
x=a:deltax:b;
%
uu(1:dimx+2)=ur;
uu(1:(dimx+2)/2)=ul;
pp(1:dimx+2)=pr;
pp(1:(dimx+2)/2)=pl;
w1=uu;
w2=pp;
%
u=zeros(1,dimx+2);
p=zeros(1,dimx+2);

%Plot the initial condition  
movegui(figure(1),[120 120])
hold off
plot(x,uu,'.-black')
xlabel('x')
title('Initial condition u0(x)');
ylabel('u_0(x)')
xlim([a b]);
%
movegui(figure(2),[684 120])
hold off
plot(x,pp,'.-black')
xlabel('x')
title('Initial condition p0(x)');
ylabel('p_0(x)')
xlim([a b]);
pause
for n=1:dimt+1
  %Exact solution
  for i=1:dimx+2
    if (x(i)-L(1,1)*n*deltat)*(x(i)-L(2,2)*n*deltat)<0
      w1(i)=wm(1);
      w2(i)=wm(2);
    elseif x(i)-L(1,1)*n*deltat<0
      w1(i)=wl(1);
      w2(i)=wl(2);
    elseif x(i)-L(2,2)*n*deltat>0
      w1(i)=wr(1);
      w2(i)=wr(2);
    end
  end
  % 
  %Godunov solution
  for i=2:dimx+1
    u(i)=uu(i)-(deltat/(2*deltax))*((1/r0)*(pp(i+1)-pp(i-1)) - c0*(uu(i+1)-2*uu(i)+uu(i-1)));
    p(i)=pp(i)-(deltat/(2*deltax))*((r0*c0^2)*(uu(i+1)-uu(i-1)) - c0*(pp(i+1)-2*pp(i)+pp(i-1)));
  %Transmissive conditions
    u(1)=u(2);
    u(dimx+2)=u(dimx+1);
    p(1)=p(2);
    p(dimx+2)=p(dimx+1);
  end
  %Redefinition
    uu=u;
    pp=p;   
  %Plot the solution   
  movegui(figure(3),[120 120])
  hold off
  plot(x,w1,'-black',x,uu,'.-blue')
  xlabel('x')
  title(['t = ',num2str(n*deltat)]);
  ylabel('u(x,t)')
  %
  movegui(figure(4),[684 120])
  hold off
  plot(x,w2,'-black',x,pp,'.-red')
  xlabel('x')
  title(['t = ',num2str(n*deltat)]);
  ylabel('p(x,t)')
  %
%   %Error plot
%   figure(5)
%   plot(x,abs(w1-u),'-+',x,abs(w2-p),'-*')
%   xlabel('x'); ylabel('Error');
%   xlim([a b])
%   title(['Acoustic equations, t = ',num2str(n*deltat)]);
%   legend('Velocity error','Pressure error')
%   %
  drawnow
end

