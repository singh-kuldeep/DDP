clear all
close all
clc
clear
% One gas plots 
Array=csvread('premetive.csv');
extra=csvread('extra.csv');
x = Array(:, 1);
rho = Array(:, 2);
u = Array(:, 3);
p = Array(:, 4);
T = p./(287.15*rho);
[m,n] = size(Array);

dx = extra(1,1);
t = extra(1,2);

% p2sim = 
for i = 2:m-1
	d2u(i) = (u(i+1)-(2*u(i))+u(i-1))/(dx*dx);
end

[Max,I] = max(d2u);

rho21Ratiosim = rho(I-5)/rho(I+5) ; 
u21Ratiosim = u(I-5)/u(I+5) ;
p21Ratiosim = p(I-5)/p(I+5) 	
T21Ratiosim = (p(I-5)*rho(I+5))/(rho(I-5)*p(I+5)) ;

% Theory
rho4 = extra(1,3) ; 
u4 = extra(1,4) ;
p4 = extra(1,5) ;

rho1 = extra(1,6) ; 
u1 = extra(1,7) ;
p1 = extra(1,8) ;

T41Ratiotheo = p4*rho1/(p1*rho4);

fun = @(p21Ratiotheo) (p4/p1) - p21Ratiotheo*(1-(0.4*(sqrt(1/T41Ratiotheo))*(p21Ratiotheo-1))/sqrt(2.8*(2.8+(2.4)*(p21Ratiotheo-1))))^(-7.0); % function
p21Ratiotheo0 = [0 10]; % initial interval

% options = optimset('Display','iter'); % show iterations
% [p21Ratiotheo fval exitflag output] = fzero(fun,p21Ratiotheo0,options);

[p21Ratiotheo fval exitflag output] = fzero(fun,p21Ratiotheo0);

T21Ratiotheo = p21Ratiotheo*((6.0+p21Ratiotheo)/(1.0+6.0*p21Ratiotheo));
rho21Ratiotheo = p21Ratiotheo/T21Ratiotheo;

% simulation Shock Spped
ussim = ((I-floor(m/2))*dx)/t
% theoretical Shock speed
usTheo = sqrt( (1.4*p1/rho1)*( 1 + ((1.4+1)/(2*1.4))*(p21Ratiotheo -1) ) )	

p21RatioErr = (p21Ratiotheo- p21Ratiosim)*100/p21Ratiotheo
T21RatioErr = (T21Ratiotheo- T21Ratiosim)*100/T21Ratiotheo
rho21RatioErr = (rho21Ratiotheo- rho21Ratiosim)*100/rho21Ratiotheo
usErr = (usTheo- ussim)*100/usTheo



% figure('units','normalized','outerposition',[0 0 1 1])
% subplot(2,2,1)
% plot(x, rho, 'b--o')
% title('Density v/s x');

% subplot(2,2,2)
% plot(x, u, 'r--o')
% title('Velocity v/s x');

% subplot(2,2,3)
% plot(x, p, 'g--o')
% title('Pressure v/s x')

% subplot(2,2,4)
% plot(x, T, 'b--o')
% title('Temperature v/s x')

% d2u(m) = 0;
% figure('units','normalized','outerposition',[0 0 1 1])
% plot(x, d2u, 'b--o')
% title('d2ubydx2 v/s x')

% % saving the plots
% print('d2ubydx2','-dpng')	
