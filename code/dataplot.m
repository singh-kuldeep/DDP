% One gas plots 
Array=csvread('premetive.csv');
x = Array(:, 1);
rho = Array(:, 2);
u = Array(:, 3);
p = Array(:, 4);
T = p./(287.15*rho);

% Plots
% figure(1)
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
plot(x, rho, 'b--o')
title('Density v/s x')

subplot(2,2,2)
plot(x, u, 'r--o')
title('Velocity v/s x')

subplot(2,2,3)
plot(x, p, 'g--o')
title('Pressure v/s x')

subplot(2,2,4)
plot(x, T, 'b--o')
title('Temperature v/s x')

% One gas plots 
while(1)
	Array=csvread('premetive.csv');
	x = Array(:, 1);
	rho = Array(:, 2);
	u = Array(:, 3);
	p = Array(:, 4);
	T = p./(287.15*rho);
	
	pause(3)
	refresh()
end

% saving the plots
% print('plots','-dpng')	