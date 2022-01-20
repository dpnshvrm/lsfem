% Solving Stiff ODEs using lsfem
% Examples take from https://www.dam.brown.edu/people/alcyew/handouts/numODE5.pdf
%% Solving Robertson Chemical Reaction using lsfem
% ODE: x' = -0.04 x+ 1e4 y z)
%      y' = (0.04 x - 1e4 y z - 3e7 y^2)
%      z' = 3e7 y^2


clear, close all

model = 'Robertson''s Problem';
rhs = @(t,y) [ (-0.04*y(1) + 1e4*y(2)*y(3))
   (0.04*y(1) - 1e4*y(2)*y(3) - 3e7*y(2)^2)
   3e7*y(2)^2 ];

y0 = [1; 0; 0];

tspan = [0, 4e6];

% solve ode
fdmSol = ode15s(rhs,tspan,y0); % use ode15s
% fdmSol = ode45(ivp.odefun,ivp.tspan,ivp.y0); % use ode45
femSol = lsfem(rhs,tspan,y0); % use lsfem

% evaluate ODE solution
% t = linspace(tspan(1),tspan(end),1e4);
t = [0 4*logspace(-6,6,1000)];
yfdm = deval(fdmSol,t);
yfem = femSol.eval(t);

% uplift y 
yfdm(2,:) = 1e4*yfdm(2,:);
yfem(2,:) = 1e4*yfem(2,:);

% plot solutions
figure(1)
for i = 1:size(yfdm,1)
    subplot(3,1,i)
    semilogx(t,yfdm(i,:)), hold on
    semilogx(t,yfem(i,:))
    leg = {fdmSol.solver,femSol.solver};
    legend(leg, 'Location','Best')
    ylabel('$y^h$', 'Interpreter','Latex')
    xlabel('$t$', 'Interpreter','Latex')
%     axis([0 inf -1 1 ])
end
sgtitle(model) 

figure(2)
subplot(2,1,1)
semilogx(t,yfdm)
title([model,' using ',fdmSol.solver])
subplot(2,1,2)
semilogx(t,yfem)
title([model,' using ',femSol.solver])


%% Flame growth model
% ODE: y' = y.^2 - y.^3
clear, close all

model = 'Flame growth model';
rhs = @(t,y) y.^2 - y.^3;

y0 = 0.001;

tspan = [0, 2/y0];

% solve ode
fdmSol = ode15s(rhs,tspan,y0); % use ode15s
% fdmSol = ode45(ivp.odefun,ivp.tspan,ivp.y0); % use ode45
femSol = lsfem(rhs,tspan,y0); % use lsfem

% evaluate ODE solution
% t = linspace(tspan(1),tspan(end),1e4);
t = [0 2*logspace(-6,2,100)];
yfdm = deval(fdmSol,t);
yfem = femSol.eval(t);

% uplift y 
% yfdm(2,:) = 1e4*yfdm(2,:);
% yfem(2,:) = 1e4*yfem(2,:);

% plot solutions
figure(1)
plot(t,yfdm), hold on
plot(t,yfem)
leg = {fdmSol.solver,femSol.solver};
legend(leg, 'Location','Best')
ylabel('$y^h$', 'Interpreter','Latex')
xlabel('$t$', 'Interpreter','Latex')
title([model])





