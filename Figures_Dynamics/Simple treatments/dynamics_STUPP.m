clear variables;
close all;

% 1. Fixed Parameters and Constants
g1 = 10^10; g2 = 10^10; g3 = 2*10^9; 
v = 0.5*10^9; 
K = 5*10^12; 
mu = 8.32; E0 = 1; D = 2;
V = v/K; 
alpha2 = 2.5*10^(-10); a2 = alpha2*K;  
gamma1 = g1/K; gamma2 = g2/K; gamma3 = g3/K;

% Simulation settings
Tfinal = 200;
grayColor = [.7 .7 .7];
tic

% TMZ cycle configuration: time interval in days between TMZ applications. 
% L11+1 represents the total number of TMZ applications per cycle. 
Tr = 7; Lr = 6; Tr1 = 1; Lr1 = 5; 

% Number of stored time points during ODE integration
Np = 100; 
Sy = (Lr*(Lr1))*Np + Lr*Np + Np; 

% 2. Load Virtual Patient Cohort
load('Params_defenitve_-2026-8-25-17-35_N=10000.mat');

% Extract median values to simulate the Median Virtual Patient (MVP)
Tc0 = median(T0val); fs2 = median(delta2val); frc = median(delta1val); Ki = median(Kival);

% Initial conditions
S10 = Ki*Tc0*(1-frc-fs2); S20 = Ki*Tc0*frc; S30 = Ki*Tc0*fs2;  
Q10 = (1-Ki)*Tc0*(1-frc-fs2); Q20 = (1-Ki)*Tc0*frc; Q30 = (1-Ki)*Tc0*fs2;

% Kinetic rates
r1 = median(r1val); r2 = median(r2val); 
alpha1 = median(alpha1val); alpha3 = median(alpha3val); 
rho1 = median(rho1val); rho2 = median(rho2val); rho3 = median(rho3val); rho4 = median(rho4val); 
epsilon1 = median(epsilon1val); beta1 = median(beta1val); beta2 = median(beta2val); tau = median(tauval);
gamma = median(gammaval); A = median(Aval); B = median(Bval);

% 3. Run Simulation (Stupp Protocol: RT + TMZ)
[tt, Y] = calc_1cycle(r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,V,K,beta1,beta2,tau,A,B,gamma,D,Np,Tr, ...
                      Lr,Tr1,Lr1,S10,S20,S30,Q10,Q20,Q30,E0,Tfinal);
toc

% 4. Post-processing
% Re-dimensionalize variables
YY = K .* Y;

% Total tumor burden
Tc = K .* (Y(:,1) + Y(:,2) + Y(:,3) + Y(:,4) + Y(:,5) + Y(:,6) + Y(:,7));

Ynew = Y; YYnew = YY; ttnew = tt; 
save('new.mat','Ynew','YYnew','ttnew');

% 5. Plotting Results
% Tumor Subpopulations Dynamics
f = figure();
plot(tt, YY(:,1), '-', 'color', '#77AC30', 'LineWidth', 1.5); hold on;
plot(tt, YY(:,3), '-', 'color', 'green', 'LineWidth', 1.5);
plot(tt, YY(:,5), '-', 'color', 'red', 'LineWidth', 1.5);
plot(tt, YY(:,7), '-', 'color', '#f5cb42', 'LineWidth', 1.5);
plot(tt, Tc, '-black', 'LineWidth', 1.5); hold on;

% Dynamic y-axis tick calculation
ymax1 = max(Tc);
axis([min(tt) max(tt) min(YY(:,6)) 1.05*ymax1]);
xticks([0 100 200]);

% Extract magnitude and determine maximum
mag1 = 10^floor(log10(ymax1));
max_tick1 = round(ymax1 / mag1, 1) * mag1;

% Generate 4 ticks rounded to the first decimal
ticks1 = round(linspace(0, max_tick1, 4) / mag1, 1) * mag1;
yticks(ticks1);

xlabel('t');
ylabel('cells');
legend('S', 'R_{C}', 'R_{E}', 'D', 'T', 'Location', 'northeast');
fontsize(f, 14, 'point');
fontname(f, "Arial");
f.WindowState = 'maximized';
saveas(f, 'RT_TMZ_tumor.eps', 'epsc');

% Quiescent Populations Dynamics
f2 = figure();
plot(tt, YY(:,2), '-', 'color', '#77AC30', 'LineWidth', 1.5); hold on;
plot(tt, YY(:,4), '-', 'color', 'green', 'LineWidth', 1.5);
plot(tt, YY(:,6), '-', 'color', 'red', 'LineWidth', 1.5); hold on;

% Dynamic y-axis tick calculation
ymax2 = max(YY(:,2));
axis([min(tt) max(tt) min(YY(:,6)) ymax2]);
xticks([0 100 200]);

% Extract magnitude and determine maximum
mag2 = 10^floor(log10(ymax2));
max_tick2 = round(ymax2 / mag2, 1) * mag2;

% Generate and round intermediate ticks
ticks2 = round(linspace(0, max_tick2, 4) / mag2, 1) * mag2;
yticks(ticks2);

xlabel('t');
ylabel('cells');
legend('Q', 'Q_{C}', 'Q_{E}', 'Location', 'northeast');
fontsize(f2, 14, 'point');
fontname(f2, "Arial");
f2.WindowState = 'maximized';
saveas(f2, 'RT_TMZ_quiescent.eps', 'epsc');

% RT & TMZ Schedule & Total Tumor Dynamics
RTdays = makeRTschedule(Lr, Lr1, Tr);
f3 = figure();
colororder({'black','black'});
yyaxis left;
plot(tt, Y(:,9), '-.', 'color', '#EDB120', 'LineWidth', 1.5); hold on;

% Plot vertical lines for each RT application
for i = 1:length(RTdays)
    if i == 1
        plot([RTdays(i) RTdays(i)], [0 1], '-', 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'DisplayName', 'RT');
    else
        plot([RTdays(i) RTdays(i)], [0 1], '-', 'Color', [.7 .7 .7], 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
end
xticks([0 100 200]);
yticks([0 1]);
ylim([0 1]);
xlabel('t');
ylabel('RT, TMZ (E) applications');

yyaxis right;
plot(tt, Tc, '-black', 'LineWidth', 1.5, 'DisplayName', 'T');
ymax1 = max(Tc);
axis([min(tt) max(tt) min(YY(:,6)) 1.05*ymax1]);
xticks([0 100 200]);

% Extract magnitude and determine maximum
mag1 = 10^floor(log10(ymax1));
max_tick1 = round(ymax1 / mag1, 1) * mag1;

% Generate 4 ticks rounded to the first decimal
ticks1 = round(linspace(0, max_tick1, 4) / mag1, 1) * mag1;
yticks(ticks1);

ylabel('T');
legend('E', 'RT', 'T', 'Location', 'northeast');
fontsize(f3, 14, 'point');
fontname(f3, "Arial");
f3.WindowState = 'maximized';
saveas(f3, 'RT_TMZ_treatment.eps', 'epsc');

% =========================================================================
% FUNCTIONS
% =========================================================================

% ODE solver wrapper
function [t, y] = calc(ics, T1, T2, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np)
    tspan = linspace(T1, T2, Np);
    yy0 = ics;
    options = odeset('RelTol', 1e-10, 'NonNegative', [1,2,3,4,5]); 
    [t, y] = ode78(@(t, y) mODE(t, y, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu), tspan, yy0, options);
end

% Simulates the full treatment sequence (RT + TMZ)
function [tt, Y] = calc_1cycle(r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, mu, V, K, beta1, beta2, tau, A, B, gamma, D, Np, Tr, Lr, Tr1, Lr1, S10, S20, S30, Q10, Q20, Q30, E0, Tfinal)
    
    ics(1) = S10/K; ics(2) = Q10/K; ics(3) = S20/K; ics(4) = Q20/K; ics(5) = S30/K; ics(6) = Q30/K; ics(7) = 0; ics(8) = 0; ics(9) = 0;
    
    % Generate treatment schedules
    RTdays = makeRTschedule(Lr, Lr1, Tr);
    [days1, type1] = scheduleCTMZ(RTdays);
    [days2, type2] = scheduleATMZ(RTdays, 3);
    
    % Combine and sort all clinical events chronologically
    events.time = [days1 days2];
    events.type = [type1 type2];
    
    [events.time, idx] = sort(events.time);
    events.type = events.type(idx);
    
    % disp(events.time);
    % disp(events.type);
    
    Sy = (length(events.time) + 1) * Np;
    Y = zeros(Sy, 9); 
    tt = zeros(Sy, 1);
    
    M = 0; t0 = 0;
    
    % Piecewise integration between treatment events
    for i = 1:length(events.time)
        [t, y] = calc(ics, t0, events.time(i), r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np);
        ics = applyTreatment(y(end,:), events.type{i}, A, B, D, gamma, E0);
        
        tt(M*Np+1 : M*Np+Np) = t;
        Y(M*Np+1 : M*Np+Np, :) = y;
        M = M + 1;
        t0 = events.time(i);
    end
    
    % Final integration step until Tfinal
    [t, y] = calc(ics, t0, Tfinal, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np);
    tt(M*Np+1 : M*Np+Np) = t;
    Y(M*Np+1 : M*Np+Np, :) = y;
end

% Applies Radiation Therapy (RT) effects
function y = applyR(y, A, B, D, gamma)
    SF = exp(-A*D - B*D^2);
    
    y(7) = y(7) + (1-SF)*y(1) + (1-SF)*y(3) + (1-SF)*y(5);
    
    y(1) = SF*y(1) + gamma*y(2);
    y(2) = (1-gamma)*y(2);
    
    y(3) = SF*y(3) + gamma*y(4);
    y(4) = (1-gamma)*y(4);
    
    y(5) = SF*y(5) + gamma*y(6);
    y(6) = (1-gamma)*y(6);
    
    y(8) = SF*y(8);
end

% System of Ordinary Differential Equations
function dy = mODE(~, y, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu)
    dy = zeros(9, 1);
    
    % Total tumor burden logic
    T = y(1) + y(2) + y(3) + y(4) + y(5) + y(6) + y(7);
    
    dy(1) = r1*y(1)*(1 - T) - beta1*y(1) + beta2*y(2) - (alpha1 + epsilon1)*y(9)*y(1) - a2*y(8)*y(1);
    dy(2) = beta1*y(1) - beta2*y(2);
    dy(3) = r1*y(3)*(1 - T) - beta1*y(3) + beta2*y(4) - (alpha1 + epsilon1)*y(9)*y(3);
    dy(4) = beta1*y(3) - beta2*y(4);
    dy(5) = r2*y(5)*(1 - T) - beta1*y(5) + beta2*y(6) - a2*y(8)*y(5) + epsilon1*(y(1) + y(3))*y(9);
    dy(6) = beta1*y(5) - beta2*y(6);
    dy(7) = alpha1*(y(1) + y(3))*y(9) + a2*(y(1) + y(5))*y(8) - tau*y(7);
    dy(8) = -rho1*y(8) + (rho2.*y(1)*y(8))./(gamma1 + y(1)) + (rho3*y(5)*y(8))/(gamma2 + y(5)) ...
          - rho4*T*y(8)/(gamma3 + y(8)) - alpha3*y(9)*y(8);
    dy(9) = -mu*y(9);
end

% Generates Radiotherapy schedule
function RTdays = makeRTschedule(Lr, Lr1, Tr)
    RTdays = [];
    for week = 0:Lr-1
        RTdays = [RTdays, week*Tr + (1:Lr1)];
    end
end

% Evaluates treatment effects at discrete time points
function y = applyTreatment(y, treatment, A, B, D, gamma, E0)
    switch treatment
        case 'RT_TMZ'
            y = applyR(y, A, B, D, gamma);
            y(9) = y(9) + E0/3;
        case "TMZ_RT"
            y(9) = y(9) + E0/3;
        case "TMZ_ADJ"
            y(9) = y(9) + 2*E0/3;
    end
end

% Schedules Concomitant TMZ
function [days, type] = scheduleCTMZ(RTdays)
    Tmax = RTdays(end);
    days = 1:Tmax;
    type = strings(size(days));
    
    % Default: TMZ only
    type(:) = "TMZ_RT";
    % RT weekdays become RT+TMZ
    type(ismember(days, RTdays)) = "RT_TMZ";
end

% Schedules Adjuvant TMZ
function [days, type] = scheduleATMZ(RTdays, GapTMZ)
    startAdj = RTdays(end) + GapTMZ;
    days = [];
    for cycle = 0:5
        cycleStart = startAdj + 28*cycle;
        days = [days, cycleStart + (0:4)];
    end
    type = repmat("TMZ_ADJ", size(days));
end