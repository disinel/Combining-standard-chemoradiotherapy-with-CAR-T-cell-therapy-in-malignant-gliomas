clear variables;
close all;

% 1. Fixed Parameters and Constants
g1 = 10^10; g2 = 10^10; g3 = 2*10^9; 
K = 5*10^12; 
mu = 8.32; E0 = 1; D = 2;
alpha2 = 2.5*10^(-10); a2 = alpha2*K;  
gamma1 = g1/K; gamma2 = g2/K; gamma3 = g3/K;

% CAR-T dosing
v1 = 3*10^8; v2 = 3*10^8; v3 = 4*10^8;
V1 = v1/K; V2 = v2/K; V3 = v3/K;

% Simulation settings
Tfinal = 300;
grayColor = [.7 .7 .7];
tic

% TMZ cycle configuration: time interval in days between TMZ applications. 
% L11+1 represents the total number of TMZ applications per cycle. 
Tr = 7; Lr = 6; Tr1 = 1; Lr1 = 5; 
TgapSC = 7; TgapC = 7;

% Number of stored time points during ODE integration
Np = 100; 

% 2. Load Virtual Patient Cohort
load('Params_defenitve_with_alphaT_-2026-9-24-17-21_N=10000.mat');

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
gamma = median(gammaval); A = median(Aval); B = median(Bval); alphaT = median(alphaTval);

% 3. Run Simulation (Stupp(C)-3C-Stupp(A) protocol)
[tt, Y] = calc_1cycle(r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,V1,V2,V3,K,beta1,beta2,tau,A,B,gamma,D,Np,Tr, ...
                      Lr,Tr1,Lr1,TgapSC,TgapC,S10,S20,S30,Q10,Q20,Q30,E0,Tfinal,alphaT);
toc

% 4. Post-processing
% Re-dimensionalize variables
YY = K .* Y;

% Total tumor burden
Tc = K .* (Y(:,1) + Y(:,2) + Y(:,3) + Y(:,4) + Y(:,5) + Y(:,6) + Y(:,7));

% 5. Plotting Results
% Tumor Subpopulations Dynamics
f1 = figure();
plot(tt, YY(:,1), '-', 'color', '#77AC30', 'LineWidth', 1.5); hold on;
plot(tt, YY(:,3), '-', 'color', 'green', 'LineWidth', 1.5);
plot(tt, YY(:,5), '-', 'color', 'red', 'LineWidth', 1.5);
plot(tt, YY(:,7), '-', 'color', '#f5cb42', 'LineWidth', 1.5);
plot(tt, Tc, '-k', 'LineWidth', 1.5);

% Dynamic y-axis tick calculation
ymax1 = 5 * max(YY(:,5));
axis([min(tt) max(tt) min(YY(:,6)) ymax1]);
xticks([0 150 300]);

% Extract magnitude and determine maximum
mag1 = 10^floor(log10(ymax1));
max_tick1 = round(ymax1 / mag1, 1) * mag1;

% Generate 4 ticks rounded to the first decimal
ticks1 = round(linspace(0, max_tick1, 4) / mag1, 1) * mag1;
yticks(ticks1);

xlabel('t');
legend('S', 'R_{C}', 'R_{E}', 'D', 'T', 'Location', 'northwest');
fontsize(f1, 14, 'point');
fontname(f1, "Arial");
f1.WindowState = 'maximized';
saveas(f1, 'Stupp(C)_3C_Stupp(A)_tumor.eps', 'epsc');

% Quiescent Populations Dynamics
f2 = figure();
plot(tt, YY(:,2), '-', 'color', '#EDB120', 'LineWidth', 1.5); hold on;
plot(tt, YY(:,4), '-', 'color', '#EDB120', 'LineWidth', 1.5);
plot(tt, YY(:,6), '-', 'color', grayColor, 'LineWidth', 1.5);
axis([min(tt) max(tt) min(YY(:,6)) max(YY(:,2))]);
xlabel('t');
legend('Q_{s}', 'Q_{C}', 'Q_{E}', 'Location', 'best');
fontsize(f2, 14, 'point');
fontname(f2, "Arial");

% Treatments & Total Tumor
f3 = figure();
colororder({'black','black'});
plot(tt, Y(:,9), '-.', 'color', '#EDB120', 'LineWidth', 1.5); hold on;
plot(tt, YY(:,8)./10^(9), '-.', 'color', 'blue', 'LineWidth', 1.5);
xlim([0 Tfinal]);
xticks([0 150 300]);
yticks([0 1 2]);
ylim([0 2]);
xlabel('t');
ylabel('Normalized E, C', 'Interpreter', 'latex');

yyaxis right;
plot(tt, Tc, '-k', 'LineWidth', 1.5);
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
legend('E', 'C', 'T', 'Location', 'northeast', 'Interpreter', 'latex');
fontsize(f3, 14, 'point');
fontname(f3, "Arial");
f3.WindowState = 'maximized';
saveas(f3, 'Stupp(C)_3C_Stupp(A)_treatment.eps', 'epsc');

% Proliferative Fraction (Ki-67 tracking)
f4 = figure();
plot(tt, (YY(:,1) + YY(:,3) + YY(:,5))./Tc, '-k', 'LineWidth', 1.5); hold on;
yline(Ki, '-r', 'LineWidth', 1.5);
xlim([0 max(tt)]);
xlabel('t');
ylabel('U');
legend('U=X_{p}/T', 'Ki-67', 'Location', 'best');
fontsize(f4, 14, 'point');
fontname(f4, "Arial");

% =========================================================================
% FUNCTIONS
% =========================================================================

% ODE solver wrapper
function [t, y] = calc(ics, T1, T2, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np)
    tspan = linspace(T1, T2, Np);
    yy0 = ics;
    options = odeset('RelTol', 1e-10, 'NonNegative', [1,2,3,4,5,6,7,8,9]); 
    [t, y] = ode78(@(t, y) mODE(t, y, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu), tspan, yy0, options);
end

% Simulates the treatment sequence: Concomitant Stupp -> 3C -> Adjuvant Stupp
function [tt, Y] = calc_1cycle(r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, mu, V1, V2, V3, K, beta1, beta2, tau, A, B, gamma, D, Np, Tr, Lr, Tr1, Lr1, TgapSC, TgapC, S10, S20, S30, Q10, Q20, Q30, E0, Tfinal, alphaT)
    
    ics(1) = S10/K; ics(2) = Q10/K; ics(3) = S20/K; ics(4) = Q20/K; ics(5) = S30/K; ics(6) = Q30/K; ics(7) = 0; ics(8) = 0; ics(9) = 0;
    v = [V1 V2 V3];
    
    % Concomitant phase
    RTdays = makeRTschedule(Lr, Lr1, Tr);
    [days1, type1] = scheduleCTMZ(RTdays);
    TendCTMZ = max(days1);
    
    % Intermediate CAR-T block
    [days3, type3, dose3] = scheduleCART(TendCTMZ, TgapSC, TgapC, v);
    StartAdj = days3(end) + TgapSC;
    
    % Adjuvant phase
    [days2, type2] = scheduleATMZ(StartAdj);
    
    % Combine and sort all clinical events chronologically
    events.time = [days1 days2 days3];
    events.type = [type1 type2 type3];
    events.dose = [zeros(size(days1)) zeros(size(days2)) dose3];
    
    [events.time, idx] = sort(events.time);
    events.type = events.type(idx);
    events.dose = events.dose(idx);

    Sy = (length(events.time) + 1) * Np;
    Y = zeros(Sy, 9); 
    tt = zeros(Sy, 1);
    
    M = 0; t0 = 0;
    
    % Piecewise integration between treatment events
    for i = 1:length(events.time)
        [t, y] = calc(ics, t0, events.time(i), r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np);
        ics = applyTreatment(y(end,:), events.type{i}, events.dose(i), A, B, D, gamma, E0, alphaT);
        
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
function y = applyR(y, A, B, D, gamma, alphaT)
    SF = exp(-A*D - B*D^2);
    SFC = exp(-alphaT*D);
    
    y(7) = y(7) + (1-SF)*y(1) + (1-SF)*y(3) + (1-SF)*y(5);
    
    y(1) = SF*y(1) + gamma*y(2);
    y(2) = (1-gamma)*y(2);
    
    y(3) = SF*y(3) + gamma*y(4);
    y(4) = (1-gamma)*y(4);
    
    y(5) = SF*y(5) + gamma*y(6);
    y(6) = (1-gamma)*y(6);
    
    y(8) = SFC*y(8);
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
function y = applyTreatment(y, treatment, dose, A, B, D, gamma, E0, alphaT)
    switch treatment
        case 'RT_TMZ'
            y = applyR(y, A, B, D, gamma, alphaT);
            y(9) = y(9) + E0/3;
        case "TMZ_RT"
            y(9) = y(9) + E0/3;
        case "TMZ_ADJ"
            y(9) = y(9) + 2*E0/3;
        case "CART"
            y(8) = y(8) + dose;
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
function [days, type] = scheduleATMZ(StartDay)
    days = [];
    for cycle = 0:5
        cycleStart = StartDay + 28*cycle;
        days = [days, cycleStart + (0:4)];
    end
    type = repmat("TMZ_ADJ", size(days));
end

% Schedules CAR-T infusions
function [days, type, dose] = scheduleCART(TendStupp, TgapSC, TgapC, v)
    days = TendStupp + TgapSC + (0:length(v)-1)*TgapC;
    type = repmat("CART", size(days));
    dose = v;
end