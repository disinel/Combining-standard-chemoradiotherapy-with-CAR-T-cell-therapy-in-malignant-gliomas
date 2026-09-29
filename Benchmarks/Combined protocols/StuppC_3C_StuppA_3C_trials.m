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

% CAR-T dosing (half dose for split blocks)
v1 = 3*10^8 / 2; v2 = 3*10^8 / 2; v3 = 4*10^8 / 2;
V1 = v1/K; V2 = v2/K; V3 = v3/K;

% Simulation settings
Tfinal = 150000; % Extended time for survival tracking

% TMZ cycle configuration: time interval in days between TMZ applications. 
% L11+1 represents the total number of TMZ applications per cycle. 
Tr = 7; Lr = 6; Tr1 = 1; Lr1 = 5; 
TgapSC = 7; TgapC = 7;

% ODE integration settings
Np = 10; 

% 2. Load Virtual Patient Cohort (10,000 patients)
load('Params_defenitve_with_alphaT_-2026-9-24-17-21_N=10000.mat');

% Preallocate arrays for parallel computing
Tsurv = zeros(1, N); 
MinT = zeros(1, N);  

% Uncomment to initialize parallel pool if not already running
% delete(gcp('nocreate'));
% pool = parpool("Processes", 10);

% 3. Run Virtual Clinical Trial (Parallelized)
tic
parfor n = 1:N
    % Extract patient-specific parameters
    Tc0 = T0val(n); fs2 = delta2val(n); frc = delta1val(n); Ki = Kival(n);
    
    % Initial conditions for patient n
    S10 = Ki*Tc0*(1-frc-fs2); S20 = Ki*Tc0*frc; S30 = Ki*Tc0*fs2;  
    Q10 = (1-Ki)*Tc0*(1-frc-fs2); Q20 = (1-Ki)*Tc0*frc; Q30 = (1-Ki)*Tc0*fs2;
    
    % Kinetic rates for patient n
    r1 = r1val(n); r2 = r2val(n); 
    alpha1 = alpha1val(n); alpha3 = alpha3val(n); 
    rho1 = rho1val(n); rho2 = rho2val(n); rho3 = rho3val(n); rho4 = rho4val(n); 
    epsilon1 = epsilon1val(n); beta1 = beta1val(n); beta2 = beta2val(n); tau = tauval(n);
    gamma = gammaval(n); A = Aval(n); B = Bval(n); alphaT = alphaTval(n);
    
    % Compute the solution of the ODE system corresponding to the Stupp(C)-3C-Stupp(A)-3C protocol
    [te, minT] = calc_1cycle(r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,V1,V2,V3,K,beta1,beta2,tau,A,B,gamma,D,Np, ...
                             Tr,Lr,Tr1,Lr1,TgapSC,TgapC,S10,S20,Q10,Q20,S30,Q30,E0,Tfinal,alphaT);
                         
    % Store survival time and minimum tumor burden
    Tsurv(n) = te; 
    MinT(n) = minT * K;
end
toc

% 4. Post-processing and Statistics
% Identify patients achieving deep responses/minimal residual disease
[~, I1, Val1] = find(MinT <= 1);
[~, I10, Val10] = find(MinT <= 10);
[~, I100, Val100] = find(MinT <= 100);
[~, I1000, Val1000] = find(MinT <= 1000);
[~, I10000, Val10000] = find(MinT <= 10000);

N1 = length(I1);
N10 = length(I10);
N100 = length(I100);
N1000 = length(I1000);
N10000 = length(I10000);

% Save cohort survival results
save("Stupp(C)_3C_Stupp(A)_3C_Tsurv-" + year(datetime) + "-" + month(datetime) + "-" + day(datetime) + "-" + hour(datetime) + "-" + minute(datetime) + "_N=" + N + ".mat");

% Display Median Overall Survival (Months) safely omitting NaNs
disp('Median OS (months):');
disp(median(Tsurv, 'omitnan') / 30);

% Calculate and display 95% Confidence Intervals via Bootstrapping
nBoot = 10000;
bootstat = bootstrp(nBoot, @(x) median(x, 'omitnan'), Tsurv);
ci = prctile(bootstat, [2.5 97.5]);
disp('95% CI (months):');
disp(ci ./ 30);

% =========================================================================
% FUNCTIONS
% =========================================================================

% ODE solver wrapper with event detection (Lethal tumor burden threshold)
function [t, y, te, ye] = calc_event(ics, T1, T2, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np, K)
    
    tspan = linspace(T1, T2, Np);
    yy0 = ics;
    
    options = odeset('RelTol', 1e-10, 'AbsTol', 1e-16, 'Events', @events, 'NonNegative', [1,2,3,4,5,6,7,8,9]);
    
    function [position, isterminal, direction] = events(~, y)
        % Integration terminates when tumor burden reaches the lethal threshold (K/5)
        position = K.*(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)) - K/5;                                           
        isterminal = 1;                                                        
        direction = 0;                                                         
    end

    [t, y, te, ye] = ode78(@(t, y) mODE(t, y, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu), tspan, yy0, options);
end

% Simulates the full treatment sequence: Concomitant Stupp -> 3C -> Adjuvant Stupp -> 3C
function [te, minT] = calc_1cycle(r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, mu, V1, V2, V3, K, beta1, beta2, tau, A, B, gamma, D, Np, Tr, Lr, Tr1, Lr1, TgapSC, TgapC, S10, S20, Q10, Q20, S30, Q30, E0, Tfinal, alphaT)
    
    ics(1) = S10/K; ics(2) = Q10/K; ics(3) = S20/K; ics(4) = Q20/K; ics(5) = S30/K; ics(6) = Q30/K; ics(7) = 0; ics(8) = 0; ics(9) = E0;
    v = [V1 V2 V3];
    
    % Concomitant phase schedule
    RTdays = makeRTschedule(Lr, Lr1, Tr);
    [days1, type1] = scheduleCTMZ(RTdays);
    TendCTMZ = max(days1);
    
    % First split CAR-T block
    [days3, type3, dose3] = scheduleCART(TendCTMZ, TgapSC, TgapC, v);
    StartAdj = days3(end) + TgapSC;
    
    % Adjuvant phase schedule
    [days2, type2] = scheduleATMZ(StartAdj);
    TendAdj = max(days2);
    
    % Second split CAR-T block
    [days4, type4, dose4] = scheduleCART(TendAdj, TgapSC, TgapC, v);
    
    % Combine and sort all clinical events chronologically
    events.time = [days1 days2 days3 days4];
    events.type = [type1 type2 type3 type4];
    events.dose = [zeros(size(days1)) zeros(size(days2)) dose3 dose4];
    
    [events.time, idx] = sort(events.time);
    events.type = events.type(idx);
    events.dose = events.dose(idx);
    
    Sy = (length(events.time) + 1) * Np;
    Y = zeros(Sy, 9); 
    tt = zeros(Sy, 1);
    
    M = 0; t0 = 0;
    minT = inf;
    
    % Piecewise integration between treatment events
    for i = 1:length(events.time)
        [t, y, te] = calc_event(ics, t0, events.time(i), r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np, K);
        
        minT = min(minT, min(sum(y(:,1:6), 2)));
        if ~isempty(te); return; end
        
        ics = applyTreatment(y(end,:), events.type{i}, events.dose(i), A, B, D, gamma, E0, alphaT);
        tt(M*Np+1 : M*Np+Np) = t;
        Y(M*Np+1 : M*Np+Np, :) = y;
        M = M + 1;
        t0 = events.time(i);
    end
    
    % Final integration step until Tfinal or lethal event
    [t, y, te] = calc_event(ics, t0, Tfinal, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np, K);
    minT = min(minT, min(sum(y(:,1:6), 2)));
    
    if ~isempty(te); return; end
    
    tt(M*Np+1 : M*Np+length(t)) = t;
    Y(M*Np+1 : M*Np+length(t), :) = y;
    
    % Robustness check for patients surviving beyond Tfinal
    if isempty(te)
        te = NaN;
    end
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