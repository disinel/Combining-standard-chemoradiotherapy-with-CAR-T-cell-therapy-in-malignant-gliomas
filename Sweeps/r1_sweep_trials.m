clear variables;
close all;

% 1. Fixed Parameters
g1 = 10^10; g2 = 10^10; g3 = 2*10^9; 
K = 5*10^12; mu = 8.32; E0 = 1; D = 2;
alpha2 = 2.5*10^(-10); a2 = alpha2*K;  
gamma1 = g1/K; gamma2 = g2/K; gamma3 = g3/K;

% Time Configuration
Tfinal = 150000;
Tr = 7; Lr = 6; Tr1 = 1; Lr1 = 5; 
Np = 10; 

% CAR-T Dose Protocol 1: 3C-Stupp-3C (Half dose)
v1_P1 = 3*10^8/2; v2_P1 = 3*10^8/2; v3_P1 = 4*10^8/2;
V_P1 = [v1_P1/K, v2_P1/K, v3_P1/K];

% CAR-T Dose Protocol 2: CAR-T Only (Full dose)
v1_P2 = 3*10^8; v2_P2 = 3*10^8; v3_P2 = 4*10^8;
V_P2 = [v1_P2/K, v2_P2/K, v3_P2/K];

% 2. Load SINGLE Base Cohort
% (This file must contain the r1-independent variables and/or base seeds/vectors)
load('Params_defenitve_with_alphaT_-2026-9-24-17-21_N=10000.mat'); 

% r1_max values to explore
r1_test_values = 0.15 : 0.01 : 0.25;
num_tests = length(r1_test_values);

% Data Matrices Preallocation for the 3 Protocols
All_Tsurv_3C = zeros(num_tests, N); All_MinT_3C = zeros(num_tests, N);
All_Tsurv_CART = zeros(num_tests, N); All_MinT_CART = zeros(num_tests, N);
All_Tsurv_NoTreat = zeros(num_tests, N); All_MinT_NoTreat = zeros(num_tests, N);

disp('Starting unified massive simulation...');
disp('-------------------------------------------------------------------------');

% =========================================================================
% 3. EXTERNAL LOOP (r1 Sweep)
% =========================================================================
for idx = 1:num_tests
    r1_max = r1_test_values(idx);
    disp(['--> Processing r1_max = ', num2str(r1_max), ' (Iteration ', num2str(idx), '/', num2str(num_tests), ')']);
    
    % ---------------------------------------------------------------------
    r1_distribucion_base = r1val / max(r1val);
    
    r1val_current = r1_distribucion_base * r1_max;
    
    r2val_current = 0.5 * r1val_current;
    
    beta2val_current = Kival .* (beta1val ./ (1 - Kival) - r1val_current);
    tauval_current = unifrnd(r1val_current ./ 100, 2 .* r1val_current, [1, N]);

    % Temporary matrices for parfor
    Tsurv_3C = zeros(1, N); MinT_3C = zeros(1, N);
    Tsurv_CART = zeros(1, N); MinT_CART = zeros(1, N);
    Tsurv_No = zeros(1, N); MinT_No = zeros(1, N);
    
    % =====================================================================
    % 4. INTERNAL PARALLEL LOOP (Patients)
    % =====================================================================
    tic
    parfor n = 1:N
        % Initial conditions for patient n
        Tc0 = T0val(n); fs2 = delta2val(n); frc = delta1val(n); Ki = Kival(n);
        S10 = Ki*Tc0*(1-frc-fs2); S20 = Ki*Tc0*frc; S30 = Ki*Tc0*fs2;  
        Q10 = (1-Ki)*Tc0*(1-frc-fs2); Q20 = (1-Ki)*Tc0*frc; Q30 = (1-Ki)*Tc0*fs2;
        ics = [S10/K, Q10/K, S20/K, Q20/K, S30/K, Q30/K, 0, 0, E0];
        
        % Dynamic and kinetic parameters of the patient
        r1 = r1val_current(n); r2 = r2val_current(n); 
        alpha1 = alpha1val(n); alpha3 = alpha3val(n); alphaT = alphaTval(n);
        rho1 = rho1val(n); rho2 = rho2val(n); rho3 = rho3val(n); rho4 = rho4val(n); 
        epsilon1 = epsilon1val(n); beta1 = beta1val(n); 
        beta2 = beta2val_current(n); tau = tauval_current(n);
        
        gamma = gammaval(n); A = Aval(n); B = Bval(n);
        % --- Protocol 1: 3C - Stupp - 3C ---
        [te_3c, mt_3c] = calc_3C_Stupp_3C(ics, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, ...
            epsilon1, gamma1, gamma2, gamma3, mu, V_P1, K, beta1, beta2, tau, A, alphaT, B, gamma, D, Np, ...
            Tr, Lr, Tr1, Lr1, 7, 7, E0, Tfinal);
        Tsurv_3C(n) = te_3c; MinT_3C(n) = mt_3c * K;
        
        % --- Protocol 2: CAR-T Only ---
        [te_cart, mt_cart] = calc_CART_Only(ics, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, ...
            epsilon1, gamma1, gamma2, gamma3, mu, V_P2, K, beta1, beta2, tau, A, alphaT, B, gamma, D, Np, ...
            E0, Tfinal, 0, 30);
        Tsurv_CART(n) = te_cart; MinT_CART(n) = mt_cart * K;
        
        % --- Protocol 3: No Treatment ---
        [te_no, mt_no] = calc_No_Treat(ics, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, ...
            epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np, K, Tfinal);
        Tsurv_No(n) = te_no; MinT_No(n) = mt_no * K;
    end
    toc
    
    % Save results in global matrices
    All_Tsurv_3C(idx, :) = Tsurv_3C; All_MinT_3C(idx, :) = MinT_3C;
    All_Tsurv_CART(idx, :) = Tsurv_CART; All_MinT_CART(idx, :) = MinT_CART;
    All_Tsurv_NoTreat(idx, :) = Tsurv_No; All_MinT_NoTreat(idx, :) = MinT_No;
    
    % Print metrics
    disp(['    Median OS (months) -> 3C-S-3C: ', num2str(median(Tsurv_3C)/30), ...
          ' | CAR-T: ', num2str(median(Tsurv_CART)/30), ...
          ' | No Treatment: ', num2str(median(Tsurv_No)/30)]);
end

% =========================================================================
% 5. UNIFIED FINAL SAVE
% =========================================================================
save("Unified_r1_Sweep_Tsurv-"+year(datetime)+"-"+month(datetime)+"-"+day(datetime)+"-"+hour(datetime)+"-"+minute(datetime)+".mat", ...
    "All_Tsurv_3C", "All_Tsurv_CART", "All_Tsurv_NoTreat", ...
    "All_MinT_3C", "All_MinT_CART", "All_MinT_NoTreat", "r1_test_values", "N");
disp('Sweep completed and saved successfully!');

% =========================================================================
%                        UNIFIED AUXILIARY FUNCTIONS
% =========================================================================

function [t, y, te, ye] = calc_event(ics, T1, T2, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np, K)
    tspan = linspace(T1, T2, Np);
    options = odeset('RelTol', 1e-10, 'Events', @events, 'NonNegative', [1,2,3,4,5,6,7,8,9]); 
    function [position, isterminal, direction] = events(~, y)
        position = K.*(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)) - K/5;                                           
        isterminal = 1;                                                        
        direction = 0;                                                         
    end
    [t, y, te, ye] = ode78(@(t, y) mODE(t, y, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu), tspan, ics, options);
end

% Simulator Protocol 1
function [te, minT] = calc_3C_Stupp_3C(ics, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, mu, v, K, beta1, beta2, tau, A, alphaT, B, gamma, D, Np, Tr, Lr, Tr1, Lr1, TgapSC, TgapC, E0, Tfinal)
    [days3, type3, dose3] = scheduleCART(1, 0, TgapC, v);
    StuppStart = days3(end) + TgapSC;
    RTdays = makeRTschedule(Lr, Lr1, Tr) + (StuppStart - 1);
    [days1, type1] = scheduleCTMZ(RTdays);
    [days2, type2] = scheduleATMZ(RTdays, 3);
    TendStupp = max([days1 days2]);
    [days4, type4, dose4] = scheduleCART(TendStupp, TgapSC, TgapC, v);
    
    events.time = [days1 days2 days3 days4];
    events.type = [type1 type2 type3 type4];
    events.dose = [zeros(size(days1)) zeros(size(days2)) dose3 dose4];
    [events.time, idx] = sort(events.time);
    events.type = events.type(idx);
    events.dose = events.dose(idx);
    
    [te, minT] = run_integration_cycle(ics, events, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, mu, K, beta1, beta2, tau, A, alphaT, B, gamma, D, Np, E0, Tfinal);
end

% Simulator Protocol 2
function [te, minT] = calc_CART_Only(ics, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, mu, v, K, beta1, beta2, tau, A, alphaT, B, gamma, D, Np, E0, Tfinal, TgapSC, TgapC)
    [days3, type3, dose3] = scheduleCART(1, TgapSC, TgapC, v);
    events.time = days3;
    events.type = type3;
    events.dose = dose3;
    
    [te, minT] = run_integration_cycle(ics, events, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, mu, K, beta1, beta2, tau, A, alphaT, B, gamma, D, Np, E0, Tfinal);
end

% Simulator Protocol 3
function [te, minT] = calc_No_Treat(ics, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np, K, Tfinal)
    [t, y, te] = calc_event(ics, 0, Tfinal, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np, K);
    minT = min(sum(y(:,1:6), 2));
end

% Shared Integration Engine
function [te, minT] = run_integration_cycle(ics, events, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, mu, K, beta1, beta2, tau, A, alphaT, B, gamma, D, Np, E0, Tfinal)
    t0 = 0; minT = inf;
    for i = 1:length(events.time)
        if events.time(i) > t0
            [t, y, te] = calc_event(ics, t0, events.time(i), r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np, K);
            minT = min(minT, min(sum(y(:,1:6), 2)));
            if ~isempty(te); return; end
            
            ics = y(end,:); 
        end
        
        ics = applyTreatment(ics, events.type(i), events.dose(i), A, alphaT, B, D, gamma, E0);
        t0 = events.time(i);
    end
    
    if Tfinal > t0
        [t, y, te] = calc_event(ics, t0, Tfinal, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np, K);
        minT = min(minT, min(sum(y(:,1:6), 2)));
    end
end

% Ordinary Differential Equations
function dy = mODE(~, y, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu)
    dy = zeros(9, 1);
    dy(1) = r1*y(1)*(1-y(1)-y(2)-y(3)-y(4)-y(5)-y(6)-y(7))-beta1*y(1)+beta2*y(2)-(alpha1+epsilon1)*y(9)*y(1)-a2*y(8)*y(1);
    dy(2) = beta1*y(1)-beta2*y(2);
    dy(3) = r1*y(3)*(1-y(1)-y(2)-y(3)-y(4)-y(5)-y(6)-y(7))-beta1*y(3)+beta2*y(4)-(alpha1+epsilon1)*y(9)*y(3);
    dy(4) = beta1*y(3)-beta2*y(4);
    dy(5) = r2*y(5)*(1-y(1)-y(2)-y(3)-y(4)-y(5)-y(6)-y(7))-beta1*y(5)+beta2*y(6)-a2*y(8)*y(5)+epsilon1*(y(1)+y(3))*y(9);
    dy(6) = beta1*y(5)-beta2*y(6);
    dy(7) = alpha1*(y(1)+y(3))*y(9)+a2*(y(1)+y(5))*y(8)-tau*y(7);
    dy(8) = -rho1*y(8)+(rho2.*y(1)*y(8))./(gamma1+y(1))+(rho3*y(5)*y(8))/(gamma2+y(5)) ...
          - rho4*(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7))*y(8)/(gamma3+y(8))-alpha3*y(9)*y(8);
    dy(9) = -mu*y(9);
end

% Event and Protocol Generators
function RTdays = makeRTschedule(Lr, Lr1, Tr)
    RTdays = [];
    for week = 0:Lr-1; RTdays = [RTdays, week*Tr + (1:Lr1)]; end
end

function [days, type] = scheduleCTMZ(RTdays)
    days = RTdays(1) : RTdays(end);
    type = strings(size(days));
    type(:) = "TMZ_RT";
    type(ismember(days, RTdays)) = "RT_TMZ";
end

function [days, type] = scheduleATMZ(RTdays, GapTMZ)
    startAdj = RTdays(end) + GapTMZ;
    days = [];
    for cycle = 0:5; cycleStart = startAdj + 28*cycle; days = [days, cycleStart + (0:4)]; end
    type = repmat("TMZ_ADJ", size(days));
end

function [days, type, dose] = scheduleCART(TendStupp, TgapSC, TgapC, v)
    days = TendStupp + TgapSC + (0:length(v)-1)*TgapC;
    type = repmat("CART", size(days));
    dose = v;
end

function y = applyR(y, A, alphaT, B, D, gamma)
    SF = exp(-A*D - B*D^2);
    SFC = exp(-alphaT*D);
    y(7) = y(7) + (1-SF)*y(1) + (1-SF)*y(3)+ (1-SF)*y(5);
    y(1) = SF*y(1) + gamma*y(2); y(2) = (1-gamma)*y(2);
    y(3) = SF*y(3) + gamma*y(4); y(4) = (1-gamma)*y(4);
    y(5) = SF*y(5) + gamma*y(6); y(6) = (1-gamma)*y(6);
    y(8) = SFC*y(8);
end

function y = applyTreatment(y, treatment, dose, A, alphaT, B, D, gamma, E0)
    switch treatment
        case 'RT_TMZ'; y = applyR(y, A, alphaT, B, D, gamma); y(9) = y(9) + E0/3;
        case "TMZ_RT"; y(9) = y(9) + E0/3;
        case "TMZ_ADJ"; y(9) = y(9) + 2*E0/3;
        case "CART"; y(8) = y(8) + dose;
    end
end