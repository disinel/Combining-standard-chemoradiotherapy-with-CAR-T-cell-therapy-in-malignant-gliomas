clear variables;
close all;

% 1. Fixed Parameters and Constants
g1 = 10^10; g2 = 10^10; g3 = 2*10^9; 
K = 5*10^12; 
mu = 8.32; E0 = 1; D = 2;
alpha2 = 2.5*10^(-10); a2 = alpha2*K;  
gamma1 = g1/K; gamma2 = g2/K; gamma3 = g3/K;

% End time of the integration
Tfinal = 150000; % Extended time for survival tracking
% Number of stored time points during ODE integration
Np = 10; 

% 2. Load Virtual Patient Cohort (10,000 patients)
load('Params_defenitve_with_alphaT_-2026-9-24-17-21_N=10000.mat');

% Uncomment to initialize parallel pool if not already running
% delete(gcp('nocreate'));
% pool = parpool("Processes", 10);

% =========================================================================
% 3. DEFINITION OF STUPP - MC / MC - STUPP PROTOCOLS (M = 1 to 10)
% =========================================================================
M_max = 10;
num_protocols = 2 * M_max; 
protocols_time = cell(num_protocols, 1);
protocols_type = cell(num_protocols, 1);
protocols_dose = cell(num_protocols, 1);
protocol_names = strings(num_protocols, 1);
idx = 1;

for m = 1:M_max
    % Stupp - MC Protocol
    [pt, pty, pd] = build_protocol(m, 'Stupp-MC', K);
    protocols_time{idx} = pt; 
    protocols_type{idx} = pty; 
    protocols_dose{idx} = pd;
    protocol_names(idx) = "Stupp - " + m + "C";
    idx = idx + 1;
    
    % MC - Stupp Protocol
    [pt, pty, pd] = build_protocol(m, 'MC-Stupp', K);
    protocols_time{idx} = pt; 
    protocols_type{idx} = pty; 
    protocols_dose{idx} = pd;
    protocol_names(idx) = m + "C - Stupp";
    idx = idx + 1;
end

% Data matrix preallocation
Tsurv = zeros(num_protocols, N); 
MinT = zeros(num_protocols, N);  

% =========================================================================
% 4. MAIN VIRTUAL CLINICAL TRIAL (Parallelized)
% =========================================================================
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
    
    ics = [S10/K, Q10/K, S20/K, Q20/K, S30/K, Q30/K, 0, 0, E0];
    
    t_surv_temp = zeros(num_protocols, 1);
    min_t_temp = zeros(num_protocols, 1);
    
    for p = 1:num_protocols
        evt_time = protocols_time{p};
        evt_type = protocols_type{p};
        evt_dose = protocols_dose{p};
        
        [te, min_val] = calc_1cycle(ics, evt_time, evt_type, evt_dose, ...
                                   r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, ...
                                   epsilon1, gamma1, gamma2, gamma3, mu, K, beta1, beta2, ...
                                   tau, A, alphaT, B, gamma, D, Np, E0, Tfinal);
                               
        t_surv_temp(p) = te;
        min_t_temp(p)  = min_val * K;
    end
    
    Tsurv(:, n) = t_surv_temp;
    MinT(:, n)  = min_t_temp;
end
toc

% =========================================================================
% 5. STATISTICS AND SAVING RESULTS
% =========================================================================
median_Tsurv = median(Tsurv, 2) / 30;
disp("Median Overall Survival (months) per protocol:");
for p = 1:num_protocols
    fprintf('%-15s: %.2f months\n', protocol_names(p), median_Tsurv(p));
end

save("Stupp_CART_Combinations_2V_Tsurv-" + year(datetime) + "-" + month(datetime) + "-" + day(datetime) + ...
    "-" + hour(datetime) + "-" + minute(datetime) + "_N=" + N + ".mat", ...
    "Tsurv", "MinT", "protocol_names", "M_max");

% =========================================================================
% FUNCTIONS
% =========================================================================

% Protocol builder for fractionated CAR-T (M doses) and Stupp combinations
function [evt_time, evt_type, evt_dose] = build_protocol(M, order, K)
    Tr = 7; Lr = 6; Lr1 = 5; 
    Tgap = 7; 
    
    % Total scaled CAR-T dose adjusted and divided into M fractions
    total_dose_scaled = 2 * 10^9 / K;
    dose_m = total_dose_scaled / M;
    
    if strcmp(order, 'Stupp-MC')
        % Stupp schedule starting on day 1
        offset = 0;
        RTdays = offset + make_base_RT(Lr, Lr1, Tr);
        [d_cTMZ, t_cTMZ] = scheduleCTMZ(RTdays);
        [d_aTMZ, t_aTMZ] = scheduleATMZ(RTdays, 3);
        
        TendStupp = RTdays(end) + 3 + 6*28;
        
        % CAR-T fractions administered after the post-Stupp gap
        d_CART = TendStupp + Tgap + (0 : M-1) * Tgap;
        t_CART = repmat("CART", 1, M);
        dose_CART = repmat(dose_m, 1, M);
        
    elseif strcmp(order, 'MC-Stupp')
        % CAR-T fractions starting on day 1
        d_CART = 1 + (0 : M-1) * Tgap;
        t_CART = repmat("CART", 1, M);
        dose_CART = repmat(dose_m, 1, M);
        
        CART_end = d_CART(end);
        
        % Stupp schedule following CAR-T completion
        offset = CART_end + Tgap - 1; 
        RTdays = offset + make_base_RT(Lr, Lr1, Tr);
        [d_cTMZ, t_cTMZ] = scheduleCTMZ(RTdays);
        [d_aTMZ, t_aTMZ] = scheduleATMZ(RTdays, 3);
    end
    
    evt_time = [d_cTMZ, d_aTMZ, d_CART];
    evt_type = [t_cTMZ, t_aTMZ, t_CART];
    evt_dose = [zeros(size(d_cTMZ)), zeros(size(d_aTMZ)), dose_CART];
    
    [evt_time, idx] = sort(evt_time);
    evt_type = evt_type(idx);
    evt_dose = evt_dose(idx);
end

% Generates Radiotherapy base schedule
function RTdays = make_base_RT(Lr, Lr1, Tr)
    RTdays = [];
    for week = 0:Lr-1
        RTdays = [RTdays, week*Tr + (1:Lr1)];
    end
end

% Schedules Concomitant TMZ
function [days, type] = scheduleCTMZ(RTdays)
    days = RTdays(1) : RTdays(end);
    type = repmat("TMZ_RT", size(days));
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

% Simulates full protocol cycle for a single patient
function [te, minT] = calc_1cycle(ics, evt_time, evt_type, evt_dose, ...
                                  r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, ...
                                  epsilon1, gamma1, gamma2, gamma3, mu, K, beta1, beta2, ...
                                  tau, A, alphaT, B, gamma, D, Np, E0, Tfinal)
    M_idx = 0; t0 = 0;
    minT = inf;
    for i = 1:length(evt_time)
        [t, y, te] = calc_event(ics, t0, evt_time(i), r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np, K);
        minT = min(minT, min(sum(y(:,1:6),2)));
        
        if ~isempty(te) 
            return; 
        end
        ics = applyTreatment(y(end,:), evt_type(i), evt_dose(i), A, alphaT, B, D, gamma, E0);
        M_idx = M_idx + 1;
        t0 = evt_time(i);
    end
    
    [t, y, te] = calc_event(ics, t0, Tfinal, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np, K);
    minT = min(minT, min(sum(y(:,1:6),2)));
end

% ODE solver wrapper with event detection (Lethal tumor burden threshold)
function [t, y, te, ye] = calc_event(ics, T1, T2, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np, K)
    tspan = linspace(T1, T2, Np);
    options = odeset('RelTol', 1e-10, 'Events', @events, 'NonNegative', [1,2,3,4,5,6,7,8,9]); 
    
    function [position, isterminal, direction] = events(~, y)
        position = K.*(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)) - K/5;                                           
        isterminal = 1;                                                        
        direction = 0;                                                         
    end
    [t, y, te, ye] = ode78(@(t,y)mODE(t,y,r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,beta1,beta2,tau,mu), tspan, ics, options);
end

% System of Ordinary Differential Equations
function dy = mODE(~,y,r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,beta1,beta2,tau,mu)
    dy = zeros(9,1);
    dy(1) = r1*y(1)*(1-y(1)-y(2)-y(3)-y(4)-y(5)-y(6)-y(7))-beta1*y(1)+beta2*y(2)-(alpha1+epsilon1)*y(9)*y(1)-a2*y(8)*y(1);
    dy(2) = beta1*y(1)-beta2*y(2);
    dy(3) = r1*y(3)*(1-y(1)-y(2)-y(3)-y(4)-y(5)-y(6)-y(7))-beta1*y(3)+beta2*y(4)-(alpha1+epsilon1)*y(9)*y(3);
    dy(4) = beta1*y(3)-beta2*y(4);
    dy(5) = r2*y(5)*(1-y(1)-y(2)-y(3)-y(4)-y(5)-y(6)-y(7))-beta1*y(5)+beta2*y(6)-a2*y(8)*y(5)+epsilon1*(y(1)+y(3))*y(9);
    dy(6) = beta1*y(5)-beta2*y(6);
    dy(7) = alpha1*(y(1)+y(3))*y(9)+a2*(y(1)+y(5))*y(8)-tau*y(7);
    dy(8) = -rho1*y(8)+(rho2.*y(1)*y(8))./(gamma1+y(1))+(rho3*y(5)*y(8))/(gamma2+y(5)) - rho4*(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7))*y(8)/(gamma3+y(8))-alpha3*y(9)*y(8);
    dy(9) = -mu*y(9);
end

% Applies Radiation Therapy (RT) effects
function y = applyR(y, A, alphaT, B, D, gamma)
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
function y = applyTreatment(y,treatment,dose,A,alphaT,B,D,gamma,E0)
    switch treatment
        case 'RT_TMZ'
            y = applyR(y,A,alphaT,B,D,gamma);
            y(9) = y(9) + E0/3;
        case "TMZ_RT"
            y(9) = y(9) + E0/3;
        case "TMZ_ADJ"
            y(9) = y(9) + 2*E0/3;
        case "CART"
            y(8) = y(8) + dose;
    end
end