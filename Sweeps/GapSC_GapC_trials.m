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
Np = 10; 

% 2. Load Virtual Patient Cohort (10,000 patients) - Using the alphaT version!
load('Params_defenitve_with_alphaT_-2026-9-24-17-21_N=10000.mat');

% Uncomment to initialize parallel pool if not already running
% delete(gcp('nocreate'));
% pool = parpool("Processes", 10);

% =========================================================================
% 3. DEFINITION OF THE UNIVERSAL SWEEP GRID (270 combinations)
% =========================================================================
P_vals     = ["Stupp-3C", "3C-Stupp", "StuppC-3C-StuppA-3C"];
D_vals     = [10^9, 2*10^9];                    % Total dose
gapSC_vals = [0, 7, 14, 21, 28];                % Gap between major blocks (Stupp/CART)
gapC_vals  = [1, 2, 3, 4, 5, 6, 7, 10, 14];     % Gap between CAR-T injections

% Generate all possible combinations
[P_grid, D_grid, GSC_grid, GC_grid] = ndgrid(1:length(P_vals), D_vals, gapSC_vals, gapC_vals);
P_idx_grid = P_grid(:);
D_grid = D_grid(:);
GSC_grid = GSC_grid(:);
GC_grid = GC_grid(:);

num_protocols = length(D_grid);
protocols_time = cell(num_protocols, 1);
protocols_type = cell(num_protocols, 1);
protocols_dose = cell(num_protocols, 1);
protocol_names = strings(num_protocols, 1);

% Build schedules for each combination
for p = 1:num_protocols
    prot  = P_vals(P_idx_grid(p));
    dtot  = D_grid(p);
    t_sc  = GSC_grid(p);
    t_c   = GC_grid(p);
    
    [pt, pty, pd] = build_universal_protocol(prot, dtot, t_sc, t_c, K);
    protocols_time{p} = pt; 
    protocols_type{p} = pty; 
    protocols_dose{p} = pd;
    protocol_names(p) = sprintf("%s | D: %1.0e | GapSC: %dd | GapC: %dd", prot, dtot, t_sc, t_c);
end

% =========================================================================
% 4. DATA MATRICES PREALLOCATION
% =========================================================================
Tsurv = zeros(num_protocols, N); 
MinT  = zeros(num_protocols, N);  

% =========================================================================
% 5. MAIN VIRTUAL CLINICAL TRIAL (Parallelized Universal Evaluation)
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
    min_t_temp  = zeros(num_protocols, 1);
    
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
% 6. STATISTICS AND SAVING RESULTS
% =========================================================================
disp("=================================================================");
disp("Protocol (Gaps)                                     | Median OS (m)");
disp("-----------------------------------------------------------------");
median_Tsurv = median(Tsurv, 2) / 30;
for p = 1:num_protocols
    fprintf('%-51s | %7.2f \n', protocol_names(p), median_Tsurv(p));
end
disp("=================================================================");

save("Universal_Gaps_Sweep_Tsurv-" + year(datetime) + "-" + month(datetime) + "-" + day(datetime) + ...
    "-" + hour(datetime) + "-" + minute(datetime) + "_N=" + N + ".mat", ...
    "Tsurv", "MinT", "protocol_names", "P_idx_grid", "D_grid", "GSC_grid", "GC_grid", "P_vals");
disp("Simulation finished. Results matrix successfully saved.");

% =========================================================================
% FUNCTIONS
% =========================================================================

% Universal Protocol Builder
function [evt_time, evt_type, evt_dose] = build_universal_protocol(prot, D_tot, TgapSC, TgapC, K)
    Tr = 7; Lr = 6; Lr1 = 5; 
    
    if prot == "Stupp-3C"
        M = 3;
        dose_m = (D_tot / M) / K;
        
        RTdays = make_base_RT(Lr, Lr1, Tr);
        [d_cTMZ, t_cTMZ] = scheduleCTMZ(RTdays);
        StartAdj = RTdays(end) + 3;
        [d_aTMZ, t_aTMZ] = scheduleATMZ_from_start(StartAdj);
        
        TendStupp = max(d_aTMZ);
        d_CART = TendStupp + TgapSC + (0:M-1)*TgapC;
        t_CART = repmat("CART", 1, M);
        dose_CART = repmat(dose_m, 1, M);
        
        evt_time = [d_cTMZ, d_aTMZ, d_CART];
        evt_type = [t_cTMZ, t_aTMZ, t_CART];
        evt_dose = [zeros(size(d_cTMZ)), zeros(size(d_aTMZ)), dose_CART];
        
    elseif prot == "3C-Stupp"
        M = 3;
        dose_m = (D_tot / M) / K;
        
        d_CART = 1 + (0:M-1)*TgapC;
        TendCART = d_CART(end);
        t_CART = repmat("CART", 1, M);
        dose_CART = repmat(dose_m, 1, M);
        
        offset = TendCART + TgapSC - 1;
        RTdays = offset + make_base_RT(Lr, Lr1, Tr);
        
        [d_cTMZ, t_cTMZ] = scheduleCTMZ(RTdays);
        StartAdj = RTdays(end) + 3;
        [d_aTMZ, t_aTMZ] = scheduleATMZ_from_start(StartAdj);
        
        evt_time = [d_CART, d_cTMZ, d_aTMZ];
        evt_type = [t_CART, t_cTMZ, t_aTMZ];
        evt_dose = [dose_CART, zeros(size(d_cTMZ)), zeros(size(d_aTMZ))];
        
    elseif prot == "StuppC-3C-StuppA-3C"
        M_total = 6;
        dose_m = (D_tot / M_total) / K; % 50-50 Split between blocks implicitly
        
        RTdays = make_base_RT(Lr, Lr1, Tr);
        [d_cTMZ, t_cTMZ] = scheduleCTMZ(RTdays);
        TendCTMZ = max(d_cTMZ);
        
        d_CART1 = TendCTMZ + TgapSC + (0:2)*TgapC;
        TendCART1 = d_CART1(end);
        
        StartAdj = TendCART1 + TgapSC;
        [d_aTMZ, t_aTMZ] = scheduleATMZ_from_start(StartAdj);
        TendAdj = max(d_aTMZ);
        
        d_CART2 = TendAdj + TgapSC + (0:2)*TgapC;
        
        d_CART = [d_CART1, d_CART2];
        t_CART = repmat("CART", 1, M_total);
        dose_CART = repmat(dose_m, 1, M_total);
        
        evt_time = [d_cTMZ, d_aTMZ, d_CART];
        evt_type = [t_cTMZ, t_aTMZ, t_CART];
        evt_dose = [zeros(size(d_cTMZ)), zeros(size(d_aTMZ)), dose_CART];
    end
    
    % Sort events chronologically
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
    Tmax = RTdays(end);
    days = 1:Tmax;
    type = strings(size(days));
    type(:) = "TMZ_RT";
    type(ismember(days, RTdays)) = "RT_TMZ";
end

% Schedules Adjuvant TMZ based on start day directly
function [days, type] = scheduleATMZ_from_start(StartDay)
    days = [];
    for cycle = 0:5
        cycleStart = StartDay + 28*cycle;
        days = [days, cycleStart + (0:4)];
    end
    type = repmat("TMZ_ADJ", size(days));
end

% Simulates full treatment cycle for a single patient
% Simulates full treatment cycle for a single patient
function [te, minT] = calc_1cycle(ics, evt_time, evt_type, evt_dose, ...
                                  r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, ...
                                  epsilon1, gamma1, gamma2, gamma3, mu, K, beta1, beta2, ...
                                  tau, A, alphaT, B, gamma, D, Np, E0, Tfinal)
    M_idx = 0; t0 = 0;
    minT = inf;
    for i = 1:length(evt_time)
        % Only run the ODE solver if time has advanced
        if evt_time(i) > t0
            [t, y, te] = calc_event(ics, t0, evt_time(i), r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np, K);
            minT = min(minT, min(sum(y(:,1:6), 2)));
            
            if ~isempty(te) 
                return; 
            end
            ics = y(end,:); % Update state vector for the treatment application
        end
        
        % Apply the treatment instantly
        ics = applyTreatment(ics, evt_type(i), evt_dose(i), A, alphaT, B, D, gamma, E0);
        M_idx = M_idx + 1;
        t0 = evt_time(i);
    end
    
    % Final integration after the last treatment event
    if Tfinal > t0
        [t, y, te] = calc_event(ics, t0, Tfinal, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np, K);
        minT = min(minT, min(sum(y(:,1:6), 2)));
    end
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

    [t, y, te, ye] = ode78(@(t, y) mODE(t, y, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu), tspan, ics, options);
end

% System of Ordinary Differential Equations
function dy = mODE(~, y, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu)
    dy = zeros(9, 1);
    dy(1) = r1*y(1)*(1 - y(1) - y(2) - y(3) - y(4) - y(5) - y(6) - y(7)) - beta1*y(1) + beta2*y(2) - (alpha1 + epsilon1)*y(9)*y(1) - a2*y(8)*y(1);
    dy(2) = beta1*y(1) - beta2*y(2);
    dy(3) = r1*y(3)*(1 - y(1) - y(2) - y(3) - y(4) - y(5) - y(6) - y(7)) - beta1*y(3) + beta2*y(4) - (alpha1 + epsilon1)*y(9)*y(3);
    dy(4) = beta1*y(3) - beta2*y(4);
    dy(5) = r2*y(5)*(1 - y(1) - y(2) - y(3) - y(4) - y(5) - y(6) - y(7)) - beta1*y(5) + beta2*y(6) - a2*y(8)*y(5) + epsilon1*(y(1) + y(3))*y(9);
    dy(6) = beta1*y(5) - beta2*y(6);
    dy(7) = alpha1*(y(1) + y(3))*y(9) + a2*(y(1) + y(5))*y(8) - tau*y(7);
    dy(8) = -rho1*y(8) + (rho2.*y(1)*y(8))./(gamma1 + y(1)) + (rho3*y(5)*y(8))/(gamma2 + y(5)) - rho4*(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7))*y(8)/(gamma3+y(8)) - alpha3*y(9)*y(8);
    dy(9) = -mu*y(9);
end

% Applies Radiation Therapy (RT) effects (Includes alphaT for CAR-T survival)
function y = applyR(y, A, alphaT, B, D, gamma)
    SF = exp(-A*D - B*D^2);
    SFC = exp(-alphaT*D);
    y(7) = y(7) + (1-SF)*y(1) + (1-SF)*y(3)+ (1-SF)*y(5);
    y(1) = SF*y(1) + gamma*y(2);
    y(2) = (1-gamma)*y(2);
    y(3) = SF*y(3) + gamma*y(4);
    y(4) = (1-gamma)*y(4);
    y(5) = SF*y(5) + gamma*y(6);
    y(6) = (1-gamma)*y(6);
    y(8) = SFC*y(8);
end

% Evaluates treatment effects at discrete time points
function y = applyTreatment(y, treatment, dose, A, alphaT, B, D, gamma, E0)
    switch treatment
        case 'RT_TMZ'
            y = applyR(y, A, alphaT, B, D, gamma);
            y(9) = y(9) + E0/3;
        case "TMZ_RT"
            y(9) = y(9) + E0/3;
        case "TMZ_ADJ"
            y(9) = y(9) + 2*E0/3;
        case "CART"
            y(8) = y(8) + dose;
    end
end