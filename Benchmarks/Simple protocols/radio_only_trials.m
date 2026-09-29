clear variables;
close all;

% 1. Fixed Parameters and Constants
g1 = 10^10; g2 = 10^10; g3 = 2*10^9; 
v = 0.5*10^9; 
K = 5*10^12; 
mu = 8.32; E0 = 0; D = 2;
V = v/K; 
alpha2 = 2.5*10^(-10); a2 = alpha2*K;  
gamma1 = g1/K; gamma2 = g2/K; gamma3 = g3/K;

% Simulation settings
Tfinal = 150000; % Extended time for survival tracking

% RT cycle configuration: time interval in days between RT applications. 
% L11+1 represents the total number of RT applications per cycle. 
Tr = 7; Lr = 6; Tr1 = 1; Lr1 = 5; 

% ODE integration settings
Np = 10; 
Sy = (Lr*(Lr1))*Np + Lr*Np + Np; 

% 2. Load Virtual Patient Cohort (10,000 patients)
load('Params_defenitve_-2026-8-25-17-35_N=10000.mat');

% Uncomment to initialize parallel pool if not already running
% delete(gcp('nocreate'));
% pool = parpool("Processes", 10);

% Preallocate array for parallel computing
Tsurv = zeros(1, N);

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
    gamma = gammaval(n); A = Aval(n); B = Bval(n);
    
    % Compute the solution of the ODE system corresponding to the RT Only protocol
    [te] = calc_1cycle(r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,V,K,beta1,beta2,tau,A,B,gamma,D,Np, ...
                       Tr,Lr,Tr1,Lr1,S10,S20,Q10,Q20,S30,Q30,E0,Tfinal);
    
    % Store survival time
    Tsurv(n) = te;
end
toc

% 4. Post-processing and Statistics
% Save cohort survival results
save("Radio_only_Tsurv-" + year(datetime) + "-" + month(datetime) + "-" + day(datetime) + "-" + hour(datetime) + "-" + minute(datetime) + "_N=" + N + ".mat");

% Display Median Overall Survival (Months)
disp('Median OS (months):');
disp(median(Tsurv)/30);

% Calculate and display 95% Confidence Intervals via Bootstrapping
nBoot = 10000;
bootstat = bootstrp(nBoot, @median, Tsurv);
ci = prctile(bootstat, [2.5 97.5]);
disp('95% CI (months):');
disp(ci./30);

% =========================================================================
% FUNCTIONS
% =========================================================================

% ODE solver wrapper with event detection (Lethal tumor burden threshold)
function [t, y, te, ye] = calc_event(ics, T1, T2, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np, K)
    
    tspan = linspace(T1, T2, Np);
    yy0 = ics;
    
    options = odeset('RelTol', 1e-10, 'Events', @events, 'NonNegative', [1,2,3,4,5,6,7,8,9]); 
    
    function [position, isterminal, direction] = events(~, y)
        % Integration terminates when tumor burden reaches the lethal threshold (K/5)
        position = K.*(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)) - K/5;                                           
        isterminal = 1;                                                        
        direction = 0;                                                         
    end

    [t, y, te, ye] = ode45(@(t, y) mODE(t, y, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu), tspan, yy0, options);
end

% Simulates the full treatment sequence (RT Only)
function [te] = calc_1cycle(r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, mu, V, K, beta1, beta2, tau, A, B, gamma, D, Np, Tr, Lr, Tr1, Lr1, S10, S20, Q10, Q20, S30, Q30, E0, Tfinal)
    
    ics(1) = S10/K; ics(2) = Q10/K; ics(3) = S20/K; ics(4) = Q20/K; ics(5) = S30/K; ics(6) = Q30/K; ics(7) = 0; ics(8) = 0; ics(9) = E0;
    
    % Generate treatment schedules
    RTdays = makeRTschedule(Lr, Lr1, Tr);
    events.time = RTdays;
    events.type = repmat({'RT'}, size(RTdays));
    
    Sy = (length(events.time) + 1) * Np;
    Y = zeros(Sy, 9); 
    tt = zeros(Sy, 1);
    
    M = 0; t0 = 0;
    
    % Piecewise integration between treatment events
    for i = 1:length(events.time)
        [t, y, te] = calc_event(ics, t0, events.time(i), r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np, K);
        
        if ~isempty(te); return; end
        
        ics = applyTreatment(y(end,:), events.type{i}, A, B, D, gamma);
        tt(M*Np+1 : M*Np+Np) = t;
        Y(M*Np+1 : M*Np+Np, :) = y;
        M = M + 1;
        t0 = events.time(i);
    end
    
    % Final integration step until Tfinal or lethal event
    [t, y, te] = calc_event(ics, t0, Tfinal, r1, r2, alpha1, a2, alpha3, rho1, rho2, rho3, rho4, epsilon1, gamma1, gamma2, gamma3, beta1, beta2, tau, mu, Np, K);
    
    if ~isempty(te); return; end
    
    tt(M*Np+1 : M*Np+length(t)) = t;
    Y(M*Np+1 : M*Np+length(t), :) = y;
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

% Evaluates treatment effects at discrete time points
function y = applyTreatment(y, treatment, A, B, D, gamma)
    switch treatment
        case 'RT'
            y = applyR(y, A, B, D, gamma);
    end
end