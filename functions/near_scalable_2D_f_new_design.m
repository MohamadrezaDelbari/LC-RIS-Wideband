function [S,s] = near_scalable_2D_f_new_design(Paramphase,ParamC)

% --- Preserve exact variable names ---
Power=Paramphase.Power;
kappa=ParamC(1).kappa;
kappa2=Paramphase.kappa; 
gamma=Paramphase.gamma;
N=length(Paramphase.pris);
% s_opt=Paramphase.s_opt.';
Imax=120;
minSR=-100;
mu=0.1;

%______________________Pre-Calculation (Vectorized)________________________
% 1. Setup User Locations
prx=Paramphase.pp_mu(:,1);
pry=Paramphase.pp_mu(:,2);
ppr_x=linspace(prx(1),prx(2),(prx(2)-prx(1))/Paramphase.step+1);
ppr_y=linspace(pry(1),pry(2),(pry(2)-pry(1))/Paramphase.step+1);
[XX_u, YY_u] = meshgrid(ppr_x, ppr_y);
pos_u_all = [XX_u(:), YY_u(:), repmat(Paramphase.pp_mu(1,3), numel(XX_u), 1)];

% 2. Setup Eve Locations
prxe=Paramphase.pp_e(:,1);
prye=Paramphase.pp_e(:,2);
ppr_xe=linspace(prxe(1),prxe(2),(prxe(2)-prxe(1))/Paramphase.step+1);
ppr_ye=linspace(prye(1),prye(2),(prye(2)-prye(1))/Paramphase.step+1);
[XX_e, YY_e] = meshgrid(ppr_xe, ppr_ye);
pos_e_all = [XX_e(:), YY_e(:), repmat(Paramphase.pp_e(1,3), numel(XX_e), 1)];

% Dimensions
NumU = size(pos_u_all, 1);
NumE = size(pos_e_all, 1);
NumK = length(kappa2);
TotalScenarios = NumU * NumE * NumK;

% Pre-allocate flat arrays
a_user_flat = zeros(N, TotalScenarios);
a_eve_flat  = zeros(N, TotalScenarios);
factor_flat = zeros(1, TotalScenarios);

idx = 1;

% --- Channel Generation Loop (Fixed for Struct Arrays) ---
for ie = 1:NumE
    p_e_curr = pos_e_all(ie, :);
    ParamC_e = ParamC; 
    
    % FIX 1: Loop to assign p_mu to struct array
    for k=1:numel(ParamC_e)
        ParamC_e(k).p_mu = p_e_curr;
    end
    
    for iu = 1:NumU
        p_u_curr = pos_u_all(iu, :);
        ParamC_u = ParamC;
        
        % FIX 2: Loop to assign p_mu to struct array
        for k=1:numel(ParamC_u)
            ParamC_u(k).p_mu = p_u_curr;
        end
        
        for t = 1:NumK
            kappa_temp = kappa2(t);
            f_val = 1+2.4*(kappa_temp-kappa)/kappa;
            
            % Update kappa for User Struct Array
            for k=1:numel(ParamC_u)
                ParamC_u(k).kappa = kappa_temp;
            end
            [~,H_i_u,H_r_u,Po_u] = func_channel(ParamC_u);
            vec_u = diag(H_r_u{1})*H_i_u*Po_u(1).a_bs_irs(:,1)*sqrt(Power);
            
            % Update kappa for Eve Struct Array
            for k=1:numel(ParamC_e)
                ParamC_e(k).kappa = kappa_temp;
            end
            [~,H_i_e,H_r_e,Po_e] = func_channel(ParamC_e);
            vec_e = diag(H_r_e{1})*H_i_e*Po_e(1).a_bs_irs(:,1)*sqrt(Power);
            
            % Store
            a_user_flat(:, idx) = vec_u;
            a_eve_flat(:, idx)  = vec_e;
            factor_flat(idx)    = f_val;
            idx = idx + 1;
        end
    end
end

% Pre-compute dot products (Analytical inputs)
uu = sum(abs(a_user_flat).^2, 1);             
ee = sum(abs(a_eve_flat).^2, 1);              
ue = sum(conj(a_user_flat) .* a_eve_flat, 1); 
eu = conj(ue);

for initial_index = 1:10
% --- Optimization Setup --- 
s_0 = exp(1i*2*pi*rand(N,1));
% s_0 = exp(1i*mod(angle(mean(a_user_flat,2)),2*pi));
% s_0 = s_opt;

% --- Main Optimization Loop ---
for i_gamma=1:2
    
    % TRICK 1: Analytical Max Eigenvalue (Rank-2)
    T_val = ee - gamma * uu;
    D_val = -gamma * (ee .* uu - ue .* eu); 
    
    delta = sqrt(T_val.^2 - 4 * D_val);
    lambda_max_vec = real(0.5 * (T_val + delta));
    lambda_max_vec = lambda_max_vec + 1e-4; % PSD Safety
    
    for i=1:Imax
        s_pow = s_0 .^ factor_flat;
        % s_pow_opt = s_opt.^ factor_flat;
        
        % TRICK 2: Implicit Matrix Multiplication
        adot_u = sum(a_user_flat .* s_pow, 1); 
        adot_e = sum(a_eve_flat .* s_pow, 1);  
        % adot_u_opt = sum(a_user_flat .* s_pow_opt, 1); 
        % adot_e_opt = sum(a_eve_flat .* s_pow_opt, 1); 
        
        gain_u = abs(adot_u).^2;
        gain_e = abs(adot_e).^2;
        % gain_u_opt = abs(adot_u_opt).^2;
        % gain_e_opt = abs(adot_e_opt).^2;
        SR_vec = log2(1 + gain_u) - log2(1 + gain_e);
        % min_SR_vec_opt = min(log2(1 + gain_u_opt) - log2(1 + gain_e_opt));
        
        % Beta Calculation
        T1 = lambda_max_vec .* s_pow; 
        T2 = conj(a_eve_flat) .* adot_e - gamma * (conj(a_user_flat) .* adot_u);
        beta_all = T1 - T2;
        
        % LSE Aggregation
        min_sr_val = min(SR_vec);
        w_numer = exp(-(SR_vec - min_sr_val)/mu);
        w_denom = sum(w_numer);
        weights = w_numer / w_denom;
        
        LSE = sum(beta_all .* weights, 2);
        
        s_0 = exp(1i * mod(angle(LSE), 2*pi));
        
        current_min = min(SR_vec);
        if current_min > minSR
            minSR = current_min;
            s_final = s_0;
        end
    end
    
    gamma = 1/(2^minSR);
    fprintf('Gamma Loop %d, Best MinSR: %.4f\n', i_gamma, minSR);
end
end
s = s_final;
S = diag(s);
end