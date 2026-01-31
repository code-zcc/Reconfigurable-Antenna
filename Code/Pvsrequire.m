clear; clc; close all;
load('ChannelData_100.mat');
sim_param.fc = 2.4e9;            
sim_param.c  = 3e8;               
sim_param.lambda = sim_param.c / sim_param.fc;
sim_param.d_ant = sim_param.lambda / 2;


sim_param.M_x = 4;               
sim_param.M_y = 4;               
sim_param.M   = sim_param.M_x * sim_param.M_y; 
sim_param.K   = 5;               

sim_param.L = 8;                 
sim_param.theta_max = pi/5;
sim_param.noise_pwr_dBm = -80;   
sim_param.sigma2 = 10^((sim_param.noise_pwr_dBm)/10);
sim_param.p = 2;          
sim_param.G0 = 2*(2*sim_param.p + 1);
sim_param.XPD = 0.1;   
sim_param.K_factor = 10^(10/10); 

sim_param.mu = 10;
sim_param.lambda_sinr = 50;
sim_param.lambda_angle = 100;
sim_param.alpha = 20;

bs_pos = zeros(3, sim_param.M);
idx = 1;
for nx = 0 : sim_param.M_x - 1
    for ny = 0 : sim_param.M_y - 1

        pos_x = (nx - (sim_param.M_x-1)/2) * sim_param.d_ant;
        pos_y = (ny - (sim_param.M_y-1)/2) * sim_param.d_ant;
        bs_pos(:, idx) = [pos_x; pos_y; 0];
        idx = idx + 1;
    end
end

user_pos = zeros(3, sim_param.K);


r_min = 30; r_max = 60;
theta_bound = deg2rad(70);


bs_rotations = zeros(3, sim_param.M);
ue_rotations = zeros(3, sim_param.K);


TargetRates = 1:5; 
NumSamples = 100; 
itermax = 100; epsilon_stop = 1e-3;

P_Sch1 = zeros(length(TargetRates), NumSamples,itermax);
P_Sch2 = zeros(length(TargetRates), NumSamples,itermax);
P_Sch3 = zeros(length(TargetRates), NumSamples,itermax);
P_Sch5 = zeros(length(TargetRates), NumSamples,itermax);
P_Sch4 = zeros(length(TargetRates), NumSamples,1);%只做digital
for r_idx = 1:length(TargetRates)
    
    current_rate = TargetRates(r_idx);
    sim_param.gamma_th = 2^(current_rate) - 1;
    gamma_vec = sim_param.gamma_th * ones(sim_param.K, 1);
    P_Sch1_slice = zeros(NumSamples, itermax);
    P_Sch2_slice = zeros(NumSamples, itermax);
    P_Sch3_slice = zeros(NumSamples, itermax);
    P_Sch5_slice = zeros(NumSamples, itermax);
    P_Sch4_slice = zeros(NumSamples, 1);
    user_pos=zeros(3,sim_param.K);scatter_pos=zeros(3,sim_param.L);
    nlos_phase_bank=zeros(sim_param.L,sim_param.M,sim_param.K,2);
    parfor s = 1:NumSamples
        param_local = sim_param;
        user_pos = ChannelSamples(s).user_pos;
        scatter_pos = ChannelSamples(s).scatter_pos;
        nlos_phase_bank = ChannelSamples(s).nlos_phase_bank;


        F_analog_base = ChannelSamples(s).F_init;
        U_user_base   = ChannelSamples(s).U_init;
        

        bs_rot_base = zeros(3, sim_param.M); 
        ue_rot_base = zeros(3, sim_param.K);
        R_init_base = zeros(3, 3, sim_param.M);
        for m=1:sim_param.M
            R_init_base(:,:,m)=eye(3);
        end
        H_eff_start = construct_channel_H(bs_pos, user_pos, scatter_pos, bs_rot_base, ue_rot_base,  F_analog_base, U_user_base, sim_param,nlos_phase_bank);
        H_func_data_start = prepare_RCG_data_multipath(bs_pos, user_pos, scatter_pos, ...
            ue_rot_base, sim_param, ...
            nlos_phase_bank, U_user_base,F_analog_base);
        [~, V_current_base] = initial_ZF(H_eff_start, sim_param.sigma2, gamma_vec);
        [W_opt_base, V_opt_base, ~] = optimize_W_Scheme1_ZF(H_eff_start, V_current_base, gamma_vec, sim_param);
        

        P_init_val = trace(sum(V_opt_base, 3));
        P_Sch4_slice(s) = P_init_val;
        PP  = zeros(1, itermax);
        PP1 = zeros(1, itermax);
        PP2 = zeros(1, itermax);
        %% 方案5
        PP5 = zeros(1, itermax);
        PP5(1) = P_init_val;
        W_opt = W_opt_base;
        V_opt = V_opt_base;
        R_curr = R_init_base; % 当前的 R
        r1_fixed = R_curr(:, 1, :); 
        r2_fixed = R_curr(:, 2, :);
        r3_curr  = R_curr(:, 3, :);
        r3_curr  = squeeze(r3_curr);
    r1_fixed = squeeze(r1_fixed);
    r2_fixed = squeeze(r2_fixed);
        F_curr = F_analog_base;
        U_curr = U_user_base;
        for i=1:itermax-1
        H_func_data_S5 = prepare_RCG_data_multipath(bs_pos, user_pos, scatter_pos, ue_rot_base, ...
            sim_param, nlos_phase_bank, U_curr, F_curr);
        if sim_param.theta_max<1e-4
                R_opt=R_init_base;
            else
                for t=1:30
                    [R_opt1, min_sinr_val1]=run_RCG_step_pointing_only(r3_curr,r1_fixed,r2_fixed, W_opt, H_func_data_S5, param_local);
                    
                    z_hat=[0;0;1];
                    max_violation = 0;
                    for m=1:sim_param.M
                        r3 = R_opt1(:, 3, m);
                        proj = r3.' * z_hat;
                        if proj > 1, proj = 1; end
                        if proj < -1, proj = -1; end
                        theta_curr = acos(proj);
                        if theta_curr > sim_param.theta_max
                            viol = theta_curr - sim_param.theta_max;
                            if viol > max_violation
                                max_violation = viol;
                            end
                        end
                    end
                    sinr_target = param_local.gamma_th * 0.999; 
                    is_sinr_bad = (min_sinr_val1 <= sinr_target);
                    if (max_violation < 1e-4) && (~is_sinr_bad) 
                        R_opt=R_opt1;
                        param_local.lambda_angle =100;
                        param_local.lambda_sinr  = 50;
                        break;
                    else
                        if max_violation >= 1e-4
                            param_local.lambda_angle = param_local.lambda_angle * 10;
                        end
                        if is_sinr_bad
                            param_local.lambda_sinr = param_local.lambda_sinr * 2;
                        end

                        if (param_local.lambda_angle > 1e10) || (param_local.lambda_sinr > 1e10)
                            R_opt=R_opt;
                            param_local.lambda_angle=100;
                            param_local.lambda_sinr  = 50;
                            break;
                        end
                    end
                end
            end
            [H_val, ~, ~] = compute_H_internal(R_opt, H_func_data_S5, sim_param);
            [W_opt, V_opt, status] = optimize_W_Scheme1_ZF(H_val, V_opt, gamma_vec, sim_param);
            PP5(1+i)=trace(sum(V_opt,3));

            if (abs(PP5(i) - PP5(1+i)) / PP5(i)) < epsilon_stop
                PP5(1+i:end) = PP5(1+i);
                break;
            end
        end
        P_Sch5_slice(s, :) = PP5;
        %% 方案1
        U_new = U_user_base; 
        F_new = F_analog_base; 
        R_opt = R_init_base;
        W_opt = W_opt_base; 
        V_opt = V_opt_base;
        PP(1) = P_init_val;
        for i=1:itermax-1
            for loop=1:30
                        H_func_data = prepare_RCG_data_multipath(bs_pos, user_pos, scatter_pos, ue_rot_base, ...
                sim_param, nlos_phase_bank, U_new, F_new);
            if sim_param.theta_max<1e-4
                R_opt=R_init_base;
            else
                for t=1:30
                    [R_opt1, min_sinr_val1]=run_RCG_step(R_opt, W_opt, H_func_data, param_local);
                    z_hat=[0;0;1];
                    max_violation = 0;
                    for m=1:sim_param.M
                        r3 = R_opt1(:, 3, m);
                        proj = r3.' * z_hat;
                        if proj > 1, proj = 1; end
                        if proj < -1, proj = -1; end
                        theta_curr = acos(proj);
                        if theta_curr > sim_param.theta_max
                            viol = theta_curr - sim_param.theta_max;
                            if viol > max_violation
                                max_violation = viol;
                            end
                        end
                    end
                    sinr_target = param_local.gamma_th * 0.999; 
                    is_sinr_bad = (min_sinr_val1 <= sinr_target);
                    if (max_violation < 1e-4) && (~is_sinr_bad)
                        R_opt=R_opt1;
                        param_local.lambda_angle =100;
                        param_local.lambda_sinr  = 50;
                        break;
                    else
                        if max_violation >= 1e-4
                            param_local.lambda_angle = param_local.lambda_angle * 10;
                        end
                        if is_sinr_bad
                            param_local.lambda_sinr = param_local.lambda_sinr * 2;
                        end

                        if (param_local.lambda_angle > 1e10) || (param_local.lambda_sinr > 1e10)
                            R_opt=R_opt;
                            param_local.lambda_angle=100;
                            param_local.lambda_sinr  = 50;
                            break;
                        end
                    end
                end
            end

            for loop=1:30
            [F_candidate, min_sinr_val2] = run_F_step(F_new, R_opt, W_opt, H_func_data, sim_param);
            if min_sinr_val2 <= 0.999 * param_local.gamma_th
                    param_local.lambda_sinr = param_local.lambda_sinr * 2;
                    if param_local.lambda_sinr > 1e10
                        F_new = F_new;
                        param_local.lambda_sinr =50;
                        break;
                    end
            else
                    F_new = F_candidate;
                    param_local.lambda_sinr =50;
                    break;
                end
            end
            H_func_data.f_analog = F_new;
            H_raw_UE = compute_raw_channel_for_U(bs_pos, user_pos, scatter_pos, R_opt, ...
                ue_rot_base, F_new, sim_param, nlos_phase_bank);
                
            [U_candidate, min_sinr_U] = run_U_step(U_new, H_raw_UE, W_opt, sim_param);
            if min_sinr_U <= 0.999 * param_local.gamma_th

                    param_local.lambda_sinr = param_local.lambda_sinr * 2;

                    if param_local.lambda_sinr > 1e10
                        U_new = U_new;
                        param_local.lambda_sinr =50;
                        break;
                    end
                else

                    U_new = U_candidate;
                    param_local.lambda_sinr =50;
                    break;
                end
            end
            
            H_func_data = prepare_RCG_data_multipath(bs_pos, user_pos, scatter_pos, ue_rot_base, ...
                sim_param, nlos_phase_bank, U_new, F_new);
            [H_val, ~, ~] = compute_H_internal(R_opt, H_func_data, sim_param);
            [W_opt, V_opt, status] = optimize_W_Scheme1_ZF(H_val, V_opt, gamma_vec, sim_param);
            PP(1+i)=trace(sum(V_opt,3));

            if (abs(PP(i) - PP(1+i)) / PP(i)) < epsilon_stop
                PP(1+i:end) = PP(1+i);
                break;
            end
        end
        P_Sch1_slice(s, :) = PP;
        %% 方案2
        param_local = sim_param;
        R_opt = R_init_base;
        W_opt = W_opt_base; 
        V_opt = V_opt_base;       
        PP1(1) = P_init_val;
        H_func_data = H_func_data_start;
        for i=1:itermax-1
                for t=1:40
                    [R_opt1, min_sinr_val1]=run_RCG_step(R_opt, W_opt, H_func_data, param_local);
                    z_hat=[0;0;1];
                    max_violation = 0;
                    for m=1:sim_param.M
                        r3 = R_opt1(:, 3, m);
                        proj = r3.' * z_hat;
                        if proj > 1, proj = 1; end
                        if proj < -1, proj = -1; end
                        theta_curr = acos(proj);
                        if theta_curr > sim_param.theta_max
                            viol = theta_curr - sim_param.theta_max;
                            if viol > max_violation
                                max_violation = viol;
                            end
                        end
                    end
                    sinr_target = param_local.gamma_th * 0.999; 
                    is_sinr_bad = (min_sinr_val1 <= sinr_target);
                    if (max_violation < 1e-5) && (~is_sinr_bad) 

                        R_opt=R_opt1;
                        param_local.lambda_angle =100;
                        param_local.lambda_sinr  = 50;
                        break;
                    else
                        if max_violation >= 1e-5
                            param_local.lambda_angle = param_local.lambda_angle * 10;
                        end
                        if is_sinr_bad
                            param_local.lambda_sinr = param_local.lambda_sinr * 2;
                        end

                        if (param_local.lambda_angle > 1e10) || (param_local.lambda_sinr > 1e10)
                            R_opt=R_opt;
                            param_local.lambda_angle=100;
                            param_local.lambda_sinr  = 50;
                            break;
                        end
                    end
                end

            [H_val, ~, ~] = compute_H_internal(R_opt, H_func_data, sim_param);
            [W_opt, V_opt, status] = optimize_W_Scheme1_ZF(H_val, V_opt, gamma_vec, sim_param);
            PP1(1+i)=trace(sum(V_opt,3));

            if (abs(PP1(i) - PP1(1+i)) / PP1(i)) < epsilon_stop
                PP1(1+i:end) = PP1(1+i);
                break;
            end
        end
        P_Sch2_slice(s, :) = PP1;
        %% 方案3
        U_new = U_user_base; 
        F_new = F_analog_base; 
        W_opt = W_opt_base; 
        V_opt = V_opt_base;

        R_opt = R_init_base; 
        
        PP2(1) = P_init_val;
        for i=1:itermax-1
            H_raw_UE = compute_raw_channel_for_U(bs_pos, user_pos, scatter_pos, R_opt, ...
                ue_rot_base, F_new, sim_param, nlos_phase_bank);

            for loop=1:30
            [U_candidate, min_sinr_U] = run_U_step(U_new, H_raw_UE, W_opt, sim_param);
            if min_sinr_U <= 0.999 * param_local.gamma_th

                    param_local.lambda_sinr = param_local.lambda_sinr * 2;
                    

                    if param_local.lambda_sinr > 1e10
                        U_new = U_new;
                        param_local.lambda_sinr =50;
                        break;
                    end
                else

                    U_new = U_candidate;
                    param_local.lambda_sinr =50;
                    break;
                end
            end
            H_func_data = prepare_RCG_data_multipath(bs_pos, user_pos, scatter_pos, ue_rot_base, ...
                sim_param, nlos_phase_bank, U_new, F_new);
            for loop=1:30
            [F_candidate, min_sinr_val2] = run_F_step(F_new, R_opt, W_opt, H_func_data, sim_param);
            if min_sinr_val2 <= 0.999 * param_local.gamma_th

                    param_local.lambda_sinr = param_local.lambda_sinr * 2;
                    
                  
                    if param_local.lambda_sinr > 1e10
                        F_new = F_new;
                        param_local.lambda_sinr =50;
                        break;
                    end
            else
                    F_new = F_candidate;
                    param_local.lambda_sinr =50;
                    break;
                end
            end
            H_func_data.f_analog = F_new;
            [H_val, ~, ~] = compute_H_internal(R_opt, H_func_data, sim_param);
            [W_opt, V_opt, status] = optimize_W_Scheme1_ZF(H_val, V_opt, gamma_vec, sim_param);
            PP2(1+i)=trace(sum(V_opt,3));

            if (abs(PP2(i) - PP2(1+i)) / PP2(i)) < epsilon_stop
                PP2(1+i:end) = PP2(1+i);
                break;
            end
        end
        P_Sch3_slice(s, :) = PP2;


    end
    P_Sch1(r_idx, :, :) = P_Sch1_slice;
    P_Sch2(r_idx, :, :) = P_Sch2_slice;
    P_Sch3(r_idx, :, :) = P_Sch3_slice;
    P_Sch4(r_idx, :, :) = P_Sch4_slice;
    P_Sch5(r_idx, :, :) = P_Sch5_slice;
end
save('Pvsrequire_long1.mat', 'P_Sch1', 'P_Sch2', 'P_Sch3', 'P_Sch4','P_Sch5', 'TargetRates');