function [J, sinr_vals] = calc_objective_fixed_r12(r3_list, r1_fixed, r2_fixed, H_data, W, param)
% 目的 (Purpose):
%   在 "Pointing Only" 优化阶段，固定极化轴 (r1, r2)，仅优化视轴 r3。
%   计算系统的总代价 J (SoftMin SINR + Penalties)。
% 输入 (Inputs):
%   r3_list   - [3 x M]   当前的视轴向量 (优化变量)。通常定义在 Sphere 流形上。
%   r1_fixed  - [3 x M]   固定的本地 X 轴 (H极化方向)。
%   r2_fixed  - [3 x M]   固定的本地 Y 轴 (V极化方向)。
%   H_data    - Struct
%               预计算的信道数据。
%   W         - [M x K]   发射波束赋形矩阵。
%   param     - Struct
%               仿真参数。
% 输出 (Outputs):
%   J         - Double
%               总目标函数值。
%   sinr_vals - [K x 1]
%               用户 SINR 值。
    [~, M] = size(r3_list);
    
    R = zeros(3, 3, M);
    for m = 1:M
        R(:, 1, m) = r1_fixed(:, m);
        R(:, 2, m) = r2_fixed(:, m);
        R(:, 3, m) = r3_list(:, m); 
    end
    

    [~, K] = size(W);
    [H_val, ~, ~] = compute_H_internal(R, H_data, param);
    

    S_vals = zeros(K,1); I_vals = zeros(K,1); sinr_vals = zeros(K,1);
    for k=1:K
        S_vals(k) = abs(W(:,k)' * H_val(:,k))^2;
        intf = 0;
        for j=1:K, if j~=k, intf = intf + abs(W(:,j)' * H_val(:,k))^2; end, end
        I_vals(k) = intf + param.sigma2;
        sinr_vals(k) = S_vals(k) / I_vals(k);
    end
    
    % 1. Main Objective: Soft-min (-t_mu)
    % Maximize Min SINR => Minimize -(-1/mu * log(sum(exp(-mu*gamma))))
    % => Minimize (1/mu) * log(...)
    mu = param.mu;
    % Log-Sum-Exp trick for stability
    max_exp = -mu * max(sinr_vals); 
    % Note: We want to maximize t_mu, so we minimize -t_mu.
    % t_mu = -1/mu * log(sum(exp(-mu * gamma)))
    % Obj = - t_mu = 1/mu * log(sum(exp(-mu * gamma)))
    term_main = (1/mu) * log(sum(exp(-mu * sinr_vals)));
    
    % 2. Threshold Penalty: Sum_k ( softplus(gamma_th - gamma_k) )^2
    lam_sinr = param.lambda_sinr;
    gamma_th = param.gamma_th;
    alpha = param.alpha;
    pen_sinr = 0;
    for k=1:K
        delta = gamma_th - sinr_vals(k);
        if delta > 50/alpha
            sp = delta;
        else
            sp = (1/alpha)*log(1+exp(alpha*delta));
        end
        pen_sinr = pen_sinr + lam_sinr * (sp^2);
    end
    
    % 3. Angle Penalty
    lam_ang = param.lambda_angle;
    cos_max = cos(param.theta_max);
    z_hat = [0;0;1];
    pen_ang = 0;
    M = size(R,3);
    for m=1:M
        r3 = R(:,3,m);
        delta = cos_max - r3'*z_hat;
        if delta > 50/alpha
            sp = delta;
        else
            sp = (1/alpha)*log(1+exp(alpha*delta));
        end
        pen_ang = pen_ang + lam_ang * (sp^2);
    end
    
    J = term_main + pen_sinr + pen_ang;
end
