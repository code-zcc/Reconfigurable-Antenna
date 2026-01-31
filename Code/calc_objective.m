function [J, sinr_vals] = calc_objective(R, H_data, W, param)
% 功能:
%   计算在给定天线旋转矩阵 R 的情况下的系统总代价 J。
% 输入 (Inputs):
%   R       - [3 x 3 x M] 当前的旋转矩阵 (优化变量)。
%   H_data  - Struct
%             预计算的信道几何数据 
%   W       - [M x K]  当前固定的发射波束赋形矩阵。
%   param   - Struct
%             仿真参数 
% 输出 (Outputs):
%   J       -  总目标函数值 (Minimize this)。
%   sinr_vals - [K x 1] 当前所有用户的实际 SINR 线性值。
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
