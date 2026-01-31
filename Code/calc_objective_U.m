function [J, sinr_vals] = calc_objective_U(U, H_raw, W, param)
% 目的 (Purpose):
%   计算在给定接收波束赋形矩阵 U 的情况下的系统总代价 J。
% 输入 (Inputs):
%   U      - [2 x K] 
%            当前的接收波束赋形矩阵 (优化变量)。
%   H_raw  - [2 x M x K]
%            原始信道张量。
%   W      - [M x K]
%            发射波束赋形矩阵。
%   param  - Struct
%            仿真参数。
% 输出 (Outputs):
%   J     总目标函数值。
%   sinr_vals - [K x 1] 当前所有用户的实际 SINR 线性值。
    [~, K] = size(U);
    sigma2 = param.sigma2;
    
    V_arr = zeros(2, K, K); 
    for k = 1:K
        H_k = H_raw(:,:,k); 
        V_arr(:, :, k) = H_k * conj(W);
    end
    
    S_vals = zeros(K,1); I_vals = zeros(K,1); sinr_vals = zeros(K,1);
    for k = 1:K
        u_k = U(:, k);
        v_sig = V_arr(:, k, k);
        S_vals(k) = abs(u_k' * v_sig)^2;
        intf = 0;
        for j = 1:K
            if j ~= k
                v_int = V_arr(:, j, k);
                intf = intf + abs(u_k' * v_int)^2;
            end
        end
        I_vals(k) = intf + sigma2;
        sinr_vals(k) = S_vals(k) / I_vals(k);
    end
    
    mu = param.mu;
    term_main = 1/mu * log(sum(exp(-mu * sinr_vals)));
    
    lam2 = param.lambda_sinr;
    gamma_th = param.gamma_th;
    alpha = param.alpha;
    pen = 0;
    for k=1:K
        Tk = gamma_th - sinr_vals(k);
        if Tk > 50/alpha, sp=Tk; else, sp=log(1+exp(alpha*Tk))/alpha; end
        pen = pen + lam2 * sp^2;
    end
    J = term_main + pen;
end