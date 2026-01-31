function [J, sinr_vals] = calc_objective_F(F, R, D, W, param)
% 功能:
%   计算在给定发射模拟波束 F 的情况下的系统总代价 J。
% 输入 (Inputs):
%   F      - [2 x M] 
%            当前的发射模拟波束矩阵
%   R      - [3 x 3 x M]
%            当前固定的旋转矩阵。
%   D      - Struct
%            预计算的信道几何数据
%   W      - [M x K] 
%            当前固定的数字波束赋形矩阵。
%   param  - Struct
%            仿真参数 
% 输出 (Outputs):
%   J      -  总目标函数值 (Minimize this)。
%   sinr_vals -   当前所有用户的实际 SINR 线性值。
    [H, ~] = compute_H_with_F(R, F, D, param);
    

    [~, K] = size(W);
    S_vals = zeros(K,1); I_vals = zeros(K,1); sinr_vals = zeros(K,1);
    for k=1:K
        S_vals(k) = abs(W(:,k)' * H(:,k))^2;
        intf = 0;
        for j=1:K, if j~=k, intf = intf + abs(W(:,j)' * H(:,k))^2; end, end
        I_vals(k) = intf + param.sigma2;
        sinr_vals(k) = S_vals(k) / I_vals(k);
    end
    

    mu = param.mu;
    term_main = 1/mu * log(sum(exp(-mu * sinr_vals)));

    lam2 = param.lambda_sinr;
    gamma_th = param.gamma_th;
    alpha = param.alpha;
    pen_sinr = 0;
    for k=1:K
        Tk = gamma_th - sinr_vals(k);
        if Tk > 50/alpha, sp=Tk; else, sp=log(1+exp(alpha*Tk))/alpha; end
        pen_sinr = pen_sinr + lam2 * sp^2;
    end
    
    J = term_main + pen_sinr;
end