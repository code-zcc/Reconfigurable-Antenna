function [init_power,V_current]=initial_ZF(H_matrix,sigma2,gamma_vec)
% 功能:
%   计算满足 SINR 约束的最小功率 ZF 波束赋形，作为后续复杂迭代算法
%   的初始可行点。
% 输入 (Inputs):
%   H_matrix   - [M x K] 
%                等效信道矩阵 
%   sigma2     - [Scalar] 
%                噪声功率 
%   gamma_vec  - [K x 1] Double
%                每个用户的目标 SINR 阈值。
% 输出 (Outputs):
%   init_power - [Scalar] 
%                满足所有约束所需的总发射功率 (sum(|w_k|^2))。
%   V_current  - [M x M x K] 
%                每个用户的发射协方差矩阵 (V_k = w_k * w_k')。
%                用于初始化基于 SDP 的算法。
W_zf_raw = H_matrix * inv(H_matrix' * H_matrix);
    [M, K] = size(H_matrix);
    W_init = zeros(M, K);
    V_current = zeros(M, M, K);
    
    for k = 1:K

        w_dir = W_zf_raw(:, k);
        w_dir = w_dir / norm(w_dir); 
        channel_gain = abs(H_matrix(:, k)' * w_dir)^2;
        p_required = gamma_vec(k) * sigma2 / channel_gain;
        W_init(:, k) = sqrt(p_required) * w_dir;
        V_current(:, :, k) = W_init(:, k) * W_init(:, k)';
    end
    init_power = sum(vecnorm(W_init).^2);
end