function [H, G_vals_all] = compute_H_with_F(R, F, D, sim)
% 作用:
%   在优化发射模拟波束 F 时，利用固定的旋转矩阵 R 和预计算数据 D，
%   得到有效信道 H。
% 输入 (Inputs):
%   R   - [3 x 3 x M] 当前固定的旋转矩阵。
%   F   - [2 x M]     当前的发射模拟波束 (优化变量)。
%   D   - Struct      预计算的信道几何数据 (d_vecs, psi_vecs, consts)。
%   sim - Struct      仿真参数 (包含指向性因子 p)。
% 输出 (Outputs):
%   H          - [M x K] 有效信道矩阵。
%   G_vals_all - [M x K x L] 所有路径的增益值 
    [~, M] = size(F);
    [~, K, ~, L_total] = size(D.d_vecs);
    H = zeros(M, K);
    G_vals_all = zeros(M, K, L_total);
    
    for m=1:M
        
        r1 = R(:,1,m);
        r2 = R(:,2,m);
        r3 = R(:,3,m);
        
        f_vec = F(:, m);
        
        for k=1:K
            h_sum = 0;
            for l=1:L_total
                d_vec = D.d_vecs(:,k,m,l);
                psi   = D.psi_vecs(:,k,m,l);
                C0    = D.consts(k,m,l);
                
                u_dot = d_vec' * r3;
                G = (max(0, u_dot))^sim.p;
                

                
                coeff_v = psi' * r1; 
                coeff_h = psi' * r2; 
                

                P = coeff_v * f_vec(1) + coeff_h * f_vec(2);
                
                h_sum = h_sum + C0 * G * P;
                G_vals_all(m,k,l) = G;
            end
            H(m,k) = h_sum;
        end
    end
end