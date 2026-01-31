function [H, G_all, P_all] = compute_H_internal(R, D, sim)
% 功能:
%   在 AO 优化循环内部，利用预计算的信道和当前的旋转矩阵，
%   计算出有效信道
% 输入 (Inputs):
%   R   - [3 x 3 x M] 
%         当前的基站天线旋转矩阵。
%   D   - 信道预计算内容
%   sim - Struct
%         仿真参数
% 输出 (Outputs):
%   H     - [M x K]
%   G_all - [M x K x L_total] 
%   P_all - [M x K x L_total] 
    M = size(R,3); 
    [~, K, ~, L_total] = size(D.d_vecs);
    H = zeros(M, K);
    G_all = zeros(M, K, L_total);
    P_all = zeros(M, K, L_total);
    
    for m=1:M
        r1=R(:,1,m); r2=R(:,2,m); r3=R(:,3,m);
        fv = D.f_analog(1,m); fh = D.f_analog(2,m);
        for k=1:K
            h_sum = 0;
            for l=1:L_total
                d_vec = D.d_vecs(:,k,m,l);
                psi   = D.psi_vecs(:,k,m,l);
                C0    = D.consts(k,m,l);
                
                u_dot = d_vec' * r3;
                G = (max(0, u_dot))^sim.p;
                P = psi' * (fv*r1 + fh*r2);
                
                h_sum = h_sum + C0 * G * P;
                G_all(m,k,l) = G;
                P_all(m,k,l) = P;
            end
            H(m,k) = h_sum;
        end
    end
end