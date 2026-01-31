function [W_opt, V_opt, status] = optimize_W_Scheme1_ZF(H_eff, V_current, gamma_vec, sim_param)
% 功能:
%   在给定等效信道 H_eff 下，最小化总发射功率，同时满足用户的 SINR 约束。
% 输入 (Inputs):
%   H_eff       - [M x K] 
%                 等效信道矩阵
%   V_current   - [M x M x K]
%                 初始协方差矩阵
%   gamma_vec   - [K x 1]
%                 目标 SINR 阈值
%   sim_param   - Struct
% 输出 (Outputs):
%   W_opt       - [M x K] 
%                 优化后的波束赋形向量。
%   V_opt       - [M x M x K] 
%                 优化后的协方差矩阵。
%   status      - String
%                 优化状态
    [M, K] = size(H_eff);
    sigma2 = sim_param.sigma2;
    

    MAX_ITER = 10;
    rho = 0.5;
    TOL_RANK = 1e-4; 
    for iter = 1:MAX_ITER
        U_ref = zeros(M, K);
        for k = 1:K
            [eig_vec, eig_val] = eig(V_current(:,:,k));
            [~, idx] = max(diag(real(eig_val)));
            U_ref(:,k) = eig_vec(:, idx);
        end
        
        cvx_begin quiet
        cvx_solver mosek
            variable V(M, M, K) hermitian semidefinite
            
            obj_power = 0;
            obj_penalty = 0;
            for k = 1:K
                obj_power = obj_power + trace(V(:,:,k));
                obj_penalty = obj_penalty + trace(V(:,:,k)) - real(trace(U_ref(:,k)*U_ref(:,k)' * V(:,:,k)));
            end
            
            minimize( obj_power + rho * obj_penalty )
            
            subject to
                for k = 1:K
                    % SINR 约束
                    sig = real(trace(H_eff(:,k)*H_eff(:,k)' * V(:,:,k)));
                    intf = 0;
                    for j = 1:K
                        if j ~= k
                            intf = intf + real(trace(H_eff(:,k)*H_eff(:,k)' * V(:,:,j)));
                        end
                    end
                    sig >= gamma_vec(k) * (intf + sigma2);
                end
        cvx_end
        
        if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
            V_current = full(V);
            
            max_violation = 0;
            for k = 1:K
                eig_vals = eig(V_current(:,:,k));
                max_val = max(eig_vals);
                sum_val = sum(eig_vals);
                max_violation = max(max_violation, abs(sum_val - max_val));
            end
            
            fprintf('  Iter %d: Power=%.2f dBm, RankViol=%.2e\n', iter, 10*log10(cvx_optval), max_violation);
            
            if max_violation < TOL_RANK
                fprintf('  >> Rank-1 achieved. Stopping early.\n');
                status = 'Success';
                break;
            end
        else
            status = 'Failed';
            break;
        end
    end
    
    W_opt = zeros(M, K);
    for k = 1:K
        [u,s] = eig(V_current(:,:,k));
        [v,i] = max(diag(s));
        W_opt(:,k) = sqrt(v) * u(:,i);
    end
    V_opt = V_current;
end