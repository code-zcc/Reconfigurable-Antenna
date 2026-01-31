function [F_opt, min_sinr] = run_F_step(F_init, R_curr, W, H_func_data, sim_param)
% 目的 (Purpose):
%   在给定旋转矩阵 R 和数字波束 W 的情况下，优化发射端的模拟极化波束 F。
% 输入 (Inputs):
%   F_init      - [2 x M] 
%                 F 的初始值。每一列必须满足单位范数约束。
%   R_curr      - [3 x 3 x M]
%                 当前固定的旋转矩阵。
%   W           - [M x K] 
%                 当前固定的数字波束赋形矩阵。
%   H_func_data - Struct
%                 包含预计算的信道几何数据
%   sim_param   - Struct
%                 仿真参数。
% 输出 (Outputs):
%   F_opt       - [2 x M] 优化后的模拟波束矩阵。
%   min_sinr    -  优化后的最小用户 SINR ，用于收敛判断。
    problem.M = obliquecomplexfactory(2, sim_param.M);
    
    problem.cost = @(F) calc_objective_F(F, R_curr, H_func_data, W, sim_param);
    cost_handle=problem.cost;
    problem.egrad = @(F) egrad_F_RCG(F, R_curr, H_func_data, W, sim_param); 

    options.maxiter = 20;      
    options.tolgradnorm = 1e-3;
    options.verbosity = 0;      
    options.beta_type = 'P-R';

    grad_euclidean = problem.egrad(F_init);
    grad_riemannian = problem.M.egrad2rgrad(F_init, grad_euclidean);
    fprintf('Euclidean Grad Norm: %e\n', norm(grad_euclidean(:)));
    fprintf('Riemannian Grad Norm: %e\n', norm(grad_riemannian(:)));

    [F_opt, cost_min, ~] = conjugategradient(problem, F_init, options);
 
    [~, sinr_vals] = calc_objective_F(F_opt, R_curr, H_func_data, W, sim_param);
    min_sinr = min(sinr_vals);

end