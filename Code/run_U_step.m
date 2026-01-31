function [U_opt, min_sinr] = run_U_step(U_init, H_raw, W, sim_param)
% 功能:
%   利用RCG算法来更新U
%Inputs:
%  U_init: [2 x K] 初始值
%  H_raw:  [2 x M x K] 信道
%  W:      [M x K] 当前发射波束
%Output:
%  U_opt:  [2 x K] 优化后的接收波束
    problem.M = complexcirclefactory(2, size(U_init, 2));
    

    problem.cost = @(U) calc_objective_U(U, H_raw, W, sim_param);
    problem.egrad = @(U) egrad_U_RCG(U, H_raw, W, sim_param);

    options.maxiter = 20;
    options.tolgradnorm = 1e-3;
    options.verbosity = 0; 
    [U_opt, cost_min, ~] = conjugategradient(problem, U_init, options);
    [~, sinr_vals] = calc_objective_U(U_opt, H_raw, W, sim_param);
    min_sinr = min(sinr_vals);
end

% function [U_opt, min_sinr] = run_U_step(U_init, H_raw, W, sim_param)
% % RUN_U_STEP 使用 RCG 迭代优化接收波束 U
% % 
% % Inputs:
% %   U_init: [2 x K] 初始值
% %   H_raw:  [2 x M x K] 原始到达信道 (由 compute_raw_channel_for_U 计算)
% %   W:      [M x K] 当前发射波束
% %
% % Output:
% %   U_opt:  [2 x K] 优化后的接收波束
% 
%     % 1. 定义流形: K 个独立的 2维复数圆 (Complex Circle)
%     %    这意味着 U 的每一列 (每个用户) 模长必须为 1
%     problem.M = complexcirclefactory(2, size(U_init, 2));
%     
%     % 2. 定义目标函数 (Soft-min SINR + Penalty)
%     problem.cost = @(U) calc_objective_U(U, H_raw, W, sim_param);
%     cost_handle=problem.cost;
%     % 3. 定义欧式梯度
%     problem.egrad = @(U)  get_numeric_gradient(cost_handle, U);
%     
%     % 4. 求解器配置
%     options.maxiter = 20;           % 不需要太多迭代，AO外层会循环
%     options.tolgradnorm = 1e-3;
%     options.verbosity = 0;          % 静默模式
%     
%     % 5. 求解
%     [U_opt, ~, ~] = conjugategradient(problem, U_init, options);
%     
%     % 6. 返回结果
%     [~, sinr_vals] = calc_objective_U(U_opt, H_raw, W, sim_param);
%     min_sinr = min(sinr_vals);
%     
%     % fprintf('  U-Step (RCG) Done. Min SINR: %.2f dB\n', 10*log10(min_sinr));
% end