function H_raw = compute_raw_channel_for_U(bs_pos, user_pos, scatter_pos, bs_input, ue_rotations, F_analog, sim_param, nlos_phase_bank)
% 功能:
%   计算用于接收机优化的信道
% 输入 (Inputs):
%   bs_input   - [3xM] 欧拉角 或 [3x3xM] 旋转矩阵
%   F_analog   - [2xM] 发射端极化
% 输出 (Outputs):
%   H_raw      - [2 x M x K]

    [~, M] = size(bs_pos);
    [~, K] = size(user_pos);
    [~, L] = size(scatter_pos);

    is_matrix_input = (ndims(bs_input) == 3 && size(bs_input, 1) == 3 && size(bs_input, 2) == 3);

    H_raw = zeros(2, M, K);
    

    w_los = 1;
    w_nlos = (L > 0) * 1; 

    for k = 1:K
        u_pos = user_pos(:, k);
        u_rot = ue_rotations(:, k);
        R_rx = get_rotation_matrix_3d(u_rot(1), u_rot(2), u_rot(3));
        
        for m = 1:M
            b_pos = bs_pos(:, m);
            f_vec = F_analog(:, m);
            
            if is_matrix_input
                R_tx = bs_input(:, :, m);
            else
                rot_m = bs_input(:, m);
                R_tx = get_rotation_matrix_3d(rot_m(1), rot_m(2), rot_m(3));
            end
            

            vec_path = u_pos - b_pos;
            d_los = norm(vec_path);
            dir_los = vec_path / d_los;
            

            theta_path = acos(dir_los(3));
            phi_path   = atan2(dir_los(2), dir_los(1));
            [e_theta, e_phi] = get_spherical_basis(theta_path, phi_path);
            
            P_los = get_P_matrix(R_tx, e_theta, e_phi);
            Q_los = get_Q_matrix(R_rx, e_theta, e_phi);
            
            phase_los = exp(-1j * 2*pi * d_los / sim_param.lambda);
            M_mat = eye(2) * phase_los;
            
            n_tx = R_tx(:, 3);
            cos_off = dot(n_tx, dir_los);
            if cos_off > 0
                G_val = sim_param.G0 * (cos_off)^(2*sim_param.p);
            else
                G_val = 0;
            end
            
            pl_amp_los = (sim_param.lambda / (4*pi*d_los));

            h_los_2x1 = pl_amp_los * w_los * sqrt(G_val) * (Q_los * M_mat * P_los * f_vec);
            

            h_nlos_sum = zeros(2,1);
            if L > 0
                for l = 1:L
                    s_pos = scatter_pos(:, l);
                    d1 = norm(s_pos - b_pos);
                    d2 = norm(u_pos - s_pos);
                    d_total = d1 + d2;
                    
                    vec_1 = s_pos - b_pos; dir_1 = vec_1 / d1;
                    theta_dep = acos(dir_1(3)); phi_dep = atan2(dir_1(2), dir_1(1));
                    [e_th_dep, e_ph_dep] = get_spherical_basis(theta_dep, phi_dep);
                    
                    vec_2 = u_pos - s_pos; dir_2 = vec_2 / d2;
                    theta_arr = acos(dir_2(3)); phi_arr = atan2(dir_2(2), dir_2(1));
                    [e_th_arr, e_ph_arr] = get_spherical_basis(theta_arr, phi_arr);
                    

                    P_nlos = get_P_matrix(R_tx, e_th_dep, e_ph_dep);
                    Q_nlos = get_Q_matrix(R_rx, e_th_arr, e_ph_arr);
                    

                    phase_prop = exp(-1j * 2*pi * d_total / sim_param.lambda);
                    sqrt_co = sqrt(1 - sim_param.XPD);
                    sqrt_cx = sqrt(sim_param.XPD);
                    rand_phase_1 = nlos_phase_bank(l, m, k, 1);
                    rand_phase_2 = nlos_phase_bank(l, m, k, 2);
                    J_mat = [sqrt_co, sqrt_cx*rand_phase_1; sqrt_cx*rand_phase_2, sqrt_co] * phase_prop;
                    

                    cos_g_nlos = dot(n_tx, dir_1);
                    if cos_g_nlos > 0
                        G_nlos = sim_param.G0 * (cos_g_nlos)^(2*sim_param.p);
                    else
                        G_nlos = 0;
                    end
                    
                    pl_1 = sim_param.lambda / (4*pi*d1);
                    pl_2 = sim_param.lambda / (4*pi*d2);
                    pl_amp_nlos = pl_1 * pl_2;
                    

                    term_vec_nlos = (Q_nlos * J_mat * P_nlos * f_vec);
                    h_nlos_sum = h_nlos_sum + pl_amp_nlos * w_nlos * sqrt(G_nlos) * term_vec_nlos;
                end
            end
            
            H_raw(:, m, k) = h_los_2x1 + h_nlos_sum;
        end
    end
end