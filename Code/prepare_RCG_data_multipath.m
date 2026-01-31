function data = prepare_RCG_data_multipath(bs_pos, user_pos, scatter_pos, ue_rotations, sim_param, nlos_phase_bank, U_user, F_analog)
    [~, M] = size(bs_pos);
    [~, K] = size(user_pos);
    [~, L] = size(scatter_pos);
    
    total_paths = 1 + (L > 0) * L;
    
    data.d_vecs = zeros(3, K, M, total_paths);
    data.psi_vecs = zeros(3, K, M, total_paths);
    data.consts = zeros(K, M, total_paths);
    data.f_analog = F_analog; 
    data.num_paths = total_paths;
    
    w_los = 1; w_nlos = (L>0)*1;
    
    amp_G0 = sqrt(sim_param.G0); 
    
    for k=1:K
        u_p = user_pos(:,k);
        u_vec = U_user(:,k);
        u_rot = ue_rotations(:,k);
        R_rx = get_rotation_matrix_3d(u_rot(1), u_rot(2), u_rot(3));
        
        for m=1:M
            b_p = bs_pos(:,m);
            
            vec = u_p - b_p; d = norm(vec); dir = vec/d;
            data.d_vecs(:,k,m,1) = dir;
            
            [e_th, e_ph] = get_spherical_basis_vec(dir);
            Q = get_Q_matrix(R_rx, e_th, e_ph);
            
            M_mat = eye(2) * exp(-1j*2*pi*d/sim_param.lambda);
            
            vec_rx = u_vec' * Q * M_mat;
            psi = vec_rx(1)*e_th + vec_rx(2)*e_ph;
            data.psi_vecs(:,k,m,1) = conj(psi);
            
            pl = (sim_param.lambda/(4*pi*d));
            

            data.consts(k,m,1) = pl * w_los * amp_G0;
            

            if L>0
                for l=1:L
                    s_p = scatter_pos(:,l);
                    d1 = norm(s_p - b_p); d2 = norm(u_p - s_p);
                    vec1 = s_p - b_p; dir1 = vec1/d1;
                    
                    data.d_vecs(:,k,m,1+l) = dir1;
                    
                    [e_th1, e_ph1] = get_spherical_basis_vec(dir1);
                    vec2 = u_p - s_p; dir2 = vec2/d2;
                    [ e_th2, e_ph2] = get_spherical_basis_vec(dir2);
                    
                    Q_n = get_Q_matrix(R_rx, e_th2, e_ph2);
                    
                    pp = exp(-1j*2*pi*(d1+d2)/sim_param.lambda);
                    sq_co = sqrt(1-sim_param.XPD); sq_cx = sqrt(sim_param.XPD);
                    rp1 = nlos_phase_bank(l,m,k,1); rp2 = nlos_phase_bank(l,m,k,2);
                    J = [sq_co, sq_cx*rp1; sq_cx*rp2, sq_co] * pp;
                    
                    vec_rx_n = u_vec' * Q_n * J;
                    psi_n = vec_rx_n(1)*e_th1 + vec_rx_n(2)*e_ph1;
                    data.psi_vecs(:,k,m,1+l) = conj(psi_n);
                    
                    pl1 = sim_param.lambda/(4*pi*d1);
                    pl2 = sim_param.lambda/(4*pi*d2);
                    
                    data.consts(k,m,1+l) = (pl1*pl2) * w_nlos * amp_G0;
                end
            end
        end
    end
end