function P = get_P_matrix(R_tx, e_theta, e_phi)
% 功能:
%   计算从发射天线本地极化端口到波传播坐标系(e_theta, e_phi)的投影矩阵。
% 输入:
%   R_tx    - [3 x 3] 发射天线的旋转矩阵
%   e_theta - [3 x 1] 出射波的 theta 分量基向量
%   e_phi   - [3 x 1] 出射波的 phi 分量基向量
% 输出:
%   P       - [2 x 2] 投影矩阵
    u_loc_1 = [1; 0; 0]; 
    u_loc_2 = [0; 1; 0];
    

    u_glob_1 = R_tx * u_loc_1;
    u_glob_2 = R_tx * u_loc_2;
    

    P = [dot(e_theta, u_glob_1), dot(e_theta, u_glob_2);
         dot(e_phi,   u_glob_1), dot(e_phi,   u_glob_2)];
end