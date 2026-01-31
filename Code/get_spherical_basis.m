function [e_theta, e_phi] = get_spherical_basis(theta, phi)
% 功能:
%   根据给定的(theta)和(phi)，计算波传播方向上的正交基
% 输入:
%   theta - [Scalar]
%   phi   - [Scalar]
% 输出:
%   e_theta - [3 x 1] 
%   e_phi   - [3 x 1] 
    e_theta = [cos(theta)*cos(phi); cos(theta)*sin(phi); -sin(theta)];
    e_phi   = [-sin(phi); cos(phi); 0];
end