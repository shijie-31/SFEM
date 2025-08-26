function [phi, dphi] = rpim_shape_function(v, x)
% 计算 RPIM（径向点插值法）下的形函数值和梯度 (包含二次多项式基函数)
% 输入参数：
% v : [x1 y1; x2 y2; x3 y3]，三角形单元的三个顶点坐标
% x : [x, y]，待计算形函数的点的坐标
% 输出参数：
% phi : 形函数值向量 [phi1; phi2; phi3]
% dphi : 形函数的梯度矩阵 [dphi1_dx, dphi1_dy; dphi2_dx, dphi2_dy; dphi3_dx, dphi3_dy]

% 参数设置
ac = 0.35;
q = 1.03;

% 计算 dc，节点间最大距离
dc1 = sqrt((v(1,1) - v(2,1))^2 + (v(1,2) - v(2,2))^2);
dc2 = sqrt((v(2,1) - v(3,1))^2 + (v(2,2) - v(3,2))^2);
dc3 = sqrt((v(3,1) - v(1,1))^2 + (v(3,2) - v(1,2))^2);
dc = max([dc1, dc2, dc3]);

% 定义径向基函数及其导数
% RBF: R(r) = (r^2 + (ac*dc)^2)^q
% r^2 = (xx - xi)^2 + (yy - yi)^2
rbf = @(xx, xi, yy, yi) ((xx - xi)^2 + (yy - yi)^2 + (ac * dc)^2)^q;
% Derivative of RBF w.r.t r^2: dR/d(r^2) = q * (r^2 + (ac*dc)^2)^(q-1)
% So, dR/dx = dR/d(r^2) * d(r^2)/dx = q * (...) * 2 * (xx - xi)
% And, dR/dy = dR/d(r^2) * d(r^2)/dy = q * (...) * 2 * (yy - yi)
drbf_common_term = @(xx, xi, yy, yi) 2 * q * ((xx - xi)^2 + (yy - yi)^2 + (ac * dc)^2)^(q - 1);


% 节点数
wl = size(v, 1);

% 构建 RBF 矩阵 Rq 和多项式矩阵 Pm
% Pm 的列数变为 6，对应 [1, x, y, x^2, y^2, xy]
num_poly_terms = 6; % 1 + 2 + 3 = 6 (常数项 + 线性项 + 二次项)
Rq = zeros(wl, wl);
Pm = zeros(wl, num_poly_terms);
for k = 1:wl
    for j = 1:wl
        Rq(k, j) = rbf(v(k,1), v(j,1), v(k,2), v(j,2));
    end
    % 二次多项式基函数: [1, x, y, x^2, y^2, xy]
    Pm(k, :) = [1, v(k,1), v(k,2), v(k,1)^2, v(k,2)^2, v(k,1)*v(k,2)];
end

% 检查矩阵可逆性，避免潜在问题
if rcond(Rq) < eps % 检查条件数，如果太小，表示矩阵接近奇异
    warning('RBF矩阵Rq接近奇异，可能导致计算不稳定。请检查节点分布或参数设置。');
end

% 计算形函数系数 Sa 和 Sb
% 注意：这里假设 (Pm' * inv(Rq) * Pm) 是可逆的。
% 在节点数量较少时，可能出现 Pm' * inv(Rq) * Pm 奇异的情况。
% 此时需要更多节点或者使用伪逆。
try
    invRq = inv(Rq);
    Sb = inv(Pm' * invRq * Pm) * Pm' * invRq;
    Sa = invRq - invRq * Pm * Sb;
catch ME
    if strcmp(ME.identifier, 'MATLAB:singularMatrix')
        error('多项式基函数矩阵 Pm 引起的矩阵奇异。这通常发生在节点数量不足以唯一确定二次多项式时 (例如，只有3个节点)。对于二次多项式，通常需要至少6个节点。');
    else
        rethrow(ME); % 重新抛出其他错误
    end
end


% 在点 x 处计算形函数值
Rqsg = zeros(1, wl);
for w = 1:wl
    Rqsg(w) = rbf(x(1), v(w,1), x(2), v(w,2));
end
% 在点 x 处计算多项式基函数值
Pmsg = [1, x(1), x(2), x(1)^2, x(2)^2, x(1)*x(2)];
phi = (Rqsg * Sa + Pmsg * Sb)';

% 计算形函数的梯度
dRqsg = zeros(wl, 2);
for w = 1:wl
    dx = x(1) - v(w,1);
    dy = x(2) - v(w,2);
    common_term = drbf_common_term(x(1), v(w,1), x(2), v(w,2));
    dRqsg(w,1) = common_term * dx; % dR/dx
    dRqsg(w,2) = common_term * dy; % dR/dy
end

% 在点 x 处计算多项式基函数的梯度
% Pmsg = [1, x, y, x^2, y^2, xy]
% dPmsg_dx = [0, 1, 0, 2x, 0, y]
% dPmsg_dy = [0, 0, 1, 0, 2y, x]
dPmsg = zeros(2, num_poly_terms);
dPmsg(1, :) = [0, 1, 0, 2*x(1), 0, x(2)]; % 对 x 的偏导
dPmsg(2, :) = [0, 0, 1, 0, 2*x(2), x(1)]; % 对 y 的偏导

dphi = zeros(wl, 2);
for dim = 1:2 % dim = 1 for dx, dim = 2 for dy
    dphi(:, dim) = (dRqsg(:, dim)' * Sa + dPmsg(dim, :) * Sb)';
end

end
