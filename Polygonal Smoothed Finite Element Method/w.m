function [phi,dphi] = wachspress2d(v,x)
% Reference: GRADIENT BOUNDS FOR WACHSPRESS COORDINATES ON POLYTOPES
%            MICHAEL S. FLOATER, ANDREW GILLETTE, N. SUKUMAR
%
% Evaluate Wachspress basis functions in a convex polygon


%评估凸多边形中的 Wachspress 基函数
% Inputs:
% v : [x1 y1; x2 y2; ...; xn yn], the n vertices of the polygon in ccw
% x : [x(1) x(2)], the point at which the basis functions are computed
% Outputs:
% phi : output basis functions = [phi_1; ...; phi_n]

n = size(v,1);%总共多少个点
w = zeros(n,1);
phi = zeros(n,1);
dphi = zeros(n,2);
for i = 1:n
    d = v(mod(i,n)+1,:) - v(i,:);%对应坐标相减
    normal_vector = [d(2) -d(1)]/norm(d);%标准化向量
    h = (v(i,:)-x) * normal_vector'; % dot product%点积（点到边的距离）
    if h==0,h=1e-10;end
    p(i,:) = normal_vector / h;  % scaled normal vectors%外法向量

end
for i = 1:n
    im1 = mod(i-2,n) + 1;
    w(i) = det([p(im1,:);p(i,:)]);
    R(i,:) = p(im1,:) + p(i,:);
end
wsum = sum(w);
phi = w/wsum;%形函数计算公式

phiR = phi' * R;
for k = 1:2
dphi(:,k) = phi .* (R(:,k) - phiR(:,k));
end

