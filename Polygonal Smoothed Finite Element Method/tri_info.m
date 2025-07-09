function [midPts, areas, edgeLens, edgeNorms] = tri_info(p, tri)
% TRI_INFO  计算二维三角网格的几何信息
% ------------------------------------------------------------
% 输入:
%   p    : m×2  节点坐标 (x, y)
%   tri  : k×3  三角形连接表 (每行为三个顶点索引)
%
% 输出:
%   midPts    : k×3×2  三条边中点坐标 (按 AB, BC, CA 顺序)
%   areas     : k×1    三角形面积 (正值)
%   edgeLens  : k×3    三条边长度
%   edgeNorms : k×3×2  三条边的外法向量 (单位向量, 按 AB, BC, CA 顺序)
% ------------------------------------------------------------
% 最后更新: 2025‑05‑07

% ---------- 预分配 ----------
k = size(tri, 1);
midPts    = zeros(k, 3, 2);    % 中点
areas     = zeros(k, 1);       % 面积
edgeLens  = zeros(k, 3);       % 边长
edgeNorms = zeros(k, 3, 2);    % 外法向
epsLen    = 1e-12;             % 边长下限 (避免除 0)

for i = 1:k
    % ---- 顶点 ----
    v = tri(i, :);
    A = p(v(1), :);
    B = p(v(2), :);
    C = p(v(3), :);

    % ---- 边向量 (AB, BC, CA) ----
    e1 = B - A;   % AB
    e2 = C - B;   % BC
    e3 = A - C;   % CA

    % ---- 边长 ----
    L1 = max(norm(e1), epsLen);
    L2 = max(norm(e2), epsLen);
    L3 = max(norm(e3), epsLen);
    edgeLens(i, :) = [L1, L2, L3];

    % ---- 面积 (带符号) ----
    twiceAreaSigned = e1(1)*e2(2) - e1(2)*e2(1);   % 2*area (AB×BC)
    areas(i)        = 0.5 * abs(twiceAreaSigned);

    % ---- 边中点 ----
    midPts(i, 1, :) = (A + B)/2;   % AB
    midPts(i, 2, :) = (B + C)/2;   % BC
    midPts(i, 3, :) = (C + A)/2;   % CA

    % ---- 初始法向 (右手法向, 未保证外向) ----
    n1 = [ e1(2), -e1(1) ] / L1;
    n2 = [ e2(2), -e2(1) ] / L2;
    n3 = [ e3(2), -e3(1) ] / L3;

    % ---- 统一指向外侧 ----
    % 以三角形质心为参考, 若 n·(质心-中点) > 0 则说明 n 指向内侧, 需翻转
    G  = (A + B + C) / 3;          % 质心
    M1 = squeeze(midPts(i, 1, :)).';   % AB 中点
    M2 = squeeze(midPts(i, 2, :)).';   % BC 中点
    M3 = squeeze(midPts(i, 3, :)).';   % CA 中点

    if dot(n1, G - M1) > 0, n1 = -n1; end
    if dot(n2, G - M2) > 0, n2 = -n2; end
    if dot(n3, G - M3) > 0, n3 = -n3; end

    edgeNorms(i, 1, :) = n1;
    edgeNorms(i, 2, :) = n2;
    edgeNorms(i, 3, :) = n3;
end
end
