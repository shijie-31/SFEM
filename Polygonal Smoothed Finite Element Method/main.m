%Polygonal Smoothed Finite Element Method (CS‑FEM)--RPIM
%   --------------------------------------------------------------------
%   (c) 2025 Shijie Zhao <shijie31123@163.com>
%   --------------------------------------------------------------------

%% Load mesh data
load('100.mat')
 ele_nods=elem;
   gcoord=node;

%% Boundary conditions
nodeL = find(node(:,1)<0+0.001);
nodeR  = find(node(:,1)>1-0.001);
[bcdof,bcval]=get_boundary_condition(gcoord,nodeL,nodeR);
[K]=cal_K(ele_nods,gcoord);
Q=sparse( length(node(:,1)),1);

%% Apply boundary conditions & solve
[K,Q]=apply(K,Q,bcdof,bcval);
U=K\Q;

%% Visualisation
figure(1)
showsolution(node, ele_nods,U);
c = colorbar;
 axis off;
 c.FontSize = 16; % Replace 14 with your desired font size
c.FontWeight = 'bold'; % Make the font bold
c.FontName = 'Arial';  % Change the font type

figure(2)
showmesh(node,elem)
