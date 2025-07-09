function [K,ff]=apply_bcdof(K,ff,bcdof,bcval)

%----------------------------------------------------------
%  Purpose:
%     Apply constraints to matrix equation [K]{disp}={ff}��Լ��Ӧ���ھ��󷽳�[k]{disp}={ff}
%
%  Synopsis:ժҪ
%     [K,f]=feaplybc(K,ff,bcdof,bcval)
%
%  Variable Description:
%     K - system matrix before applying constraints Ӧ��Լ��ǰ��ϵͳ����
%     ff - system vector before applying constraintsӦ��Լ��ǰ��ϵͳ����
%     bcdof - a vector containging constrained d.o.f
%     bcval - a vector containing constrained value ����Լ��ֵ������ 
%
%     For example, there are constraints at d.o.f=2 and 10
%     and their constrained values are 0.0 and 2.5, 
%     respectively.  Then, bcdof(1)=2 and bcdof(2)=10; and
%     bcval(1)=1.0 and bcval(2)=2.5.
%--------------------------------------------------------------------------
% Coded by Dr. Nguyen Thoi Trung (Nguyen-Thoi T or Nguyen T.T)            %
% University of Science - Vietnam National University ?HCMC, Vietnam     %
% National University of Singapore (NUS)                                  %
% email: thoitrung76@gmail.com                                            %
% Last modified: December 2009                                            %
%--------------------------------------------------------------------------
 
% Important note: The authors decided to release the source codes free of charge with the hope that the S-FEM technique 
% can be applied to more problems and can be further developed into even more powerful methods. 
% The authors are not be able to provide any services or user-guides, but appreciate very much your feedback on errors and suggestions, 
% so that we can improve these codes/methods in the future. If the idea, method, and any part of these codes are used in anyway,
% the users are required to cite the book and the following related original papers of the authors:

% Liu GR, The Finite element method ?a practical course, Elsevier (BH), UK.  2003.
% Liu, G. R. and Nguyen Thoi Trung, Smoothed Finite Element Method, CRC press, Boca Raton, USA, 2010.
%--------------------------------------------------------------------------

n=length(bcdof);
 sdof=size(K);%size���ص���һ����������һ��Ԫ���Ǿ�����������ڶ���Ԫ���Ǿ��������
%���Ŀ���ǰ���Լ��������������һ,��ȫ���0Ȼ��K��c,c����λ�ñ��1���޸�ff��ֵ�������ȴ����о���ȷ�����Ƶ��¶�
 for i=1:n
    c=bcdof(i);
    for j=1:sdof
       K(c,j)=0;
    end
    K(c,c)=1;
    ff(c)=bcval(i);
 end

