function[bcdof,bcval]=get_boundary_conditions(gcoord,nodeL,nodeR)
%------------------------------------------------------------------------
%  Purpose:
%     Create initial data of boundary condition 
%
%  Synopsis:
%     [bcdof,bcval,ff]=get_initdata
%
%  Variable Description:  
%   bcdof      = vector containing dofs associated with boundary conditions������߽�������ص�dofs������
%   bcval      = vector containing boundary condition values associated with dofs in bcdof����bcdof����dofs������ı߽�����ֵ������

%--------------------------------------------------------------------------
     bcdof=[];
     bcval=[];

     
        %----------------------------------
        %input data for boundary conditions%�߽���������
        %----------------------------------
       
        for i=1:length(nodeL)
            bcdof(1,i)=nodeL(i,1);%x=0��һ����
            bcval(1,i)=0;   
        end
        for i=1:length(nodeR)
            bcdof(1,length(nodeL)+i)=nodeR(i,1);%x=xn��һ����
            bcval(1,length(nodeL)+i)=1000;
        end
%         for i=1:length(nodeB)
%             bcdof(1,length(nodeL)+length(nodeR)+i)=nodeB(i,1);%y=0��һ����
%             bcval(1,length(nodeL)+length(nodeR)+i)=0;
%             
%             
%       end
%         for i=1:length(nodeT)
%             bcdof(1,length(nodeL)+length(nodeR)+length(nodeB)+i)=nodeT(i,1);%y=ym��һ����
%             bcval(1,length(nodeL)+length(nodeR)+length(nodeB)+i)=0;
%         end
   
        %ÿһ���ߵı߽�ֵ��Ҳ���Ǳ߽�����
end