function[bcdof,bcval]=get_boundary_conditions(gcoord,nodeL,nodeR)
%------------------------------------------------------------------------
%  Purpose:
%     Create initial data of boundary condition 
%
%  Synopsis:
%     [bcdof,bcval,ff]=get_initdata
%
%  Variable Description:  
%   bcdof      = vector containing dofs associated with boundary conditions包含与边界条件相关的dofs的向量
%   bcval      = vector containing boundary condition values associated with dofs in bcdof包含bcdof中与dofs相关联的边界条件值的向量

%--------------------------------------------------------------------------
     bcdof=[];
     bcval=[];

     
        %----------------------------------
        %input data for boundary conditions%边界条件！！
        %----------------------------------
       
        for i=1:length(nodeL)
            bcdof(1,i)=nodeL(i,1);%x=0那一条边
            bcval(1,i)=0;   
        end
        for i=1:length(nodeR)
            bcdof(1,length(nodeL)+i)=nodeR(i,1);%x=xn那一条边
            bcval(1,length(nodeL)+i)=1000;
        end
%         for i=1:length(nodeB)
%             bcdof(1,length(nodeL)+length(nodeR)+i)=nodeB(i,1);%y=0那一条边
%             bcval(1,length(nodeL)+length(nodeR)+i)=0;
%             
%             
%       end
%         for i=1:length(nodeT)
%             bcdof(1,length(nodeL)+length(nodeR)+length(nodeB)+i)=nodeT(i,1);%y=ym那一条边
%             bcval(1,length(nodeL)+length(nodeR)+length(nodeB)+i)=0;
%         end
   
        %每一条边的边界值，也就是边界条件
end