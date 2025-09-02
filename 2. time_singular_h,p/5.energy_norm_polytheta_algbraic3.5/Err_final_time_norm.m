function [Final_t_L2_L2_part,Final_t_L2_grad_part, Final_t_L2_H1_part,...
          Final_t_H1_L2_part,Final_t_H1_grad_part, Final_t_H1_H1_part]...
        = Err_final_time_norm(Space_Node,Time_Node,BDbox , coef ,S_Po,T_Po  ,FEM_index,u,grad_u,grad_ut,a,zeta)


% information about the bounding box

h = (BDbox(2,:)-BDbox(1,:))./2;  

m = 0.5.*sum(BDbox);

dim_elem =size(FEM_index,1); % number of basis for each element.



%Quadrature rules over rectangle

[Space_P_Qpoints, Space_weights] = quad_rect(Space_Node,ceil((S_Po+1)*0.5));



%% generating the points along t for sampling

t_sampling = linspace(Time_Node(1),Time_Node(2),T_Po)';

% for each sample point in time tm, 
%  we calculate the L2 and H1 norm error over space


%% generator the 3D quad points

% each time we are using different time

Time_P_Qpoints = Time_Node(2); %%本来是样本插值，我们改为最终时间就好。。。。

quad_x = Space_P_Qpoints(:,1);
quad_y = Space_P_Qpoints(:,2);
quad_t = Time_P_Qpoints.*ones(size(Space_P_Qpoints,1),1); 

weights = Space_weights;

 P_Qpoints = [quad_x,quad_y,quad_t];
    
 

%%Calculating the DG norm error based on the old way
 
 % data for quadrature
     zeta_val = zeta(P_Qpoints);  
     a_val = a(P_Qpoints);   
     u_val = u(P_Qpoints);
     grad_u_val = grad_u(P_Qpoints);
     grad_ut_val = grad_ut(P_Qpoints);
     
    % construct the matrix for all the local basis function
    
     P = zeros(size(P_Qpoints,1) ,dim_elem);
    Px = zeros(size(P_Qpoints,1) ,dim_elem);
    Py = zeros(size(P_Qpoints,1) ,dim_elem);
    Pz = zeros(size(P_Qpoints,1) ,dim_elem);
    
    %shift_leg_derivative(x,m,h,order,0)
    
    for i =1:dim_elem
        
        P(:,i)= FEM2D_DG_basis(P_Qpoints,BDbox,m,h,FEM_index(i,:));
        
        t = grad_FEM2D_DG_basis(P_Qpoints,BDbox,m,h,FEM_index(i,:));
        v = vt_grad_FEM2D_DG_basis(P_Qpoints,BDbox,m,h,FEM_index(i,:));
        
        Px(:,i) = t(:,1); Py(:,i) = t(:,2); Pz(:,i) = t(:,3); 
        Qx(:,i) = v(:,1); Qy(:,i) = v(:,2); Qz(:,i) = v(:,3); 
    end
  
    u_DG_val = P*coef;   %DG solution;
    
    grad_u_DG = [Px*coef , Py*coef, Pz*coef];   %gradient of DG
    grad_ut_DG = [Qx*coef , Qy*coef, Qz*coef];
    
     % Part 1 DG L_2 norm error  int_\kappa a(u - u_DG)^2 dx
         
     %  norm error  int_\kappa (u - u_DG)^2 dx over space    
     
     t1 = (u_val - u_DG_val).^2;          
     Final_t_L2_L2_part  = dot((t1),weights);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     % L2 norm of the H_1 semi-norm    
     grad1 = grad_u_val - grad_u_DG;    grad2 = grad_u_val - grad_u_DG;
     
     grad = [grad1(:,1).*grad2(:,1) , grad1(:,1).* grad2(:,2), grad1(:,1).* grad2(:,3)...
             grad1(:,2).*grad2(:,1) , grad1(:,2).* grad2(:,2), grad1(:,2).* grad2(:,3)...
             grad1(:,3).*grad2(:,1) , grad1(:,3).* grad2(:,2), grad1(:,3).* grad2(:,3)];
           
     t3 = sum((grad).*zeta_val,2);   
     Final_t_L2_grad_part  = dot((t3),weights);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Final_t_L2_H1_part = Final_t_L2_grad_part  ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Final_t_H1_L2_part = dot( grad(:,9),weights);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     % L2 norm of the H_1 semi-norm    
     grad11 = grad_ut_val - grad_ut_DG;    grad21 = grad_ut_val - grad_ut_DG;
     
     grad0 = [grad11(:,1).*grad21(:,1) , grad11(:,1).* grad21(:,2), grad11(:,1).* grad21(:,3)...
              grad11(:,2).*grad21(:,1) , grad11(:,2).* grad21(:,2), grad11(:,2).* grad21(:,3)...
              grad11(:,3).*grad21(:,1) , grad11(:,3).* grad21(:,2), grad11(:,3).* grad21(:,3)];
           
     t4 = sum((grad0).*a_val,2);   
     Final_t_H1_grad_part  = dot((t4),weights);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     
     % H1=L2+grad
     
     Final_t_H1_H1_part   =   Final_t_H1_grad_part  + Final_t_H1_L2_part ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

     
end
