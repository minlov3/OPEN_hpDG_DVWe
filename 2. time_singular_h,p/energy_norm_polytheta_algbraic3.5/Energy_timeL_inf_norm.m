function [L_inf_L2_part,L_inf_grad_part,L_inf_H1_part] = Energy_timeL_inf_norm(Space_Node,Time_Node,BDbox , coef ,S_Po,T_Po  ,FEM_index,u,   grad_u,a,zeta)


% information about the bounding box

h = (BDbox(2,:)-BDbox(1,:))./2;  

m = 0.5.*sum(BDbox);

dim_elem =size(FEM_index,1); % number of basis for each element.



%Quadrature rules over rectangle

[Space_P_Qpoints, Space_weights] = quad_rect(Space_Node,ceil((S_Po+1)*0.5));



%% generating the points along t for sampling

t_sampling = linspace(Time_Node(1),Time_Node(2),T_Po)';

L_inf_L2_part = NaN(T_Po,1);
L_inf_grad_part = NaN(T_Po,1);
L_inf_H1_part = NaN(T_Po,1);

% for each sample point in time tm, 
%  we calculate the L2 and H1 norm error over space

for k = 1 : T_Po


%% generator the 3D quad points

% each time we are using different time

Time_P_Qpoints = t_sampling(k);

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
    
    % construct the matrix for all the local basis function
    
     P = zeros(size(P_Qpoints,1) ,dim_elem);
    Px = zeros(size(P_Qpoints,1) ,dim_elem);
    Py = zeros(size(P_Qpoints,1) ,dim_elem);
    Pz = zeros(size(P_Qpoints,1) ,dim_elem);
    
    %shift_leg_derivative(x,m,h,order,0)
    
    for i =1:dim_elem
        
        P(:,i)= FEM2D_DG_basis(P_Qpoints,BDbox,m,h,FEM_index(i,:));
        
        t = grad_FEM2D_DG_basis(P_Qpoints,BDbox,m,h,FEM_index(i,:));
        
        Px(:,i) = t(:,1); Py(:,i) = t(:,2); Pz(:,i) = t(:,3); 
    end
  
    u_DG_val = P*coef;   %DG solution;
    
    grad_u_DG = [Px*coef , Py*coef, Pz*coef];   %gradient of DG
                       
     % Part 1 DG L_2 norm error  int_\kappa a(u - u_DG)^2 dx
         
     % L_2 norm error  int_\kappa (u - u_DG)^2 dx over space    
     
     t1 = (0 - u_DG_val).^2;  
     L_inf_L2_part(k) = dot((t1),weights);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
     % L2 norm of the H_1 semi-norm    
     grad1 = 0 - grad_u_DG;    grad2 = 0 - grad_u_DG;
     
     grad = [grad1(:,1).*grad2(:,1) , grad1(:,1).* grad2(:,2), grad1(:,1).* grad2(:,3)...
             grad1(:,2).*grad2(:,1) , grad1(:,2).* grad2(:,2), grad1(:,2).* grad2(:,3)...
             grad1(:,3).*grad2(:,1) , grad1(:,3).* grad2(:,2), grad1(:,3).* grad2(:,3)];
           
     t2 = sum((grad).*zeta_val,2);   
     L_inf_grad_part(k) = dot((t2),weights);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     % H1=L2+grad
     
     L_inf_H1_part(k)  =   L_inf_grad_part(k) + L_inf_L2_part(k);%%%%%%%%%%%%%%%%
 
 
end
     
end







