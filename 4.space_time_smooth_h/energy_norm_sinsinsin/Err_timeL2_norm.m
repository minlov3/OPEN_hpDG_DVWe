function [L2_L2_part,L2_grad_part,L2_H1_part ] = Err_timeL2_norm(Space_Node,Time_Node,BDbox , coef ,S_Po,T_Po  ,FEM_index,u,   grad_u,a)


% information about the bounding box

h = (BDbox(2,:)-BDbox(1,:))./2;  

m = 0.5.*sum(BDbox);

dim_elem =size(FEM_index,1); % number of basis for each element.


%Quadrature rules over rectangle

[Space_P_Qpoints, Space_weights] = quad_rect(Space_Node,ceil((S_Po+1)*0.5));



%% time quadracture points

 [w_t,t_points] = quad_GL(ceil((T_Po+2)*0.5));

% map the reference interval [-1,1] to [t_k, t_(k+1)]

Time_weights = w_t.*(Time_Node(2)-Time_Node(1))*0.5;

Time_P_Qpoints = t_points.*(Time_Node(2)-Time_Node(1))*0.5+ 0.5*sum(Time_Node);

%% generator the 3D quad points


quad_x = kron(Space_P_Qpoints(:,1),ones(size(Time_weights,1),1));
quad_y = kron(Space_P_Qpoints(:,2),ones(size(Time_weights,1),1));
quad_t = kron(ones(size(Space_P_Qpoints,1),1),Time_P_Qpoints); 

weights = kron(Space_weights,Time_weights);


 P_Qpoints = [quad_x,quad_y,quad_t];
    
 
%%Calculating the DG norm error based on the old way
 
 % data for quadrature
    
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
         
     % L_2 norm error  int_\kappa (u - u_DG)^2 dx
    
     t1 = (u_val - u_DG_val).^2;
     L2_L2_part = dot((t1),weights);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
     % L2 norm of the H_1 semi-norm
     
     grad1 = grad_u_val - grad_u_DG;    grad2 = grad_u_val - grad_u_DG;
       
     grad = [grad1(:,1).*grad2(:,1) , grad1(:,1).* grad2(:,2), grad1(:,1).* grad2(:,3)...
             grad1(:,2).*grad2(:,1) , grad1(:,2).* grad2(:,2), grad1(:,2).* grad2(:,3)...
             grad1(:,3).*grad2(:,1) , grad1(:,3).* grad2(:,2), grad1(:,3).* grad2(:,3)];
              
     t3 = sum((grad).*a_val,2);  
     L2_grad_part =  dot((t3),weights);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     % H1=L2+grad
     
     L2_H1_part  =  L2_grad_part + L2_L2_part;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
end







