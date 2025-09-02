function [Time_jump_u_part, Time_jump_gradu_part, Time_jump_dotu_part]...
    = Err_Jump_time_norm(Space_Node,Time_Node,Time_Node1,BDbox,BDbox1 , coef,coef1 ,S_Po,T_Po  ,FEM_index,a)


% information about the bounding box

h = (BDbox(2,:)-BDbox(1,:))./2;  
h1 = (BDbox1(2,:)-BDbox1(1,:))./2;  
m = 0.5.*sum(BDbox);
m1 = 0.5.*sum(BDbox1);
dim_elem =size(FEM_index,1); % number of basis for each element.



%Quadrature rules over rectangle

[Space_P_Qpoints, Space_weights] = quad_rect(Space_Node,ceil((S_Po+1)*0.5));



%% generating the points along t for sampling

t_sampling = linspace(Time_Node(1),Time_Node(2),T_Po)';

% for each sample point in time tm, 
%  we calculate the L2 and H1 norm error over space


%% generator the 3D quad points

% each time we are using different time

Time_P_Qpoints = Time_Node(1); %%本来是样本插值，我们改为最终时间就好。。。。


quad_x = Space_P_Qpoints(:,1);
quad_y = Space_P_Qpoints(:,2);
quad_t = Time_P_Qpoints.*ones(size(Space_P_Qpoints,1),1); 

weights = Space_weights;

 P_Qpoints = [quad_x,quad_y,quad_t];
    
 

%%Calculating the DG norm error based on the old way
 
 % data for quadrature
    
     a_val = a(P_Qpoints);   

     
    % construct the matrix for all the local basis function
    
     P = zeros(size(P_Qpoints,1) ,dim_elem);
     P1 = zeros(size(P_Qpoints,1) ,dim_elem);
    Px = zeros(size(P_Qpoints,1) ,dim_elem);
    P1x = zeros(size(P_Qpoints,1) ,dim_elem);
    Py = zeros(size(P_Qpoints,1) ,dim_elem);
     P1y = zeros(size(P_Qpoints,1) ,dim_elem);
    Pz = zeros(size(P_Qpoints,1) ,dim_elem);
     P1z = zeros(size(P_Qpoints,1) ,dim_elem);
    %shift_leg_derivative(x,m,h,order,0)
    
    for i =1:dim_elem
        
        P(:,i)= FEM2D_DG_basis(P_Qpoints,BDbox,m,h,FEM_index(i,:));
        t = grad_FEM2D_DG_basis(P_Qpoints,BDbox,m,h,FEM_index(i,:));
        v = vt_grad_FEM2D_DG_basis(P_Qpoints,BDbox,m,h,FEM_index(i,:));
        
        P1(:,i)= FEM2D_DG_basis(P_Qpoints,BDbox1,m1,h1,FEM_index(i,:));
        t1 = grad_FEM2D_DG_basis(P_Qpoints,BDbox1,m1,h1,FEM_index(i,:));
        v1 = vt_grad_FEM2D_DG_basis(P_Qpoints,BDbox1,m1,h1,FEM_index(i,:));
        
        Px(:,i) = t(:,1); Py(:,i) = t(:,2); Pz(:,i) = t(:,3); 
        Qx(:,i) = v(:,1); Qy(:,i) = v(:,2); Qz(:,i) = v(:,3); 
        
        P1x(:,i) = t1(:,1); P1y(:,i) = t1(:,2); P1z(:,i) = t1(:,3); 
        Q1x(:,i) = v1(:,1); Q1y(:,i) = v1(:,2); Q1z(:,i) = v1(:,3); 
    end
  
    u_DG_val_1 = P*coef;   %DG solution;
    u_DG_val_2 = P1*coef1;   %DG solution;
    grad_u_DG_1 = [Px*coef , Py*coef, Pz*coef];   %gradient of DG
    grad_u_DG_2 = [P1x*coef1 , P1y*coef1, P1z*coef1];   %gradient of DG
    grad_ut_DG_1 = [Qx*coef , Qy*coef, Qz*coef];
    grad_ut_DG_1 = [Q1x*coef1 , Q1y*coef1, Q1z*coef1];
    
     % Part 1 DG L_2 norm error  int_\kappa a(u - u_DG)^2 dx
         
     %  norm error  int_\kappa (u - u_DG)^2 dx over space    
     
     t1 = (u_DG_val_1 - u_DG_val_2).^2;          
     Time_jump_u_part  = dot((t1),weights);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     % L2 norm of the H_1 semi-norm    
     grad1 = grad_u_DG_1 - grad_u_DG_2;    grad2 = grad1;
     
     grad = [grad1(:,1).*grad2(:,1) , grad1(:,1).* grad2(:,2), grad1(:,1).* grad2(:,3)...
             grad1(:,2).*grad2(:,1) , grad1(:,2).* grad2(:,2), grad1(:,2).* grad2(:,3)...
             grad1(:,3).*grad2(:,1) , grad1(:,3).* grad2(:,2), grad1(:,3).* grad2(:,3)];
           
     t3 = sum((grad).*a_val,2);   
     Time_jump_gradu_part  = dot((t3),weights);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     

     Time_jump_dotu_part = dot(grad(:,9),weights);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     

     
end
