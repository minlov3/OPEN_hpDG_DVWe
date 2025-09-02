function z_vect = CR_vect_inflowface_previous_step(Space_Node,Time_Node,  BDbox, previous_BDbox, previous_coef, n ,S_Po,T_Po ,FEM_index, a,zeta, b)


% information about the bounding box


h =  (BDbox(2,:)-BDbox(1,:))./2 ;

m = 0.5.*sum(BDbox);

dim_elem =size(FEM_index,1); % number of basis for each element.

z_vect = zeros(dim_elem,1);  % initialize the local vector


%Quadrature rules over rectangle

[Space_P_Qpoints, Space_weights] = quad_rect(Space_Node,ceil((S_Po+1)*0.5));


%% no need to put time quadracture points, t = Time_Node(1)

weights = Space_weights;


 P_Qpoints = [Space_P_Qpoints,Time_Node(1).*ones(size(Space_P_Qpoints,1),1) ];


 % data for quadrature
    a_val = a(P_Qpoints); 
    zeta_val = zeta(P_Qpoints); 
    b_val = b(P_Qpoints);  
    
    n_vec =kron(n,ones(size(P_Qpoints,1),1));
 
 % impose the inflow boundary condition from previous step  
 % construct the matrix for all the local basis function from previous time
 % step

 
previous_h =  (previous_BDbox(2,:)-previous_BDbox(1,:))./2 ;

previous_m = 0.5.*sum(previous_BDbox);

   P1 = zeros(size(P_Qpoints,1) ,dim_elem);
   P2 = zeros(size(P_Qpoints,1) ,dim_elem);
   P3 = zeros(size(P_Qpoints,1) ,dim_elem); 
    for i =1:dim_elem
        
        P1(:,i)= vx_FEM2D_DG_basis(P_Qpoints,previous_BDbox,previous_m,previous_h,FEM_index(i,:));
        P2(:,i)= vy_FEM2D_DG_basis(P_Qpoints,previous_BDbox,previous_m,previous_h,FEM_index(i,:));
        P3(:,i)= vt_FEM2D_DG_basis(P_Qpoints,previous_BDbox,previous_m,previous_h,FEM_index(i,:));
    end
  
    u_t_val = P3*previous_coef;   %DG solution;   
    u_grad_val= [   P1*previous_coef,   P2*previous_coef,   P3*previous_coef   ];
 

%% constructing the vector    
    
first = zeros(dim_elem,1);
second= zeros(dim_elem,1);

% i is row which is basis of u and j is v. 

   
    for j = 1:dim_elem
    
      % first term b\cdotn *u_previous' *v' 
         
        t_1 =  sum(b_val.*n_vec,2).* u_t_val.*(vt_FEM2D_DG_basis(P_Qpoints,BDbox,m,h,FEM_index(j,:)));
        
        first(j) = dot((t_1),weights);
         
      % first term b\cdotn *g *v          
          
        grad2 = grad_FEM2D_DG_basis(P_Qpoints,BDbox,m,h,FEM_index(j,:));
        grad_u_grad_v =  [u_grad_val(:,1).*grad2(:,1) , u_grad_val(:,1).* grad2(:,2), u_grad_val(:,1).* grad2(:,3)...
                          u_grad_val(:,2).*grad2(:,1) , u_grad_val(:,2).* grad2(:,2), u_grad_val(:,2).* grad2(:,3)...
                          u_grad_val(:,3).*grad2(:,1) , u_grad_val(:,3).* grad2(:,2), u_grad_val(:,3).* grad2(:,3)];    
             
        t_2 = sum((grad_u_grad_v).*zeta_val,2); 
        
        second(j) = dot((t_2),weights);
        
        
    end



z_vect = z_vect + first + second;


end