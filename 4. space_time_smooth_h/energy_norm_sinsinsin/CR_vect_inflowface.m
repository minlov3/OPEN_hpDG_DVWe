function z_vect = CR_vect_inflowface(Space_Node,Time_Node, BDbox, n ,S_Po,T_Po, FEM_index,a,zeta, b, grad_u)



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
grad_u_initial = grad_u(P_Qpoints);
    
    n_vec =kron(n,ones(size(P_Qpoints,1),1));
   
 
first = zeros(dim_elem,1);
second = zeros(dim_elem,1);

% i is row which is basis of u and j is v. 

   
    for j = 1:dim_elem
    
        grad2= grad_FEM2D_DG_basis(P_Qpoints,BDbox,m,h,FEM_index(j,:));
      
        gradu_ini_gradv= [grad_u_initial(:,1).*grad2(:,1) , grad_u_initial(:,1).* grad2(:,2), grad_u_initial(:,1).* grad2(:,3)...
                          grad_u_initial(:,2).*grad2(:,1) , grad_u_initial(:,2).* grad2(:,2), grad_u_initial(:,2).* grad2(:,3)...
                          grad_u_initial(:,3).*grad2(:,1) , grad_u_initial(:,3).* grad2(:,2), grad_u_initial(:,3).* grad2(:,3)];
      
                      
        t_1 =  sum(b_val.*n_vec,2).* gradu_ini_gradv(:,9);
        
        first(j) = dot((t_1),weights);
         
        
        t_2 =  sum((gradu_ini_gradv).*zeta_val,2); 
        
        second(j) = dot((t_2),weights);
                 
    end



z_vect = z_vect + first + second;


end
    