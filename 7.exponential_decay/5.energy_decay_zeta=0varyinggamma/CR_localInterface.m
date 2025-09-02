function  z = CR_localInterface(Space_Node,Time_Node,BDbox1,BDbox2, n ,Po ,Lege_ind, b)



h1 = (BDbox1(2,:)-BDbox1(1,:))./2;  h2 = (BDbox2(2,:)-BDbox2(1,:))./2;


m1 =  0.5.*sum(BDbox1);     m2 =  0.5.*sum(BDbox2); 

dim_elem =size(Lege_ind,1); % number of basis for each element.



% generating quadrature points and weights for 2D space

[weights,ref_Qpoints] = quad_GL(ceil((Po+1)*0.5));

%change the quadrature nodes from reference domain to physical domain.

Space_mid= sum(Space_Node)./size(Space_Node,1);   tanvec = 0.5* (Space_Node(2,:)-Space_Node(1,:));

C = kron(Space_mid,ones(size(ref_Qpoints,1),1));

Space_P_Qpoints = kron(ref_Qpoints,tanvec) + C; 

Space_weights = weights.*norm((Space_Node(2,:)-Space_Node(1,:))).*0.5;


%% time quadracture points

 [w_t,t_points] = quad_GL(ceil((Po+1)*0.5));

% map the reference interval [-1,1] to [t_k, t_(k+1)]

Time_weights = w_t.*(Time_Node(2)-Time_Node(1))*0.5;

Time_P_Qpoints = t_points.*(Time_Node(2)-Time_Node(1))*0.5+ 0.5*sum(Time_Node);



%% generator the 3D quad points


quad_x = kron(Space_P_Qpoints(:,1),ones(size(Time_weights,1),1));

quad_y = kron(Space_P_Qpoints(:,2),ones(size(Time_weights,1),1));

quad_t = kron(ones(size(Space_P_Qpoints,1),1),Time_P_Qpoints); 

weights = kron(Space_weights,Time_weights);


P_Qpoints = [quad_x,quad_y,quad_t];

 % data for quadrature, function value b and normal vector n_vec
    
b_val = b(P_Qpoints); n_vec =kron(n,ones(size(P_Qpoints,1),1));
    
    
    
    
 % first is the local elememt Kappa  second is in the neighbour element;

local_1 = zeros(dim_elem,dim_elem); local_2 = zeros(dim_elem,dim_elem); 

local_3 = zeros(dim_elem,dim_elem); local_4 = zeros(dim_elem,dim_elem); 


% i is row which is basis of u and j is v. 1 is kappa and 2 is neighbour




for i = 1:dim_elem
   
    for j= 1:dim_elem
    
     
        % first term 1/2* (|b\cdot n| - (b\cdot n) )(u^+)*( v^+) 
         
        t1 =  0.5*(abs(sum(b_val.*n_vec,2)) -sum(b_val.*n_vec,2) )...
            .*tensor_leg3(P_Qpoints,m1,h1,Lege_ind(i,:))...
            .*tensor_leg3(P_Qpoints,m1,h1,Lege_ind(j,:));              
        
        local_1(j,i) = dot((t1),weights);
        
       
        % second term 1/2* ((b\cdot n)-|b\cdot n| )(u^-)*( v^+) 
        
        t2 =  0.5*(sum(b_val.*n_vec,2)-abs(sum(b_val.*n_vec,2)))...
            .*tensor_leg3(P_Qpoints,m2,h2,Lege_ind(i,:))...
            .*tensor_leg3(P_Qpoints,m1,h1,Lege_ind(j,:));  
            
            
        local_2(j,i) = dot((t2),weights);
        
       
        % third term -1/2* ((b\cdot n)+|b\cdot n| )(u^+)*( v^-) 
        
         t3 =  -0.5*(sum(b_val.*n_vec,2)+abs(sum(b_val.*n_vec,2)))...
            .*tensor_leg3(P_Qpoints,m1,h1,Lege_ind(i,:))...
            .*tensor_leg3(P_Qpoints,m2,h2,Lege_ind(j,:));   
            
            
        local_3(j,i) = dot((t3),weights);
        
         % fourth term 1/2* (|b\cdot n|+(b\cdot n) )(u^-)*( v^-) 
        
         t4 =  0.5*(abs(sum(b_val.*n_vec,2))+sum(b_val.*n_vec,2))...
            .*tensor_leg3(P_Qpoints,m2,h2,Lege_ind(i,:))...
            .*tensor_leg3(P_Qpoints,m2,h2,Lege_ind(j,:));   
            
            
        local_4(j,i) = dot((t4),weights);
        
    end
end


z = [local_1, local_2; local_3,local_4];

end