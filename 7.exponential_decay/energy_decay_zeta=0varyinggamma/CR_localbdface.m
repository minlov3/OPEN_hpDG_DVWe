function local = CR_localbdface(Space_Node,Time_Node,BDbox, n ,S_Po,T_Po,FEM_index,a,zeta, b)


% information about the bounding box

h = (BDbox(2,:)-BDbox(1,:))./2;  

m = 0.5.*sum(BDbox);

dim_elem =size(FEM_index,1); % number of basis for each element.



%Quadrature rules over rectangle

[Space_P_Qpoints, Space_weights] = quad_rect(Space_Node,ceil((S_Po+1)*0.5));




%% no need to put time quadracture points, t = Time_Node(1)

weights = Space_weights;


 P_Qpoints = [Space_P_Qpoints,Time_Node(1).*ones(size(Space_P_Qpoints,1),1) ];


 % data for quadrature
    
a_val = a(P_Qpoints); 
zeta_val = zeta(P_Qpoints);
b_val = b(P_Qpoints); n_vec =kron(n,ones(size(P_Qpoints,1),1));
   
 
first = zeros(dim_elem,dim_elem); 
second = zeros(dim_elem,dim_elem); 


% calculated



for i = 1:dim_elem
   
    for j= i:dim_elem
    
      % first term b\cdot n uv  (Symmetric)
         
        t_1 =  sum(b_val.*n_vec,2)...
                .*(vt_FEM2D_DG_basis(P_Qpoints,BDbox,m,h,FEM_index(i,:)))...
                .*(vt_FEM2D_DG_basis(P_Qpoints,BDbox,m,h,FEM_index(j,:)));
        
        first(j,i) = dot((t_1),weights);
      
        % second term a grad(u) grad(v)  (Symmetric)
        grad1 = grad_FEM2D_DG_basis(P_Qpoints,BDbox,m,h,FEM_index(i,:));
        grad2 = grad_FEM2D_DG_basis(P_Qpoints,BDbox,m,h,FEM_index(j,:));
        
        gradu_gradv= [grad1(:,1).*grad2(:,1) , grad1(:,1).* grad2(:,2), grad1(:,1).* grad2(:,3)...
                      grad1(:,2).*grad2(:,1) , grad1(:,2).* grad2(:,2), grad1(:,2).* grad2(:,3)...
                      grad1(:,3).*grad2(:,1) , grad1(:,3).* grad2(:,2), grad1(:,3).* grad2(:,3)];       
        
        t_2 =  sum((gradu_gradv).*zeta_val,2); 
        
        second(j,i) = dot((t_2),weights);
        
    end
end

local = tril(first,-1)'+tril(first)  +  tril(second,-1)'+tril(second);




end