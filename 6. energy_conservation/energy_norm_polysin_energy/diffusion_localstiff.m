function local = diffusion_localstiff(Space_Node,Time_Node,BDbox,S_Po,T_Po,FEM_index, a,zeta,eta,gamma)


% information about the bounding box

h = (BDbox(2,:)-BDbox(1,:))./2;  

m = 0.5.*sum(BDbox);  

dim_elem =size(FEM_index,1); % number of basis for each element.

local = zeros(dim_elem,dim_elem);  % initialize the localstiffness matrix


%Quadrature rules over rectangle

[Space_P_Qpoints, Space_weights] = quad_rect(Space_Node,ceil((S_Po+1)*0.5));



%% time quadracture points

 [w_t,t_points] = quad_GL(ceil((T_Po+1)*0.5));

% map the reference interval [-1,1] to [t_k, t_(k+1)]

Time_weights = w_t.*(Time_Node(2)-Time_Node(1))*0.5;

Time_P_Qpoints = t_points.*(Time_Node(2)-Time_Node(1))*0.5+ 0.5*sum(Time_Node);

%% generator the 3D quad points


quad_x = kron(Space_P_Qpoints(:,1),ones(size(Time_weights,1),1));

quad_y = kron(Space_P_Qpoints(:,2),ones(size(Time_weights,1),1));

quad_t = kron(ones(size(Space_P_Qpoints,1),1),Time_P_Qpoints); 

weights = kron(Space_weights,Time_weights);


P_Qpoints = [quad_x,quad_y,quad_t];


      
     % data for quadrature
    
    a_val = a(P_Qpoints);    
    zeta_val = zeta(P_Qpoints);
    eta_val = eta(P_Qpoints);
    gamma_val = gamma(P_Qpoints);
    
first = zeros(dim_elem,dim_elem); 
second = zeros(dim_elem,dim_elem); 
gam = zeros(dim_elem,dim_elem); 

for i = 1:dim_elem
    
    
    %%symetric term
   
    for j = 1:dim_elem
    
        % first term {a grad(u) \cdot grad(v)} is symetric, u is first 
        % component and v is 2 second component
             
        grad11 = grad_FEM2D_DG_basis(P_Qpoints,BDbox,m,h,FEM_index(i,:));
        grad12 = grad_FEM2D_DG_basis(P_Qpoints,BDbox,m,h,FEM_index(j,:));
        
        grad21 = vt_grad_FEM2D_DG_basis(P_Qpoints,BDbox,m,h,FEM_index(i,:));
        grad22 = vt_grad_FEM2D_DG_basis(P_Qpoints,BDbox,m,h,FEM_index(j,:));
        
        % % utvt
        gradu_gradv = [grad11(:,1).*grad12(:,1) , grad11(:,1).* grad12(:,2), grad11(:,1).* grad12(:,3)...
                       grad11(:,2).*grad12(:,1) , grad11(:,2).* grad12(:,2), grad11(:,2).* grad12(:,3)...
                       grad11(:,3).*grad12(:,1) , grad11(:,3).* grad12(:,2), grad11(:,3).* grad12(:,3)]; 
        
        t0 = sum((gradu_gradv).*gamma_val,2);
        gam(j,i) = dot((t0),weights);           
                    
         % %grad u grad vt   
         gradu_gradvt = [grad11(:,1).*grad22(:,1) , grad11(:,1).* grad22(:,2), grad11(:,1).* grad22(:,3)...
                         grad11(:,2).*grad22(:,1) , grad11(:,2).* grad22(:,2), grad11(:,2).* grad22(:,3)...
                         grad11(:,3).*grad22(:,1) , grad11(:,3).* grad22(:,2), grad11(:,3).* grad22(:,3)];          
                
        t1 = sum((gradu_gradvt).*zeta_val,2);
        first(j,i) = dot((t1),weights);
             
        gradut_gradvt = [grad21(:,1).*grad22(:,1) , grad21(:,1).* grad22(:,2), grad21(:,1).* grad22(:,3)...
                         grad21(:,2).*grad22(:,1) , grad21(:,2).* grad22(:,2), grad21(:,2).* grad22(:,3)...
                         grad21(:,3).*grad22(:,1) , grad21(:,3).* grad22(:,2), grad21(:,3).* grad22(:,3)];   
        
        t2 = sum((gradut_gradvt).*eta_val,2);
        second(j,i) = dot((t2),weights);
        
        
    end
end


%%symetric term

local = local + first + second + gam;


end

