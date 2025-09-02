function local = CR_vect_forcing(Space_Node,Time_Node,BDbox,S_Po,T_Po,FEM_index, f)


% information about the bounding box

h = (BDbox(2,:)-BDbox(1,:))./2;  

m = 0.5.*sum(BDbox);

dim_elem =size(FEM_index,1); % number of basis for each element.

local = zeros(dim_elem,1);  % initialize the localstiffness matrix



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
    
    f_val = f(P_Qpoints); 

first = zeros(dim_elem,1);
    
for j = 1:dim_elem
   

            
        % first term {fv} 
         
        t = f_val.*vt_FEM2D_DG_basis(P_Qpoints,BDbox,m,h,FEM_index(j,:));
        
        first(j) = dot((t),weights);
        
end

local = local + first;




