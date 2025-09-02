function Diri_edgeBasis = DiriBC_projecton_edge(Space_Node, Time_Node ,BDbox,S_Po,T_Po,edge_nodalbasis_index, Nodal_coe ,edge_function_index,S_Polydegree,T_Polydegree ,g_D)

%information for two bounding box , n is the normal vector from k1 to k2

h =  (BDbox(2,:)-BDbox(1,:))./2 ;

m = 0.5.*sum(BDbox);


%% generating quadrature points and weights in space

[Space_weights,Space_ref_Qpoints] = quad_GL(ceil((S_Po+1)*0.5)+1);

%change the quadrature nodes from reference domain to physical domain.

mid= sum(Space_Node)./2;   tanvec = 0.5* (Space_Node(2,:)-Space_Node(1,:));

C = kron(mid,ones(size(Space_ref_Qpoints,1),1));

Space_P_Qpoints = kron(Space_ref_Qpoints,tanvec) + C;   De = norm((Space_Node(2,:)-Space_Node(1,:))).*0.5;

Space_weights = De.*Space_weights;


%% time quadracture points

 [w_t,t_points] = quad_GL(ceil((T_Po+1)*0.5)+1);

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
    
 g_val = g_D(P_Qpoints);  Nodal_val = zeros(size(P_Qpoints,1),1);
 
 % the Nodal contribution
 
 for i =1 :size(edge_nodalbasis_index,1)

     Nodal_val = Nodal_val + Nodal_coe(i).*FEM2D_DG_basis(P_Qpoints,BDbox,m,h,edge_nodalbasis_index(i,:));
     
 end
 
 
 bd_val = g_val- Nodal_val;
 
 
 %% using the collocation projection method  based on inverting the local mass matrix 
 
 mass = NaN(size(edge_function_index,1),size(edge_function_index,1));
 
 F = NaN(size(edge_function_index,1),1);
 
 for j = 1: size(edge_function_index,1)
    
     t_F = bd_val.* tensor_leg3(P_Qpoints,m,h,edge_function_index(j,:));

     
     F(j) = dot(t_F,weights);
     
     
     for i  =1: size(edge_function_index,1)
        
         t =  tensor_leg3(P_Qpoints,m,h,edge_function_index(j,:)).*...
              FEM2D_DG_basis(P_Qpoints,BDbox,m,h,edge_function_index(i,:));
          
        mass(j,i) = dot(t,weights);
         
     end
     
          
 end
 
Diri_edgeBasis = mass\F;

end
    