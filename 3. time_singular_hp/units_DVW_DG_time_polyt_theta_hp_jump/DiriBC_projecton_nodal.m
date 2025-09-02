function Diri_Nodal =  DiriBC_projecton_nodal(Space_Node,Time_Node,S_Po,T_Po ,S_Polydegree,T_Polydegree,g_D)

% at any fixed aspce points, we only have 1D basis in time!

h =  (Time_Node(2,:)-Time_Node(1,:))./2 ;

m = 0.5.*sum(Time_Node);

time_FEM_index = 0:T_Polydegree; 

time_FEM_index = time_FEM_index';

%% time quadracture points

 [w_t,t_points] = quad_GL(ceil((T_Po+1)*0.5)+1);

% map the reference interval [-1,1] to [t_k, t_(k+1)]

Time_weights = w_t.*(Time_Node(2)-Time_Node(1))*0.5;

Time_P_Qpoints = t_points.*(Time_Node(2)-Time_Node(1))*0.5+ 0.5*sum(Time_Node);

% change 1D interval to 3D 

P_Qpoints = [kron(Space_Node,ones(size(Time_weights,1),1)),Time_P_Qpoints ];


weights = Time_weights ;




% data for quadrature
    
 g_val = g_D(P_Qpoints);

 
 
 %% using L2 projection method  based on inverting the local mass matrix 
 %% Not the L2 projection!!!!! RGH GAUSS-RADAU
 
 mass = NaN(T_Polydegree+1,T_Polydegree+1); F = NaN(T_Polydegree+1,1);
 
 for j = 1: T_Polydegree+1
     
     t_F = g_val.*shift_leg_derivative(Time_P_Qpoints,m,h,time_FEM_index(j),0);

     
     F(j) = dot(t_F,weights);
     
     for i  =1: T_Polydegree+1
                 
          
          
         t = shift_leg_derivative(Time_P_Qpoints,m,h,time_FEM_index(j),0).*...
             shift_leg_derivative(Time_P_Qpoints,m,h,time_FEM_index(i),0);
          
          
        mass(j,i) = dot(t,weights);
         
     end
     
 end

 Diri_Nodal = mass\F;
 
 
 
 
 %{
%% Computing the Gauss-Radua Interpolant, Which is L2 projection on the first p basis 
 
Diri_Nodal = NaN(Polydegree+1,1); value_t_end = NaN(Polydegree,1);

for i = 1 :Polydegree

     t_F = g_val.*shift_leg_derivative(Time_P_Qpoints,m,h,time_FEM_index(i),0);
    
     coe = dot(t_F,weights);
     
     % normalization
     
     Temp = shift_leg_derivative(Time_P_Qpoints,m,h,time_FEM_index(i),0).^2;
     
     nor = dot(Temp,weights);
     
    Diri_Nodal(i) = coe/nor;
    
    value_t_end(i) = shift_leg_derivative(Time_Node(2,:),m,h,time_FEM_index(i),0);
    
end

%% the last basis is to interpolate on right handside 


final_value = shift_leg_derivative(Time_Node(2,:),m,h,time_FEM_index(Polydegree+1),0);
     

Diri_Nodal(Polydegree+1) =  (g_D([Space_Node,Time_Node(2,:)])  - dot(value_t_end,Diri_Nodal(1:Polydegree)))/final_value ;

 
end
                  
%}

%shift_leg_derivative(t,m,h,order_basis,0);