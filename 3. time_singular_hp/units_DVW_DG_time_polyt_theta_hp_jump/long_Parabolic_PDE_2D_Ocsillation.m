%%%%%% The parabolic PDE, we use hp-DGFEM in time and hp-FEM in 2D space  %%%%%%%%%%%

%matlabpool open local;

clear all;

No_of_elements = 16;   

order_of_Basis = 2;

Basis_Type  = 'Se';     % Q,Se

Static_Condensation = 'on'; % on, off


for OB = 1 : length(order_of_Basis)

%% space elements

load([num2str(No_of_elements) ' rectangle Elements.mat'])


NO_space_element = NT; 


%% time elements, time grids t0,t1,t2....tn,  (x,t)

time_starting = 0;   time_end = 8;
  
% 640 time step

NO_time_step = 640;%10*time_end*sqrt(No_of_elements);

t_grids =linspace(time_starting,time_end,NO_time_step+1);


%% Define PDE struct

a = @(x)diffusive(x); %the diffusion tensor in 1D space 

b = @(x)adv(x); % the convection cooefficient

c = @(x)reac(x);% the reaction cooefficient

c0 = @(x)reac(x)-0.5.*gradient_b(x); % the auxilary function which is positive!
 
f = @(x)forcing(x); % the forcing founction in 2D space * time

u = @(x)u_true(x); % the analytic solution in 2D space * time

g_D = @(x)u_true(x); % the Dirichelet boundary condition in 2D space * time

grad_u = @(x)grad_u_true(x); % the spatial gradient of analytic solution in 1D space * time


%% Set up different boundary condition

% Neumann boudary condition is empty now! But can be defined

% Dirichele boundary condtion at the nodes

Diribd_index = unique(bdEdge(:)); Diribdnode = Node(Diribd_index,:);  


DiribdEdge=bdEdge;   Diribd_edge2elem=bd_edge2elem;  


TotalEdge = [bdEdge;intEdge]; 

%% data for FEM space

Dimension = 3; % 2D space  * 1D time

Polydegree = order_of_Basis(OB) ;  % polynomial degree of DGFEM    

Po = Polydegree*2+1;         % quadrature order 

disp([num2str(No_of_elements) ' elements with basis ' Basis_Type num2str(Polydegree) ])


%% Q basis for 2D space and 1D time

Space_FEM_index = FEM_Basis_index2D(Polydegree,Basis_Type);

time_FEM_index = 0:Polydegree;  time_FEM_index=time_FEM_index';

FEM_index = [kron(Space_FEM_index,ones(Polydegree+1,1)),...
            kron(ones(size(Space_FEM_index,1),1),time_FEM_index)];


%% iteration start k=1, end at  NO_time_step


%L_2_L_2,L_2_H_1, L_inf_L2 and L_inf_H1 norm error in different time step 

L2_L2_err_vector = NaN(NO_time_step,1);

L2_H1_err_vector = NaN(NO_time_step,1);

L_inf_L2_err_vector = NaN(NO_time_step,1);

L_inf_H1_err_vector = NaN(NO_time_step,1);



for k= 1 : NO_time_step

% 3D time space elements
    
dim_elem = size(FEM_index,1);   

% We are buliding C0 continuity space FEM space: nodal + edge mode + internal mode

NO_NodalBasis = size(Node,1);  NO_EdgeBasis = size(TotalEdge,1).*(Polydegree-1);

NO_ModalBasis = (size(Space_FEM_index,1)-4*Polydegree).*NT;


No_timeBasis = size(time_FEM_index,1);

step_space_dim_FEM = (NO_NodalBasis+ NO_EdgeBasis +NO_ModalBasis);


step_dim_FEM =step_space_dim_FEM*No_timeBasis;
    




%% Begin to calculate

% the stiffness matrix is fixed for each time!

if k==1 



%% Part 1

% contribution from element for the bilinear form

i =zeros(dim_elem.^2, NO_space_element ); j =zeros(dim_elem.^2,NO_space_element); s =zeros(dim_elem.^2,NO_space_element ); 


parfor t =1:NO_space_element
   
    % get the space 2D elements and build up to the local 3D prism elements   
        
     [elem, ~] = Elem{t,:};   Space_Node = Node(Elem{t,1},:);
    
     Time_Node = [t_grids(k),t_grids(k+1) ]';
      
    % generating the bounding box
    
     Space_BDbox = Bounding_box(Space_Node);
     
     BDbox  =  [Space_BDbox,Time_Node];
   
    diffusion_local = diffusion_localstiff(Space_Node,Time_Node,BDbox,  Po,FEM_index,@(x)a(x));
    
    CR_local = CR_localmass(Space_Node,Time_Node ,BDbox,  Po,FEM_index,@(x)b(x), @(x)c(x));
    
    
    local = diffusion_local + CR_local ;
    
    
    %% Assembling procedure  
    
    % nodal index always exists
    
    global_nodal_index = get_space_tensor_time_index(elem,No_timeBasis); 
   
   
   
   if Polydegree == 1
   
       global_ind = global_nodal_index;
       
   else
       
        % edge and  modal basis only exists for p>=2
       % find all the edge index
       
       local_edge = [elem, elem([2:end, 1])];
       
       space_EdgeBasis_ind = get_EdgeBasisindex(local_edge,TotalEdge,Polydegree);  
       
       % changle the index set from space to space-time
       
       global_EdgeBasis_ind = get_space_tensor_time_index(space_EdgeBasis_ind,No_timeBasis); 
       % find all the modal basis 
       
       NO_local_modalbasis = dim_elem - size(global_nodal_index,1) - size(global_EdgeBasis_ind,1);
       
       
       modal_index = (t-1)*NO_local_modalbasis +1  : t*NO_local_modalbasis;
       
       global_modal_index = modal_index'; 
       
       
       
       % consturct the global indices
       
       global_ind = [global_nodal_index ;...
              (global_EdgeBasis_ind + NO_NodalBasis*No_timeBasis);...
              (global_modal_index + (NO_NodalBasis+NO_EdgeBasis)*No_timeBasis)];
       
   end
                  
   
   
   i(:,t) = kron(ones(dim_elem,1),global_ind )  ;
   
   j(:,t) = kron(global_ind , ones(dim_elem,1)) ;
   
   s(:,t) = local(:);
    
end

B1 = sparse(i(:),j(:) ,s(:) ,step_dim_FEM, step_dim_FEM );

%disp('Bilinear form part 1 Element face is ready');






%% part 2 is the comtribution from union of inflow  boundary, which is the starting time for each step
 


i = zeros(dim_elem.^2,size(NO_space_element,1) ); j = zeros(dim_elem.^2,size(NO_space_element,1) );
s = zeros(dim_elem.^2,size(NO_space_element,1) );  

parfor t =1:NO_space_element
    
% get the space 2D elements and build up to the local 3D prism elements   
        
    elem_inflowface = t;

    Space_Node = Node(Elem{elem_inflowface,1},:);    

    Space_BDbox = Bounding_box(Node(Elem{elem_inflowface,1},:));            

    Space_norvec = [0,0];

% 3D bounding box and vertices
   
    BDbox = [Space_BDbox ,[t_grids(k),t_grids(k+1)]' ];    

     Time_Node = [t_grids(k),t_grids(k+1) ]';

    % % space and time are orthorgonal!! And inflow boundary should use the
    % outward normal vector [0,0,-1]!

    norvec = [Space_norvec,-1]; 

   
    CR_BDface = CR_localbdface(Space_Node,Time_Node, BDbox, norvec ,Po ,FEM_index, @(x)b(x));
    


   %% Assembling procedure  
    
    % nodal index always exists
    
    [elem, ~] = Elem{elem_inflowface,:};
    
    global_nodal_index = get_space_tensor_time_index(elem,No_timeBasis); 
   
   
   
   if Polydegree == 1
   
       global_ind = global_nodal_index;
       
   else
       
        % edge and  modal basis only exists for p>=2
       % find all the edge index
       
       local_edge = [elem, elem([2:end, 1])];
       
       space_EdgeBasis_ind = get_EdgeBasisindex(local_edge,TotalEdge,Polydegree);  
       
       % changle the index set from space to space-time
       
       global_EdgeBasis_ind = get_space_tensor_time_index(space_EdgeBasis_ind,No_timeBasis); 
       % find all the modal basis 
       
       NO_local_modalbasis = dim_elem - size(global_nodal_index,1) - size(global_EdgeBasis_ind,1);
       
       
       modal_index = (t-1)*NO_local_modalbasis +1  : t*NO_local_modalbasis;
       
       global_modal_index = modal_index'; 
       
       
       
       % consturct the global indices
       
       global_ind = [global_nodal_index ;...
              (global_EdgeBasis_ind + NO_NodalBasis*No_timeBasis);...
              (global_modal_index + (NO_NodalBasis+NO_EdgeBasis)*No_timeBasis)];
       
   end
   
   i(:,t) = kron(ones(dim_elem,1),global_ind )  ;
   
   j(:,t) = kron(global_ind , ones(dim_elem,1)) ;
   
   s(:,t) = CR_BDface(:);
    
    
end

B2 = sparse(i(:),j(:) ,s(:) ,step_dim_FEM,step_dim_FEM );

%disp('Bilinear form Part 4 union of inflow  boundary is over');



end



% the forcing term is fixed for each time




%% The right handside of the matrix, the linear functional;   

%%Part 1

% contribution from element

i =zeros(dim_elem,NO_space_element  ); j =ones(dim_elem,NO_space_element  ); s =zeros(dim_elem,NO_space_element  ); 


parfor t =1:NO_space_element
   
    % get the space 2D elements and build up to the local 3D prism elements   
        
     [elem, ~] = Elem{t,:};   Space_Node = Node(Elem{t,1},:);
    
     Time_Node = [t_grids(k),t_grids(k+1) ]';
      
    % generating the bounding box
    
     Space_BDbox = Bounding_box(Space_Node);
     
     BDbox  =  [Space_BDbox,Time_Node];
    
   
    localvec = CR_vect_forcing(Space_Node,Time_Node,BDbox,  Po, FEM_index, @(x)f(x));
    
    %% Assembling procedure  
    
    % nodal index always exists
    
    
    
    global_nodal_index = get_space_tensor_time_index(elem,No_timeBasis); 
   
   
   
   if Polydegree == 1
   
       global_ind = global_nodal_index;
       
   else
       
        % edge and  modal basis only exists for p>=2
       % find all the edge index
       
       local_edge = [elem, elem([2:end, 1])];
       
       space_EdgeBasis_ind = get_EdgeBasisindex(local_edge,TotalEdge,Polydegree);  
       
       % changle the index set from space to space-time
       
       global_EdgeBasis_ind = get_space_tensor_time_index(space_EdgeBasis_ind,No_timeBasis); 
       % find all the modal basis 
       
       NO_local_modalbasis = dim_elem - size(global_nodal_index,1) - size(global_EdgeBasis_ind,1);
       
       
       modal_index = (t-1)*NO_local_modalbasis +1  : t*NO_local_modalbasis;
       
       global_modal_index = modal_index'; 
       
       
       
       % consturct the global indices
       
       global_ind = [global_nodal_index ;...
              (global_EdgeBasis_ind + NO_NodalBasis*No_timeBasis);...
              (global_modal_index + (NO_NodalBasis+NO_EdgeBasis)*No_timeBasis)];
       
   end
    
   i(:,t) = global_ind ;
   
   s(:,t) = localvec;
    
end

L1 = sparse(i(:),j(:) ,s(:),step_dim_FEM,1 );










%% The most important part!

%Part 2

% contribution from inflow boundary over each element who share one face 
% with inflow boundary, which is the starting time for each step



i =zeros(dim_elem,NO_space_element  ); j =ones(dim_elem,NO_space_element  ); s =zeros(dim_elem,NO_space_element  ); 

% k = 1

if k==1

parfor t =1:NO_space_element

% get the space 2D elements and build up to the local 3D prism elements 
    
    elem_inflowface = t;

    Space_Node = Node(Elem{elem_inflowface,1},:);    

    Space_BDbox = Bounding_box(Node(Elem{elem_inflowface,1},:));            

    Space_norvec = [0,0];

% 3D bounding box and vertices
   
    BDbox = [Space_BDbox ,[t_grids(k),t_grids(k+1)]' ];    

    Time_Node = [t_grids(k),t_grids(k+1) ]';

    % % space and time are orthorgonal!! And inflow boundary should use the
    % outward normal vector [0,0,-1]!

    norvec = [Space_norvec,-1]; 
    
    
   inflowface = CR_vect_inflowface(Space_Node,Time_Node, BDbox, norvec ,Po ,FEM_index, @(x)b(x), @(x)g_D(x));
    
   %% Assembling procedure  
    
    % nodal index always exists
    
    [elem, ~] = Elem{elem_inflowface,:};
    
    global_nodal_index = get_space_tensor_time_index(elem,No_timeBasis); 
   
   
   
   if Polydegree == 1
   
       global_ind = global_nodal_index;
       
   else
       
        % edge and  modal basis only exists for p>=2
       % find all the edge index
       
       local_edge = [elem, elem([2:end, 1])];
       
       space_EdgeBasis_ind = get_EdgeBasisindex(local_edge,TotalEdge,Polydegree);  
       
       % changle the index set from space to space-time
       
       global_EdgeBasis_ind = get_space_tensor_time_index(space_EdgeBasis_ind,No_timeBasis); 
       % find all the modal basis 
       
       NO_local_modalbasis = dim_elem - size(global_nodal_index,1) - size(global_EdgeBasis_ind,1);
       
       
       modal_index = (t-1)*NO_local_modalbasis +1  : t*NO_local_modalbasis;
       
       global_modal_index = modal_index'; 
       
       
       
       % consturct the global indices
       
       global_ind = [global_nodal_index ;...
              (global_EdgeBasis_ind + NO_NodalBasis*No_timeBasis);...
              (global_modal_index + (NO_NodalBasis+NO_EdgeBasis)*No_timeBasis)];
       
   end 
    
   i(:,t) = global_ind;
   
   s(:,t) = inflowface;
    
end


% k>=2 we need to recall the solution of the previous time steps


else
    
    
    parfor t =1:NO_space_element

   % get the space 2D elements and build up to the local 3D prism elements 
    
    elem_inflowface = t;       

    Space_Node = Node(Elem{elem_inflowface,1},:);    

    Space_BDbox = Bounding_box(Node(Elem{elem_inflowface,1},:));            

    Space_norvec = [0,0];

    % 3D bounding box and vertices
   
    BDbox = [Space_BDbox ,[t_grids(k),t_grids(k+1)]' ];    

    previous_BDbox = [Space_BDbox ,[t_grids(k-1),t_grids(k)]' ]; 

    Time_Node = [t_grids(k),t_grids(k+1) ]';

    % % space and time are orthorgonal!! And inflow boundary should use the
    % outward normal vector [0,0,-1]!

    norvec = [Space_norvec,-1]; 
    
    
    %% Assembling procedure  
    
    % nodal index always exists
    
    [elem, ~] = Elem{elem_inflowface,:};
    
    global_nodal_index = get_space_tensor_time_index(elem,No_timeBasis); 
   
   
   
   if Polydegree == 1
   
       global_ind = global_nodal_index;
       
   else
       
        % edge and  modal basis only exists for p>=2
       % find all the edge index
       
       local_edge = [elem, elem([2:end, 1])];
       
       space_EdgeBasis_ind = get_EdgeBasisindex(local_edge,TotalEdge,Polydegree);  
       
       % changle the index set from space to space-time
       
       global_EdgeBasis_ind = get_space_tensor_time_index(space_EdgeBasis_ind,No_timeBasis); 
       % find all the modal basis 
       
       NO_local_modalbasis = dim_elem - size(global_nodal_index,1) - size(global_EdgeBasis_ind,1);
       
       
       modal_index = (t-1)*NO_local_modalbasis +1  : t*NO_local_modalbasis;
       
       global_modal_index = modal_index'; 
       
       
       
       % consturct the global indices
       
       global_ind = [global_nodal_index ;...
              (global_EdgeBasis_ind + NO_NodalBasis*No_timeBasis);...
              (global_modal_index + (NO_NodalBasis+NO_EdgeBasis)*No_timeBasis)];
       
   end 
    
    % the previous time step is used to impose inflow boundary condition       
    
    previous_coef = U_previous(global_ind); % c is coefficeint
        
    
    % use the   previous time step 
    
    %BDbox = Elem{elem_inflowface,2};   previous_BDbox = Elem{previous_elem_inflowface,2}; 
        
    
    inflowface = CR_vect_inflowface_previous_step(Space_Node,Time_Node, BDbox, previous_BDbox, previous_coef,  norvec ,Po ,FEM_index, @(x)b(x));
        
  
    
    i(:,t) = global_ind;
   
    s(:,t) = inflowface;
    
    end
    
    
end


L2 = sparse(i(:),j(:) ,s(:) ,step_dim_FEM,1 );

%disp('Convection Reaction Linear form Part 2 union of inflow and Dirichelt boundary is over');




%% The coefficient of global basis is U

U = zeros(size(L1));

%% the Dilichelet boundary is fixed in each time step

% pick up the nodal contribution, each spatial nodal basis
% will generate a set of edge fucntions (P_t+1), the value should be 
% determined by L_2 projection


i = zeros(Polydegree+1,size(Diribd_index,1) );  s = zeros(Polydegree+1,size(Diribd_index,1) ); 


parfor t =1:size(Diribd_index,1)
    
   % nodal basis is determined, then we just need to compute the projection 
   
   nodal_Diribd = Diribd_index(t); 
   
   Space_Node = Node(nodal_Diribd,:);   Time_Node = [t_grids(k),t_grids(k+1) ]';
   
   
   Diri_Nodal = DiriBC_projecton_nodal(Space_Node,Time_Node,Po , Polydegree ,@(x)g_D(x)); 
                      
   % Assembling    
     
      
   global_nodal_index = get_space_tensor_time_index(nodal_Diribd,No_timeBasis);      
    
   i(:,t) = global_nodal_index;
   
   s(:,t) = Diri_Nodal;
    
end

% put the nodal contribution into the global solution

U(i(:)) = s(:);

Nodal_solution_intex = i(:);







% pick up the edge contribution, each spatial edge
% will generate a set of edge fucntions (P_x-1)(P_t+1), the value should be 
% determined by collocation projection


if Polydegree >= 2

%% the NO of the edge function along each drection     
    

i = zeros((Polydegree-1)*No_timeBasis,size(DiribdEdge,1) );  

s = zeros((Polydegree-1)*No_timeBasis,size(DiribdEdge,1) );

parfor t =1:size(DiribdEdge,1)
    
   Dirielem_Diric = Diribd_edge2elem(t,:);   Diriedge_node = Node(DiribdEdge(t,:),:); 
   
    % nodal index  and bounding box
   
  [elem , Space_BDbox] = Elem{Dirielem_Diric,:}; 
    
  % Edge index
  
  local_edge = [elem, elem([2:end, 1])];
       
   
  EdgeBasis_ind = get_EdgeBasisindex(local_edge,TotalEdge,Polydegree);
  
    
  % local index for nodal and edge;  
    
  local_ind = [elem ; EdgeBasis_ind + NO_NodalBasis ];  
  
  
  global_ind = get_space_tensor_time_index(local_ind,No_timeBasis);   
  
  % find the correct index for nodal basis and pick up the correc basis for
  % computing the projection
  
  [Dirinodal_id, Diriedge_id] = findnodal_edge_time_index(Diriedge_node,Space_BDbox,FEM_index,No_timeBasis);
   
   edge_nodalbasis_index = FEM_index(Dirinodal_id,:);
   
   edge_function_index = FEM_index(Diriedge_id,:);
      
   %% calculating the projection of boundary condtion
   
   gloabl_Nodal_ind = global_ind(Dirinodal_id);
   
   % The coefficients for the nodal basis are calculated from the previous
   % step
   
   Nodal_coe = U(gloabl_Nodal_ind);     Time_Node = [t_grids(k),t_grids(k+1) ]';
   
   BDbox = [Space_BDbox,Time_Node];
      
   Diri_Modal = DiriBC_projecton_edge(Diriedge_node,Time_Node ,BDbox,Po ,edge_nodalbasis_index , Nodal_coe ,edge_function_index, Polydegree ,@(x)g_D(x));
   
   % Assembling
   
   gloabl_edge_ind = global_ind(Diriedge_id);      
    
   i(:,t) = gloabl_edge_ind;
   
   s(:,t) = Diri_Modal;
    
end

U(i(:)) = s(:);

Edge_solution_intex = i(:);

totalbd_index = [Nodal_solution_intex ; Edge_solution_intex];

end





%% Solving the linear system


B = B1-B2; 

L = L1 -L2;



% The coefficient of global basis is U

% all the nodes except the Dirichele boundary Nodes

int_index = 1:step_dim_FEM;   int_index = int_index';

if  Polydegree == 1

int_index(Nodal_solution_intex) = [];

else
   
int_index(totalbd_index) = [];    
    
end


L = L - B*U;


switch Static_Condensation

    case 'off'   
    

U(int_index) = B(int_index,int_index)\L(int_index);


  case 'on'
             
U(int_index) = StaticCondensation_FEM(B(int_index,int_index),L(int_index),NT,(size(Space_FEM_index,1)-4*Polydegree)*No_timeBasis);        

end


U_previous = U;



%% L_2_L_2, L_2_H_1, L_inf_L2 and L_inf_H1 norm 


L2_L2_err = NaN(NO_space_element,1); 

L2_H1_err = NaN(NO_space_element,1); 


L_inf_L2_err = NaN(Po,NO_space_element); 

L_inf_H1_err = NaN(Po,NO_space_element);
  
parfor t =1:NO_space_element
   
     % get the space 2D elements and build up to the local 3D prism elements   
        
     [elem, ~] = Elem{t,:};   Space_Node = Node(Elem{t,1},:);
    
     Time_Node = [t_grids(k),t_grids(k+1) ]';
      
    % generating the bounding box
    
     Space_BDbox = Bounding_box(Space_Node);
     
     BDbox  =  [Space_BDbox,Time_Node];
     
     
     %% Assembling procedure  
    
    % nodal index always exists
       
    
    global_nodal_index = get_space_tensor_time_index(elem,No_timeBasis); 
   
   
   
   if Polydegree == 1
   
       global_ind = global_nodal_index;
       
   else
       
        % edge and  modal basis only exists for p>=2
       % find all the edge index
       
       local_edge = [elem, elem([2:end, 1])];
       
       space_EdgeBasis_ind = get_EdgeBasisindex(local_edge,TotalEdge,Polydegree);  
       
       % changle the index set from space to space-time
       
       global_EdgeBasis_ind = get_space_tensor_time_index(space_EdgeBasis_ind,No_timeBasis); 
       % find all the modal basis 
       
       NO_local_modalbasis = dim_elem - size(global_nodal_index,1) - size(global_EdgeBasis_ind,1);
       
       
       modal_index = (t-1)*NO_local_modalbasis +1  : t*NO_local_modalbasis;
       
       global_modal_index = modal_index'; 
       
       
       
       % consturct the global indices
       
       global_ind = [global_nodal_index ;...
              (global_EdgeBasis_ind + NO_NodalBasis*No_timeBasis);...
              (global_modal_index + (NO_NodalBasis+NO_EdgeBasis)*No_timeBasis)];
       
   end 
    


     coef =  U(global_ind ,1); % c is coefficeint
   
   [L_2_L2_part, L_2_H1_part] = Err_L2_H1_norm(Space_Node,Time_Node,BDbox,...
                                     coef, Po ,FEM_index,@(x)u(x) ,@(x)grad_u(x) ,@(x)a(x));
         
    
     L2_L2_err(t) = L_2_L2_part;                             
    
    L2_H1_err(t) = L_2_H1_part;
    
    
    
    [L_inf_L2_part, L_inf_H1_part ] = Err_L_inf_norm(Space_Node,Time_Node,BDbox ,...
                                     coef ,Po ,FEM_index,@(x)u(x) ,@(x)grad_u(x) ,@(x)a(x));

    L_inf_L2_err(:,t) = L_inf_L2_part; 
    
    L_inf_H1_err(:,t) = L_inf_H1_part; 
    
   
   
end

% L_2 type of norm sum over all the spatial and time element

L2_L2_err = sqrt(sum(L2_L2_err));  
 

L2_H1_err = sqrt(sum(L2_H1_err)); 

 % L_inf type norm over each time points ,L_2 type of norm sum over all the spatial element


 
L_inf_L2_err= norm(sqrt(sum(L_inf_L2_err,2)),inf);

L_inf_H1_err= norm(sqrt(sum(L_inf_H1_err,2)),inf); 





L2_L2_err_vector(k) =  L2_L2_err;

L2_H1_err_vector(k) =  L2_H1_err;

L_inf_L2_err_vector(k) =  L_inf_L2_err;

L_inf_H1_err_vector(k) =  L_inf_H1_err;



end




% calculate all the dof based on the information about dof in each time
% step


dim_FEM = step_dim_FEM.*NO_time_step;

L2_L2_err = norm(L2_L2_err_vector);  % L2 type norm

L2_H1_err = norm(L2_H1_err_vector);  % L2 type norm

% L_inf type norm take maximum from different sampling points on T

L_inf_L2_err = norm(L_inf_L2_err_vector,inf);   

L_inf_H1_err = norm(L_inf_H1_err_vector,inf);  



disp(['the L2_L2 norm error is ' num2str(L2_L2_err)]);

disp(['the L2_H1 norm error is ' num2str(L2_H1_err)]);

disp(['the L_inf_L2 norm error is ' num2str(L_inf_L2_err)]);

disp(['the L_inf_H1 norm error is ' num2str(L_inf_H1_err)]);



% save the result

savefile = ['Error ' num2str(NT) ' rectangle Elements time step ' num2str(NO_time_step) ' for FEM ' Basis_Type num2str(Polydegree) ' basis.mat'];

save(savefile, 'L2_L2_err','L2_H1_err','L_inf_L2_err','L_inf_H1_err','L2_H1_err_vector','L2_L2_err_vector','L_inf_H1_err_vector','L_inf_L2_err_vector','dim_FEM','Polydegree','NO_time_step','-v7.3');

end
%matlabpool close;
