figure;

%% plotting errors for the DG for p refinement

S_Polydegree = 3; 

T_Polydegree_vec =  1:6;

NO_elem = 16;  
NO_time_step = 40;

Basis_Type  = 'Q';  % Q, Se

Norm_type = 'L2_L2'; 
%  Norm_type = 'L2_H1'; 
%  Norm_type = 'L_inf_L2';  
% Norm_type = 'L_inf_H1';
% Norm_type = 'H1_L2'; 
% Norm_type = 'H1_H1'; 
% Norm_type = 'W_1_inf_L2';  
% Norm_type = 'W_1_inf_H1';
%  Norm_type = 'Final_t_L2_L2'; 
% Norm_type = 'Final_t_L2_H1'; 
%  Norm_type = 'Final_t_H1_L2'; 
%  Norm_type = 'Final_t_H1_H1'; 


Mesh_type = 'rect';   % rect, 

errL_2 =NaN(length(T_Polydegree_vec),1); 
dof = NaN(length(T_Polydegree_vec),1);


switch  Norm_type

case 'L2_L2'
for i=1 :length(T_Polydegree_vec) 
load(['Error of ' num2str(NO_elem) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree  ) ' and time basis ' num2str(T_Polydegree_vec(i)) ' basis.mat'])

errL_2(i)=L2_L2_err ;         
dof(i) = dim_FEM;
end
switch Basis_Type 
case 'Q'
semilogy(dof.^(1/2),errL_2,'-mo','LineWidth',1,'MarkerSize',8);
case 'Se'
semilogy(dof.^(1/2),errL_2,'-go','LineWidth',1,'MarkerSize',8);
end
legend([num2str(NO_elem) ' rect FEM(' Basis_Type ') ' num2str(NO_time_step) ' time steps'])
%xlabel('Polynomial order of basis','FontSize',18);
xlabel('Dof^{1/2}','FontSize',18);
ylabel(  'Norm_type error ','FontSize',18);
title('time p refinement','FontSize',18); set(gca,'FontSize',15);    
    
case 'L2_H1'
for i=1 :length(T_Polydegree_vec) 
load(['Error of ' num2str(NO_elem) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree  ) ' and time basis ' num2str(T_Polydegree_vec(i)) ' basis.mat'])

errL_2(i)=L2_H1_err ;         
dof(i) = dim_FEM;
end
switch Basis_Type 
case 'Q'
semilogy(dof.^(1/2),errL_2,'-mo','LineWidth',1,'MarkerSize',8);
case 'Se'
semilogy(dof.^(1/2),errL_2,'-go','LineWidth',1,'MarkerSize',8);
end
legend([num2str(NO_elem) ' rect FEM(' Basis_Type ') ' num2str(NO_time_step) ' time steps'])
%xlabel('Polynomial order of basis','FontSize',18);
xlabel('Dof^{1/2}','FontSize',18);
ylabel(  'Norm_type error ','FontSize',18);
title('time p refinement','FontSize',18); set(gca,'FontSize',15);

case 'L_inf_L2'
for i=1 :length(T_Polydegree_vec) 
load(['Error of ' num2str(NO_elem) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree  ) ' and time basis ' num2str(T_Polydegree_vec(i)) ' basis.mat'])

errL_2(i)=L_inf_L2_err ;         
dof(i) = dim_FEM;
end
switch Basis_Type 
case 'Q'
semilogy(dof.^(1/2),errL_2,'-mo','LineWidth',1,'MarkerSize',8);
case 'Se'
semilogy(dof.^(1/2),errL_2,'-go','LineWidth',1,'MarkerSize',8);
end
legend([num2str(NO_elem) ' rect FEM(' Basis_Type ') ' num2str(NO_time_step) ' time steps'])
%xlabel('Polynomial order of basis','FontSize',18);
xlabel('Dof^{1/2}','FontSize',18);
ylabel(  'Norm_type error ','FontSize',18);
title('time p refinement','FontSize',18); set(gca,'FontSize',15);

case 'L_inf_H1'
for i=1 :length(T_Polydegree_vec) 
load(['Error of ' num2str(NO_elem) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree  ) ' and time basis ' num2str(T_Polydegree_vec(i)) ' basis.mat'])

errL_2(i)=L_inf_H1_err ;         
dof(i) = dim_FEM;
end
switch Basis_Type 
case 'Q'
semilogy(dof.^(1/2),errL_2,'-mo','LineWidth',1,'MarkerSize',8);
case 'Se'
semilogy(dof.^(1/2),errL_2,'-go','LineWidth',1,'MarkerSize',8);
end
legend([num2str(NO_elem) ' rect FEM(' Basis_Type ') ' num2str(NO_time_step) ' time steps'])
%xlabel('Polynomial order of basis','FontSize',18);
xlabel('Dof^{1/2}','FontSize',18);
ylabel(  'Norm_type error ','FontSize',18);
title('time p refinement','FontSize',18); set(gca,'FontSize',15);

case 'H1_L2'
for i=1 :length(T_Polydegree_vec) 
load(['Error of ' num2str(NO_elem) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree  ) ' and time basis ' num2str(T_Polydegree_vec(i)) ' basis.mat'])

errL_2(i)=H1_L2_err ;         
dof(i) = dim_FEM;
end
switch Basis_Type 
case 'Q'
semilogy(dof.^(1/2),errL_2,'-mo','LineWidth',1,'MarkerSize',8);
case 'Se'
semilogy(dof.^(1/2),errL_2,'-go','LineWidth',1,'MarkerSize',8);
end
legend([num2str(NO_elem) ' rect FEM(' Basis_Type ') ' num2str(NO_time_step) ' time steps'])
%xlabel('Polynomial order of basis','FontSize',18);
xlabel('Dof^{1/2}','FontSize',18);
ylabel(  'Norm_type error ','FontSize',18);
title('time p refinement','FontSize',18); set(gca,'FontSize',15);

case 'H1_H1'
for i=1 :length(T_Polydegree_vec) 
load(['Error of ' num2str(NO_elem) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree  ) ' and time basis ' num2str(T_Polydegree_vec(i)) ' basis.mat'])

errL_2(i)=H1_H1_err ;         
dof(i) = dim_FEM;
end
switch Basis_Type 
case 'Q'
semilogy(dof.^(1/2),errL_2,'-mo','LineWidth',1,'MarkerSize',8);
case 'Se'
semilogy(dof.^(1/2),errL_2,'-go','LineWidth',1,'MarkerSize',8);
end
legend([num2str(NO_elem) ' rect FEM(' Basis_Type ') ' num2str(NO_time_step) ' time steps'])
%xlabel('Polynomial order of basis','FontSize',18);
xlabel('Dof^{1/2}','FontSize',18);
ylabel(  'Norm_type error ','FontSize',18);
title('time p refinement','FontSize',18); set(gca,'FontSize',15);

case 'W_1_inf_L2'
for i=1 :length(T_Polydegree_vec) 
load(['Error of ' num2str(NO_elem) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree  ) ' and time basis ' num2str(T_Polydegree_vec(i)) ' basis.mat'])

errL_2(i)=W_1_inf_L2_err ;         
dof(i) = dim_FEM;
end
switch Basis_Type 
case 'Q'
semilogy(dof.^(1/2),errL_2,'-mo','LineWidth',1,'MarkerSize',8);
case 'Se'
semilogy(dof.^(1/2),errL_2,'-go','LineWidth',1,'MarkerSize',8);
end
legend([num2str(NO_elem) ' rect FEM(' Basis_Type ') ' num2str(NO_time_step) ' time steps'])
%xlabel('Polynomial order of basis','FontSize',18);
xlabel('Dof^{1/2}','FontSize',18);
ylabel(  'Norm_type error ','FontSize',18);
title('time p refinement','FontSize',18); set(gca,'FontSize',15);

case 'W_1_inf_H1'
for i=1 :length(T_Polydegree_vec) 
load(['Error of ' num2str(NO_elem) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree  ) ' and time basis ' num2str(T_Polydegree_vec(i)) ' basis.mat'])

errL_2(i)=W_1_inf_H1_err ;         
dof(i) = dim_FEM;
end
switch Basis_Type 
case 'Q'
semilogy(dof.^(1/2),errL_2,'-mo','LineWidth',1,'MarkerSize',8);
case 'Se'
semilogy(dof.^(1/2),errL_2,'-go','LineWidth',1,'MarkerSize',8);
end
legend([num2str(NO_elem) ' rect FEM(' Basis_Type ') ' num2str(NO_time_step) ' time steps'])
%xlabel('Polynomial order of basis','FontSize',18);
xlabel('Dof^{1/2}','FontSize',18);
ylabel(  'Norm_type error ','FontSize',18);
title('time p refinement','FontSize',18); set(gca,'FontSize',15);

case 'Final_t_L2_L2'
for i=1 :length(T_Polydegree_vec) 
load(['Error of ' num2str(NO_elem) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree  ) ' and time basis ' num2str(T_Polydegree_vec(i)) ' basis.mat'])

errL_2(i)=Final_t_L2_L2_err ;         
dof(i) = dim_FEM;
end
switch Basis_Type 
case 'Q'
semilogy(dof.^(1/2),errL_2,'-mo','LineWidth',1,'MarkerSize',8);
case 'Se'
semilogy(dof.^(1/2),errL_2,'-go','LineWidth',1,'MarkerSize',8);
end
legend([num2str(NO_elem) ' rect FEM(' Basis_Type ') ' num2str(NO_time_step) ' time steps'])
%xlabel('Polynomial order of basis','FontSize',18);
xlabel('Dof^{1/2}','FontSize',18);
ylabel(  'Norm_type error ','FontSize',18);
title('time p refinement','FontSize',18); set(gca,'FontSize',15);

case 'Final_t_L2_H1'
for i=1 :length(T_Polydegree_vec) 
load(['Error of ' num2str(NO_elem) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree  ) ' and time basis ' num2str(T_Polydegree_vec(i)) ' basis.mat'])

errL_2(i)=Final_t_L2_H1_err ;         
dof(i) = dim_FEM;
end
switch Basis_Type 
case 'Q'
semilogy(dof.^(1/2),errL_2,'-mo','LineWidth',1,'MarkerSize',8);
case 'Se'
semilogy(dof.^(1/2),errL_2,'-go','LineWidth',1,'MarkerSize',8);
end
legend([num2str(NO_elem) ' rect FEM(' Basis_Type ') ' num2str(NO_time_step) ' time steps'])
%xlabel('Polynomial order of basis','FontSize',18);
xlabel('Dof^{1/2}','FontSize',18);
ylabel(  'Norm_type error ','FontSize',18);
title('time p refinement','FontSize',18); set(gca,'FontSize',15);

case 'Final_t_H1_L2'
for i=1 :length(T_Polydegree_vec) 
load(['Error of ' num2str(NO_elem) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree  ) ' and time basis ' num2str(T_Polydegree_vec(i)) ' basis.mat'])

errL_2(i)=Final_t_H1_L2_err ;         
dof(i) = dim_FEM;
end
switch Basis_Type 
case 'Q'
semilogy(dof.^(1/2),errL_2,'-mo','LineWidth',1,'MarkerSize',8);
case 'Se'
semilogy(dof.^(1/2),errL_2,'-go','LineWidth',1,'MarkerSize',8);
end
legend([num2str(NO_elem) ' rect FEM(' Basis_Type ') ' num2str(NO_time_step) ' time steps'])
%xlabel('Polynomial order of basis','FontSize',18);
xlabel('Dof^{1/2}','FontSize',18);
ylabel(  'Norm_type error ','FontSize',18);
title('time p refinement','FontSize',18); set(gca,'FontSize',15);

case 'Final_t_H1_H1'
for i=1 :length(T_Polydegree_vec) 
load(['Error of ' num2str(NO_elem) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree  ) ' and time basis ' num2str(T_Polydegree_vec(i)) ' basis.mat'])

errL_2(i)=Final_t_H1_H1_err ;         
dof(i) = dim_FEM;
end
switch Basis_Type 
case 'Q'
semilogy(dof.^(1/2),errL_2,'-mo','LineWidth',1,'MarkerSize',8);
case 'Se'
semilogy(dof.^(1/2),errL_2,'-go','LineWidth',1,'MarkerSize',8);
end
legend([num2str(NO_elem) ' rect FEM(' Basis_Type ') ' num2str(NO_time_step) ' time steps'])
%xlabel('Polynomial order of basis','FontSize',18);
xlabel('Dof^{1/2}','FontSize',18);
ylabel(  'Norm_type error ','FontSize',18);
title('time p refinement','FontSize',18); set(gca,'FontSize',15);

end