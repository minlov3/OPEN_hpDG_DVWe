

figure;


%% plotting errors for the DG for p refinement

Polynomial_degree = 1:6;   NO_elem = 64;  

penalty=10;  

Basis_Type  = 'Q';  % Q,Se

Norm_type = 'L_inf_L2';  %L2_L2, L2_H1, L_inf_L2, L_inf_H1

Mesh_type = 'rect';   % rect, 


errL_2 =NaN(length(Polynomial_degree),1); dof = NaN(length(Polynomial_degree),1);

time_end = 40;

No_time_step = 3200;

switch  Norm_type

    
case 'L2_H1'


for i=1 :length(Polynomial_degree)
    

   

load(['Error ' num2str(NO_elem) ' rectangle Elements time step ' num2str(No_time_step) ' for FEM '  Basis_Type num2str(Polynomial_degree(i)) ' basis.mat'])



errL_2(i)=L2_H1_err ;         dof(i) = dim_FEM;
  
end



%semilogy(Polynomial_degree,errL_2,'-h','LineWidth',1,'MarkerSize',8);

switch Basis_Type 


case 'Q'

semilogy(dof.^(1/3),errL_2,'-mo','LineWidth',1,'MarkerSize',8);

case 'Se'

semilogy(dof.^(1/3),errL_2,'-mp','LineWidth',1,'MarkerSize',8);


end



legend([num2str(NO_elem) ' rect FEM(' Basis_Type ') ' num2str(No_time_step) ' time steps'])




%xlabel('Polynomial order of basis','FontSize',18);

xlabel('Dof^{1/3}','FontSize',18);

ylabel('|| u-u_{h}||_{L^2((0,T),H^1(\Omega))}','FontSize',18);

title('p refinement','FontSize',18); set(gca,'FontSize',15);






    
case 'L2_L2'


for i=1 :length(Polynomial_degree)
    

   

load(['Error ' num2str(NO_elem) ' rectangle Elements time step ' num2str(No_time_step) ' for FEM '  Basis_Type num2str(Polynomial_degree(i)) ' basis.mat'])



errL_2(i)=L2_H1_err ;         dof(i) = dim_FEM;
  
end



%semilogy(Polynomial_degree,errL_2,'-h','LineWidth',1,'MarkerSize',8);

switch Basis_Type 


case 'Q'

semilogy(dof.^(1/3),errL_2,'-mo','LineWidth',1,'MarkerSize',8);

case 'Se'

semilogy(dof.^(1/3),errL_2,'-mp','LineWidth',1,'MarkerSize',8);


end



legend([num2str(NO_elem) ' rect FEM(' Basis_Type ') ' num2str(No_time_step) ' time steps'])




%xlabel('Polynomial order of basis','FontSize',18);

xlabel('Dof^{1/3}','FontSize',18);

ylabel('|| u-u_{h}||_{L^2((0,T),L^2(\Omega))}','FontSize',18);

title('p refinement','FontSize',18); set(gca,'FontSize',15);






case 'L_inf_L2'


for i=1 :length(Polynomial_degree)
    


load(['Error ' num2str(NO_elem) ' rectangle Elements time step ' num2str(No_time_step) ' for FEM '  Basis_Type num2str(Polynomial_degree(i)) ' basis.mat'])




errL_2(i)=L_inf_L2_err ;         dof(i) = dim_FEM;
  
end



%semilogy(Polynomial_degree,errL_2,'-h','LineWidth',1,'MarkerSize',8);

switch Basis_Type 

case 'Q'

semilogy(dof.^(1/3),errL_2,'-mo','LineWidth',1,'MarkerSize',8);

case 'Se'

semilogy(dof.^(1/3),errL_2,'-mp','LineWidth',1,'MarkerSize',8);


end


legend([num2str(NO_elem) ' rect FEM(' Basis_Type ') ' num2str(No_time_step) ' time steps'])




%xlabel('Polynomial order of basis','FontSize',18);

xlabel('Dof^{1/3}','FontSize',18);

ylabel('|| u-u_{h}||_{L^\infty((0,T),L^2(\Omega))}','FontSize',18);

title('p refinement','FontSize',18); set(gca,'FontSize',15);



case 'L_inf_H1'


for i=1 :length(Polynomial_degree)
    

load(['Error ' num2str(NO_elem) ' rectangle Elements time step ' num2str(No_time_step) ' for FEM '  Basis_Type num2str(Polynomial_degree(i)) ' basis.mat'])


 

errL_2(i)=L_inf_H1_err ;         dof(i) = dim_FEM;
  
end



switch Basis_Type 


case 'Q'

semilogy(dof.^(1/3),errL_2,'-mo','LineWidth',1,'MarkerSize',8);

case 'Se'

semilogy(dof.^(1/3),errL_2,'-mp','LineWidth',1,'MarkerSize',8);


end


legend([num2str(NO_elem) ' rect FEM(' Basis_Type ') ' num2str(No_time_step) ' time steps'])


%xlabel('Polynomial order of basis','FontSize',18);

xlabel('Dof^{1/3}','FontSize',18);

ylabel('|| u-u_{h}||_{L^\infty((0,T),H^1(\Omega))}','FontSize',18);

title('p refinement','FontSize',18); set(gca,'FontSize',15);


end