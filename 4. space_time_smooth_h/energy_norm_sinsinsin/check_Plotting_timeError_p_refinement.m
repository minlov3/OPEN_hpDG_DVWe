clear all; close all;

%% different elements we use different order of basis
Mesh_type = 'rect';
NO_elem  = 16; 
S_Polydegree = 2;
Basis_Type  = 'Q';  % Q


color = {'b','r','g','k','y','c'};
NO_time_step_vec = [10,20];
T_Polydegree_vec = [1,2,3,4,5];
a=1; %%%%% Dof^{a}


 Norm_type = 'L2_L2';   
%  Norm_type = 'L2_grad';
%  Norm_type = 'L2_H1';
%  Norm_type = 'L_inf_L2';
%  Norm_type = 'L_inf_grad';
%  Norm_type = 'L_inf_H1';
%  Norm_type = 'H1_L2'; 
%  Norm_type = 'H1_grad';
%  Norm_type = 'H1_H1';
%  Norm_type = 'W_1_inf_L2';
%  Norm_type = 'W_1_inf_grad';
%  Norm_type = 'W_1_inf_H1';
%  Norm_type = 'Final_t_L2_L2';
%  Norm_type = 'Final_t_L2_grad';
%  Norm_type = 'Final_t_L2_H1';
%  Norm_type = 'Final_t_H1_L2';   
%  Norm_type = 'Final_t_H1_grad';
%  Norm_type = 'Final_t_H1_H1';
 Norm_type = 'Energynorm';


DG_err=figure;


for k = 1:length(NO_time_step_vec )

err_P = NaN(length(T_Polydegree_vec),1); 
dof_P = NaN(length(T_Polydegree_vec),1);

for i=1 :length(T_Polydegree_vec)  
    
     load(['Error of ' num2str(NO_elem) ' rectangle space Elements combine ' num2str(NO_time_step_vec(k)) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree) ' and time basis ' num2str(T_Polydegree_vec(i)) ' basis.mat'])

     switch  Norm_type        
     case 'L2_L2'
        err_P(i,1)=L2_L2_err; 
     case 'L2_grad'
        err_P(i,1)=L2_grad_err;    
     case 'L2_H1'   
        err_P(i,1)=L2_H1_err;      
     case 'L_inf_L2'
        err_P(i,1)= L_inf_L2_err;    
     case 'L_inf_grad'
        err_P(i,1)= L_inf_grad_err;          
     case 'L_inf_H1'
        err_P(i,1)= L_inf_H1_err;   
        
     case 'H1_L2'
        err_P(i,1)=H1_L2_err;       
     case 'H1_grad'
        err_P(i,1)=H1_grad_err;           
     case 'H1_H1'   
        err_P(i,1)=H1_H1_err;      
     case 'W_1_inf_L2'
        err_P(i,1)= W_1_inf_L2_err;       
     case 'W_1_inf_grad'
        err_P(i,1)= W_1_inf_grad_err;
     case 'W_1_inf_H1'
        err_P(i,1)= W_1_inf_H1_err;    
       
     case 'Final_t_L2_L2'
        err_P(i,1)=Final_t_L2_L2_err;        
     case 'Final_t_L2_grad'
        err_P(i,1)=Final_t_L2_grad_err; 
     case 'Final_t_L2_H1'   
        err_P(i,1)=Final_t_L2_H1_err;      
     case 'Final_t_H1_L2'
        err_P(i,1)= Final_t_H1_L2_err;  
     case 'Final_t_H1_grad'
        err_P(i,1)= Final_t_H1_grad_err;        
     case 'Final_t_H1_H1'
        err_P(i,1)= Final_t_H1_H1_err;    
     case 'Energynorm'
        err_P(i,1)= Energynorm_err;          
     end
     dof_P(i) = T_dim;  
end

switch  Basis_Type
     case 'P'
semilogy(dof_P.^a,err_P(:,1),['-.s' num2str(color{k})],'LineWidth',1.5,'MarkerSize',10);
     case 'PQ'
semilogy(dof_P.^a,err_P(:,1),['-*' num2str(color{k})],'LineWidth',1.5,'MarkerSize',10);
     case 'Q'
semilogy(dof_P.^a,err_P(:,1),['-o' num2str(color{k})],'LineWidth',1.5,'MarkerSize',10);
end

%grid on;
hold on;
end


legend([ num2str(NO_time_step_vec(1)) ' FEM rect ' Basis_Type   ],...
       [ num2str(NO_time_step_vec(2)) ' FEM rect ' Basis_Type   ],...
       'Location','SouthWest')      
xlabel(' Dof ','FontSize',18);
        
switch  Norm_type
     case 'L2_L2'
        ylabel('|| u-u_{h}||_{L^2((0,T),L^2(\Omega))}','FontSize',18);
        title('L^2((0,T),L^2(\Omega)) norm error under h-refinement','FontSize',18)
     case 'L2_grad'
        ylabel('|| \nabla (u-u_{h})||_{L^2((0,T),L^2(\Omega))}','FontSize',18);
        title('L^2((0,T),grad(\Omega)) norm error under h-refinement','FontSize',18)        
     case 'L2_H1'    
        ylabel('|| u-u_{h}||_{L^2((0,T),H^1(\Omega))}','FontSize',18);        
        title('L^2((0,T),H^1(\Omega)) norm error under h-refinement','FontSize',18)
     case 'L_inf_L2'
        ylabel('|| u-u_{h}||_{L^\infty((0,T),L^2(\Omega))}','FontSize',18);       
        title('L^\infty((0,T),L^2(\Omega)) norm error under h-refinement','FontSize',18) 
     case 'L_inf_grad'
        ylabel('|| \nabla (u-u_{h})||_{L^\infty((0,T),L^2(\Omega))}','FontSize',18);       
        title('L^\infty((0,T),grad(\Omega)) norm error under h-refinement','FontSize',18)        
     case 'L_inf_H1'
        ylabel('|| u-u_{h}||_{L^\infty((0,T),H^1(\Omega))}','FontSize',18);        
        title('L^\infty((0,T),H^1(\Omega)) norm error under h-refinement','FontSize',18) 
     
    case 'H1_L2'
        ylabel('|| u-u_{h}||_{H^1((0,T),L^2(\Omega))}','FontSize',18);
        title('H^1((0,T),L^2(\Omega)) norm error under h-refinement','FontSize',18)
    case 'H1_grad'
        ylabel('|| \nabla (u-u_{h})||_{H^1((0,T),L^2(\Omega))}','FontSize',18);
        title('H^1((0,T),grad(\Omega)) norm error under h-refinement','FontSize',18)        
     case 'H1_H1'    
        ylabel('|| u-u_{h}||_{H^1((0,T),H^1(\Omega))}','FontSize',18);        
        title('H^1((0,T),H^1(\Omega)) norm error under h-refinement','FontSize',18)
     case 'H1_inf_L2'
        ylabel('|| u-u_{h}||_{H^1,\infty((0,T),L^2(\Omega))}','FontSize',18);       
        title('H^1,\infty((0,T),L^2(\Omega)) norm error under h-refinement','FontSize',18)   
     case 'H1_inf_grad'
        ylabel('|| \nabla (u-u_{h})||_{H^1,\infty((0,T),L^2(\Omega))}','FontSize',18);       
        title('H^1,\infty((0,T),grad(\Omega)) norm error under h-refinement','FontSize',18)    
     case 'H1_inf_H1'
        ylabel('|| u-u_{h}||_{H^1,\infty((0,T),H^1(\Omega))}','FontSize',18);        
        title('H^1,\infty((0,T),H^1(\Omega)) norm error under h-refinement','FontSize',18)   
     
     case 'Final_t_L2_L2'
        ylabel('|| (u-u_{h})(t_n^-)||_{L^2(\Omega))}','FontSize',18);
        title('L^2(\Omega) norm error at the final time under h-refinement','FontSize',18)
     case 'Final_t_L2_grad'
        ylabel('|| \nabla (u-u_{h})(t_n^-)||_{L^2(\Omega))}','FontSize',18);
        title('grad(\Omega) norm error at the final time under h-refinement','FontSize',18)        
     case 'Final t_L2_H1'    
        ylabel('|| (u-u_{h})(t_n^-)||_{H^1(\Omega))}','FontSize',18);        
        title('H^1(\Omega) norm error at the final time under h-refinement','FontSize',18)     
     case 'Final_t_H1_L2'
        ylabel('|| (u-u_{h})`(t_n^-)||_{L^2(\Omega))}','FontSize',18);
        title('L^2(\Omega) norm error at the final time under h-refinement','FontSize',18)
     case 'Final_t_H1_grad'
        ylabel('|| \nabla (u-u_{h})`(t_n^-)||_{L^2(\Omega))}','FontSize',18);
        title('grad(\Omega) norm error at the final time under h-refinement','FontSize',18)        
     case 'Final t_H1_H1'    
        ylabel('|| (u-u_{h})`(t_n^-)||_{H^1(\Omega))}','FontSize',18);        
        title('H^1(\Omega) norm error at the final time under h-refinement','FontSize',18)     
     case 'Energynorm'    
        ylabel('|| u-u_{h} ||_{Energy}','FontSize',18);        
        title('Energy norm error at the final time under h-refinement','FontSize',18)         
end


set(gca,'FontSize',15)
