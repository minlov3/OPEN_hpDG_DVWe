clear all; close all;

%% different elements we use different order of basis
Mesh_type = 'rect';
NO_elem = 4; 
color = {'b','r','g','k','m','c','y'};
shape = {'-*','-o','-x'};
S_Polydegree = 2;

Basis_Type  = 'Q';  % Q
a=1/2; %%%%% Dof^{a}

Time_Step_Vector = [1,2,3,4,5,6,7,8,9];

gf_vec = [0.15,0.17,0.2,0.25,0.3,0.35,0.5];     % temporal mesh grading factor of 
mu = 2;  % polynomial degree increasing factor


%Norm_type = 'L2_L2'; 
% Norm_type = 'L2_grad';
% Norm_type = 'L2_H1';
% Norm_type = 'L_inf_L2';
      Norm_type = 'L_inf_grad';
% Norm_type = 'L_inf_H1';
% Norm_type = 'H1_L2';   
% Norm_type = 'H1_grad'; 
%        Norm_type = 'H1_H1';
 %         Norm_type = 'W_1_inf_L2';
% Norm_type = 'W_1_inf_grad';
% Norm_type = 'W_1_inf_H1';
% Norm_type = 'Final_t_H1_L2'; 
% Norm_type = 'Final_t_H1_grad'; 
% Norm_type = 'Final_t_H1_H1';
% Norm_type = 'Final_t_L2_L2';
% Norm_type = 'Final_t_L2_grad';
% Norm_type = 'Final_t_L2_H1';
 %         Norm_type = 'Energynorm';


DG_err=figure;



for OA = 1: length(gf_vec) 

gf = gf_vec(OA);
err_P = NaN(length(Time_Step_Vector),1); 
dof_P = NaN(length(Time_Step_Vector),1);

for i = 1:length(Time_Step_Vector)    
    
      load( [ 'Error of ' num2str(NO_elem) ' rectangle space Elements with ' Basis_Type num2str(S_Polydegree) ' space basis, combining '...
             num2str(Time_Step_Vector(i)) ' time step with sigma = ' num2str(gf) ' and  mu = ' num2str(mu) '.mat'])

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
        err_P(i,1)= Final_t_H1_H1_err;          
     end
dof_P(i) = dim_FEM;  
end

switch  Basis_Type
     case 'P'
loglog(dof_P.^a,err_P(:,1),['-.s' num2str(color{OA})],'LineWidth',1.5,'MarkerSize',10);
     case 'PQ'
loglog(dof_P.^a,err_P(:,1),['-*' num2str(color{OA})],'LineWidth',1.5,'MarkerSize',10);
     case 'Q'
semilogy(dof_P.^a,err_P(:,1),['-*' num2str(color{OA})],'LineWidth',1,'MarkerSize',8);
end

%grid on;
hold on;

end
 
legend([ '\sigma= ' num2str(gf_vec(1)) ],...
       [ '\sigma= ' num2str(gf_vec(2)) ],...
       [ '\sigma= ' num2str(gf_vec(3)) ],...
       [ '\sigma= ' num2str(gf_vec(4)) ],... 
       [ '\sigma= ' num2str(gf_vec(5)) ],... 
       [ '\sigma= ' num2str(gf_vec(6)) ],...
       [ '\sigma= ' num2str(gf_vec(7)) ],...
'Location','SouthWest')      
xlabel(' DoFs^{1/2} ','FontSize',18);
        
switch  Norm_type
     case 'L2_L2'
        ylabel('|| u-u_{h}||_{L^2((0,T),L^2(\Omega))}','FontSize',18);
        title('L^2((0,T),L^2(\Omega)) norm error under hp-refinement','FontSize',18)
     case 'L2_grad'
        ylabel('|| \nabla (u-u_{h})||_{L^2((0,T),L^2(\Omega))}','FontSize',18);
        title('L^2((0,T),grad(\Omega)) norm error under hp-refinement','FontSize',18)        
     case 'L2_H1'    
        ylabel('|| u-u_{h}||_{L^2((0,T),H^1(\Omega))}','FontSize',18);        
        title('L^2((0,T),H^1(\Omega)) norm error under hp-refinement','FontSize',18)
     case 'L_inf_L2'
        ylabel('|| u-u_{h}||_{L^\infty((0,T),L^2(\Omega))}','FontSize',18);       
        title('L^\infty((0,T),L^2(\Omega)) norm error under hp-refinement','FontSize',18) 
     case 'L_inf_grad'
        ylabel('$|| \nabla e||_{L^{\infty}(I;L^2(\Omega))}$','FontSize',18,'interpreter','LaTex');       
        title('$L^{\infty}(I;H^1(\Omega))$ norm error under $\tau q$-refinement','FontSize',18,'interpreter','LaTex');        
     case 'L_inf_H1'
        ylabel('|| u-u_{h}||_{L^\infty((0,T),L^2(\Omega))}','FontSize',18);        
        title('L^\infty((0,T),H^1(\Omega)) norm error under hp-refinement','FontSize',18) 
     
     case 'H1_L2'
        ylabel('|| u-u_{h}||_{H^1((0,T),L^2(\Omega))}','FontSize',18);
        title('H^1((0,T),L^2(\Omega)) norm error under hp-refinement','FontSize',18)
     case 'H1_grad'
        ylabel('|| \nabla (u-u_{h})||_{H^1((0,T),L^2(\Omega))}','FontSize',18);
        title('H^1((0,T),grad(\Omega)) norm error under hp-refinement','FontSize',18)        
     case 'H1_H1'    
        ylabel('$\|e\|_{\mathrm{E}}$','FontSize',18,'interpreter','LaTex');         
        title('$\|e\|_{\mathrm{E}}$ norm error under $\tau q$-refinement','FontSize',18,'interpreter','LaTex'); 
     case 'W_1_inf_L2'
        ylabel('$||e||_{W^{1,\infty}(I;L^2(\Omega))}$','FontSize',18,'interpreter','LaTex');       
        title('$W^{1,\infty}(I;L^2(\Omega))$ norm error under $\tau q$-refinement','FontSize',18,'interpreter','LaTex');   
     case 'W1_inf_grad'
        ylabel('|| \nabla (u-u_{h})||_{H^1,\infty((0,T),L^2(\Omega))}','FontSize',18);       
        title('H^1,\infty((0,T),grad(\Omega)) norm error under hp-refinement','FontSize',18)    
     case 'H1_inf_H1'
        ylabel('|| u-u_{h}||_{H^1,\infty((0,T),H^1(\Omega))}','FontSize',18);        
        title('H^1,\infty((0,T),H^1(\Omega)) norm error under hp-refinement','FontSize',18)   
     
     case 'Final_t_L2_L2'
        ylabel('|| (u-u_{h})(t_n^-)||_{L^2(\Omega))}','FontSize',18);
        title('L^2(\Omega) norm error at the final time under hp-refinement','FontSize',18)
     case 'Final_t_L2_grad'
        ylabel('|| \nabla (u-u_{h})(t_n^-)||_{L^2(\Omega))}','FontSize',18);
        title('grad(\Omega) norm error at the final time under hp-refinement','FontSize',18)        
     case 'Final t_L2_H1'    
        ylabel('|| (u-u_{h})(t_n^-)||_{H^1(\Omega))}','FontSize',18);        
        title('H^1(\Omega) norm error at the final time under hp-refinement','FontSize',18)     
     case 'Final_t_H1_L2'
        ylabel('|| (u-u_{h})`(t_n^-)||_{L^2(\Omega))}','FontSize',18);
        title('L^2(\Omega) norm error at the final time under hp-refinement','FontSize',18)
     case 'Final_t_H1_grad'
        ylabel('|| \nabla (u-u_{h})`(t_n^-)||_{L^2(\Omega))}','FontSize',18);
        title('grad(\Omega) norm error at the final time under hp-refinement','FontSize',18)        
     case 'Final t_H1_H1'    
        ylabel('|| (u-u_{h})`(t_n^-)||_{H^1(\Omega))}','FontSize',18);        
        title('H^1(\Omega) norm error at the final time under hp-refinement','FontSize',18)           
     case 'Energynorm'    
        ylabel('$|||e|||$','FontSize',18,'interpreter','LaTex');        
        title('DG norm error under $\tau q$-refinement','FontSize',18,'interpreter','LaTex'); 
     end

set(gca,'FontSize',10.5)