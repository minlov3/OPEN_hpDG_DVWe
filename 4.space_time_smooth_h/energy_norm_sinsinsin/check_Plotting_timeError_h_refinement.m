clear all; close all;

%% different elements we use different order of basis
Mesh_type = 'rect';
NO_elem = 64; 
S_Polydegree = 8;
Basis_Type  = 'Q';  % Q


NO_time_step_vec = [2,4,8,16,32];
 
T_Polydegree_vec = [1,2,3,4];
a=1; %%%%% Dof^{a}
color = {'b','r','g','k'};

% Norm_type = 'L2_L2'; 
% Norm_type = 'L2_grad';
% Norm_type = 'L2_H1';
% Norm_type = 'L_inf_L2';
%Norm_type = 'L_inf_grad';
% Norm_type = 'L_inf_H1';
% Norm_type = 'H1_L2';   
% Norm_type = 'H1_grad'; 
Norm_type = 'H1_H1';
%Norm_type = 'W_1_inf_L2';
% Norm_type = 'W_1_inf_grad';
% Norm_type = 'W_1_inf_H1';
 Norm_type = 'Final_t_H1_L2'; 
% Norm_type = 'Final_t_H1_H1';
% Norm_type = 'Final_t_L2_L2';
Norm_type = 'Final_t_L2_H1';%   
% Norm_type = 'Time_jump_u_err';
% Norm_type = 'Time_jump_gradu_err'; 
% Norm_type = 'Time_jump_dotu_err';
%Norm_type = 'Energynorm';

DG_err=figure;

slope_Poly = NaN(1,length(T_Polydegree_vec)); 
for k = 1:length(T_Polydegree_vec)

err_P = NaN(length(NO_time_step_vec),1); 
dof_P = NaN(length(NO_time_step_vec),1);

for i=1 :length(NO_time_step_vec)    
     load(['Error of ' num2str(NO_elem) ' rectangle space Elements combine ' num2str(NO_time_step_vec(i)) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree) ' and time basis ' num2str(T_Polydegree_vec(k)) ' basis.mat'])

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
     case 'Time_jump_u_err'
        err_P(i,1)= Time_jump_u_err;  
     case 'Time_jump_gradu_err'
        err_P(i,1)= Time_jump_gradu_err;  
     case 'Time_jump_dotu_err'
        err_P(i,1)= Time_jump_dotu_err;  
     case 'Energynorm'
        err_P(i,1)= Energynorm_err;   
     end
dof_P(i) = dim_FEM;  
end

logerr_P1 = abs(log(err_P(:,1))); 
slope_P1 = abs((logerr_P1(2:end)-logerr_P1(1:end-1))./(log(dof_P(2:end).^a)-log(dof_P(1:end-1).^a)));
slope_P1 = roundn(slope_P1,-2);

slope_Poly(k) = mean(slope_P1(end));
if k==3
slope_Poly(k) = mean(slope_P1(end)); 
elseif k==4
slope_Poly(k) = mean(slope_P1(end));
end
loglog(dof_P.^a,err_P(:,1),['-o' num2str(color{k})],'LineWidth',1.5,'MarkerSize',10);

T1=table(dof_P,err_P);
disp([ num2str(Norm_type) 'Error for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short e;
disp(T1);
h1= dof_P(2:length(NO_time_step_vec));
disp([ num2str(Norm_type) 'Rate for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short;
T2=table(h1,slope_P1);
disp(T2);   

grid on;
hold on;
end


legend(['P'  num2str(T_Polydegree_vec(1)) ' slope ' num2str(-slope_Poly(1))],...
       ['P'  num2str(T_Polydegree_vec(2)) ' slope ' num2str(-slope_Poly(2))],...
       ['P'  num2str(T_Polydegree_vec(3)) ' slope ' num2str(-slope_Poly(3))],...
       ['P'  num2str(T_Polydegree_vec(4)) ' slope ' num2str(-slope_Poly(4))],...
       'Location','SouthWest')      
xlabel(' Dof ','FontSize',30);
        
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
        ylabel('$||\nabla e||_{L^\infty(I;L^2(\Omega))}$','FontSize',18,'interpreter','LaTex');       
        title('$L^\infty(I;H^1(\Omega))$ norm error for smooth solution under $h$-refinement','FontSize',18,'interpreter','LaTex')        
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
        ylabel('$||e||_{H^1(I;H^1(\Omega))}$','FontSize',18,'interpreter','LaTex');        
        title('$H^1(I;H^1(\Omega))$ norm error for smooth solution under $h$-refinement','FontSize',18,'interpreter','LaTex')
     case 'W_1_inf_L2'
        ylabel('$||e||_{W^{1,\infty}(I;L^2(\Omega))}$','FontSize',18,'interpreter','LaTex');       
        title('$W^{1,\infty}(I;L^2(\Omega))$ norm error for smooth solution under $h$-refinement','FontSize',18,'interpreter','LaTex')   
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
        ylabel('$||\nabla e(t_n^-)||_{L^2(\Omega))}$','FontSize',18,'interpreter','LaTex');
        title('$||\nabla e(t_n^-)||_{L^2(\Omega))}$ for smooth solution under $h$-refinement','FontSize',18,'interpreter','LaTex')        
     case 'Final_t_L2_H1'    
        ylabel('$||\nabla e(t_N^-)||_{L^2(\Omega))}$','FontSize',18,'interpreter','LaTex');
        title('$||\nabla e(t_N^-)||_{L^2(\Omega))}$ for smooth solution under $h$-refinement','FontSize',18,'interpreter','LaTex')      
     case 'Final_t_H1_L2'
        ylabel('$||\dot{e}(t_N^-)||_{L^2(\Omega))}$','FontSize',18,'interpreter','LaTex');
        title('$||\dot{e}(t_N^-)||_{L^2(\Omega))}$ for smooth solution under $h$-refinement','FontSize',18,'interpreter','LaTex');
  
     case 'Final_t_H1_grad'
        ylabel('|| \nabla (u-u_{h})`(t_n^-)||_{L^2(\Omega))}','FontSize',18);
        title('grad(\Omega) norm error at the final time under h-refinement','FontSize',18)        
     case 'Final t_H1_H1'    
        ylabel('|| (u-u_{h})`(t_n^-)||_{H^1(\Omega))}','FontSize',18);        
        title('H^1(\Omega) norm error at the final time under h-refinement','FontSize',18)    
     case 'Energynorm'    
        ylabel('$||e||_{\mathrm{DG}}$','FontSize',18,'interpreter','LaTex');        
        title('DG norm error under for smooth solution $h$-refinement','FontSize',18,'interpreter','LaTex')          
end


set(gca,'FontSize',10.5)
