clear all; close all;

%% different elements we use different order of basis
Mesh_type = 'rect';

Basis_Type  = 'Q';  % Q

eta=[1,0.001,0.0001,0.00001,0.000001,0];
NO_elem_vec = [4,16,64]; 
S_Polydegree = 2;


NO_time_step = 4;
T_Polydegree_vec = [2];
a=1/2; %%%%% Dof^{a}
color = {'b','r','g','k','c','m'};

% Norm_type = 'L2_L2'; 
% Norm_type = 'L2_grad';
% Norm_type = 'L2_H1';
% Norm_type = 'L_inf_L2';
%             Norm_type = 'L_inf_grad';
% Norm_type = 'L_inf_H1';
% Norm_type = 'H1_L2';   
% Norm_type = 'H1_grad'; 
             Norm_type = 'H1_H1';
%          Norm_type = 'W_1_inf_L2';
% Norm_type = 'W_1_inf_grad';
% Norm_type = 'W_1_inf_H1';
%       Norm_type = 'Final_t_H1_L2'; 
% Norm_type = 'Final_t_H1_H1';
% Norm_type = 'Final_t_L2_L2';
%     Norm_type = 'Final_t_L2_H1';%   
% Norm_type = 'Time_jump_u_err';
% Norm_type = 'Time_jump_gradu_err'; 
% Norm_type = 'Time_jump_dotu_err';
%       Norm_type = 'Energynorm';

DG_err=figure;

slope_Poly = NaN(1,length(T_Polydegree_vec)); 

for et = 1:length(eta)
for k = 1:length(T_Polydegree_vec)

err_P = NaN(length(NO_elem_vec),1); 
dof_P = NaN(length(NO_elem_vec),1);

for i=1 :length(NO_elem_vec)    
     load([num2str(eta(et)) ' Error of ' num2str(NO_elem_vec(i)) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
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

slope_Poly(et) = mean(slope_P1(end));
 
loglog(dof_P.^a,err_P(:,1),['-o' num2str(color{et})],'LineWidth',1.5,'MarkerSize',10);

T1=table(dof_P,err_P);
disp([ num2str(Norm_type) 'Error for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short e;
disp(T1);
h1= dof_P(2:length(NO_elem_vec));
disp([ num2str(Norm_type) 'Rate for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short;
T2=table(h1,slope_P1);
disp(T2);   

grid on;
hold on;
end
end

format shortE
legend(['\eta='  num2str(eta(1),'%0.e') ' slope ' num2str(-slope_Poly(1))],...
       ['\eta='  num2str(eta(2),'%0.e') '  slope ' num2str(-slope_Poly(2))],...
       ['\eta='  num2str(eta(3),'%0.e') '  slope ' num2str(-slope_Poly(3))],...
    ['\eta='  num2str(eta(4),'%0.e') '  slope ' num2str(-slope_Poly(4))],...
      ['\eta='  num2str(eta(5),'%0.e') '  slope ' num2str(-slope_Poly(5))],...
      ['\eta='  num2str(eta(6),'%0.e') ' slope ' num2str(-slope_Poly(6))],...
       'Location','SouthWest')      
xlabel(' DoFs ','FontSize',30);
        
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
        title('$L^\infty(I;H^1(\Omega))$ norm error for problem (III)','FontSize',18,'interpreter','LaTex')        
     case 'L_inf_H1'
        ylabel('|| u-u_{h}||_{L^\infty((0,T),H^1(\Omega))}','FontSize',18);        
        title('L^\infty((0,T),H^1(\Omega)) norm error under h-refinement','FontSize',18) 
     
    case 'H1_L2'
        ylabel('|| u-u_{h}||_{H^1((0,T),L^2(\Omega))}','FontSize',18);
        title('H^1((0,T),L^2(\Omega)) norm error under h-refinement','FontSize',18)
    case 'H1_grad'
        ylabel('|| \nabla (u-u_{h})||_{H^1((0,T),L^2(\Omega))}','FontSize',18);
        title('H^1((0,T),grad(\Omega)) norm error under h-refinement with $p=2$','FontSize',18)        
     case 'H1_H1'    
        ylabel('$||e||_{\mathrm{E}}$','FontSize',18,'interpreter','LaTex');        
        title('$||e\|_{\mathrm{E}}$ norm error for spatial $h$-refinement with $p=2$','FontSize',18,'interpreter','LaTex')
     case 'W_1_inf_L2'
        ylabel('$||e||_{W^{1,\infty}(I;L^2(\Omega))}$','FontSize',18,'interpreter','LaTex');       
        title('$W^{1,\infty}(I;L^2(\Omega))$ norm error for problem (III)','FontSize',18,'interpreter','LaTex')   
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
        ylabel('$|||e|||$','FontSize',18,'interpreter','LaTex');        
        title('DG norm error under for problem (III)','FontSize',18,'interpreter','LaTex')          
end


set(gca,'FontSize',10.5)
