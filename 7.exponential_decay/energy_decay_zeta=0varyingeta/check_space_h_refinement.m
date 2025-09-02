
%% 
clear all;

S_Polydegree = 2; 
T_Polydegree = 2;
Mesh_type = 'rect';   % rect, 

NO_elem_vec = [4,16,64]; 
NO_time_step  = 20;

Basis_Type  = 'Q';  % Q, Se

Norm_type = 'L2_L2';   
%Norm_type = 'L2_H1';
%Norm_type = 'L_inf_L2';
%Norm_type = 'L_inf_H1';
%Norm_type = 'H1_L2';   
%Norm_type = 'H1_H1';
%Norm_type = 'W_1_inf_L2';
%Norm_type = 'W_1_inf_H1';
%Norm_type = 'Final_t_H1_L2';   
%Norm_type = 'Final_t_H1_H1';
%Norm_type = 'Final_t_L2_L2';
%Norm_type = 'Final_t_L2_H1';

err =NaN(length(NO_elem_vec),1); 
dof = NaN(length(NO_elem_vec),1);


switch  Norm_type
    
case 'L2_L2'
for i=1 :length(NO_elem_vec)
load(['Error of ' num2str(NO_elem_vec(i)) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree) ' and time basis ' num2str(T_Polydegree) ' basis.mat'])
err(i)=L2_L2_err ;         
dof(i) = S_dim;
end
eOrder=NaN(length(NO_elem_vec)-1,1);
for i=1:length(NO_elem_vec)-1
    eOrder(i)=log(err(i)/err(i+1))/log(2);
end
T1=table(NO_elem_vec',err);
disp([ num2str(Norm_type) 'Error for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short e;
disp(T1);
h1= NO_elem_vec(2:length(NO_elem_vec));
disp([ num2str(Norm_type) 'Rate for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short;
T2=table(h1',eOrder);
disp(T2);      
switch Basis_Type 
case 'Q'
loglog(dof.^(1/2),err,'-mo','LineWidth',1,'MarkerSize',8);
case 'Se'
loglog(dof.^(1/2),err,'-go','LineWidth',1,'MarkerSize',8);
end
legend([ num2str(NO_time_step) ' time steps with ' num2str(S_Polydegree) ' time basis '])
%xlabel('Polynomial order of basis','FontSize',18);
xlabel('log N_s^{1/2}','FontSize',18);
ylabel('|| u-u_{h}||_{L^2((0,T),L^2(\Omega))}','FontSize',18);
title('space h refinement','FontSize',18); set(gca,'FontSize',15); 





case 'L2_H1'
for i=1 :length(NO_elem_vec)
load(['Error of ' num2str(NO_elem_vec(i)) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree) ' and time basis ' num2str(T_Polydegree) ' basis.mat'])
err(i)=L2_H1_err ;         
dof(i) = dim_FEM;
end
eOrder=NaN(length(NO_elem_vec)-1,1);
for i=1:length(NO_elem_vec)-1
    eOrder(i)=log(err(i)/err(i+1))/log(2);
end
T1=table(NO_elem_vec',err);
disp([ num2str(Norm_type) 'Error for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short e;
disp(T1);
h1= NO_elem_vec(2:length(NO_elem_vec));
disp([ num2str(Norm_type) 'Rate for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short;
T2=table(h1',eOrder);
disp(T2);  

case 'L_inf_L2'
for i=1 :length(NO_elem_vec)
load(['Error of ' num2str(NO_elem_vec(i)) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree) ' and time basis ' num2str(T_Polydegree) ' basis.mat'])
err(i)=L_inf_L2_err ;         
dof(i) = dim_FEM;
end
eOrder=NaN(length(NO_elem_vec)-1,1);
for i=1:length(NO_elem_vec)-1
    eOrder(i)=log(err(i)/err(i+1))/log(2);
end
T1=table(NO_elem_vec',err);
disp([ num2str(Norm_type) 'Error for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short e;
disp(T1);
h1= NO_elem_vec(2:length(NO_elem_vec));
disp([ num2str(Norm_type) 'Rate for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short;
T2=table(h1',eOrder);
disp(T2);  

case 'L_inf_H1'
for i=1 :length(NO_elem_vec)
load(['Error of ' num2str(NO_elem_vec(i)) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree) ' and time basis ' num2str(T_Polydegree) ' basis.mat'])
err(i)=L_inf_H1_err ;         
dof(i) = dim_FEM;
end
eOrder=NaN(length(NO_elem_vec)-1,1);
for i=1:length(NO_elem_vec)-1
    eOrder(i)=log(err(i)/err(i+1))/log(2);
end
T1=table(NO_elem_vec',err);
disp([ num2str(Norm_type) 'Error for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short e;
disp(T1);
h1= NO_elem_vec(2:length(NO_elem_vec));
disp([ num2str(Norm_type) 'Rate for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short;
T2=table(h1',eOrder);
disp(T2);  

case 'H1_L2'
for i=1 :length(NO_elem_vec)
load(['Error of ' num2str(NO_elem_vec(i)) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree) ' and time basis ' num2str(T_Polydegree) ' basis.mat'])
err(i)=H1_L2_err ;         
dof(i) = dim_FEM;
end
eOrder=NaN(length(NO_elem_vec)-1,1);
for i=1:length(NO_elem_vec)-1
    eOrder(i)=log(err(i)/err(i+1))/log(2);
end
T1=table(NO_elem_vec',err);
disp([ num2str(Norm_type) 'Error for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short e;
disp(T1);
h1= NO_elem_vec(2:length(NO_elem_vec));
disp([ num2str(Norm_type) 'Rate for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short;
T2=table(h1',eOrder);
disp(T2);  

case 'H1_H1'
for i=1 :length(NO_elem_vec)
load(['Error of ' num2str(NO_elem_vec(i)) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree) ' and time basis ' num2str(T_Polydegree) ' basis.mat'])
err(i)=H1_H1_err ;         
dof(i) = dim_FEM;
end
eOrder=NaN(length(NO_elem_vec)-1,1);
for i=1:length(NO_elem_vec)-1
    eOrder(i)=log(err(i)/err(i+1))/log(2);
end
T1=table(NO_elem_vec',err);
disp([ num2str(Norm_type) 'Error for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short e;
disp(T1);
h1= NO_elem_vec(2:length(NO_elem_vec));
disp([ num2str(Norm_type) 'Rate for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short;
T2=table(h1',eOrder);
disp(T2);  

case 'W_1_inf_L2'
for i=1 :length(NO_elem_vec)
load(['Error of ' num2str(NO_elem_vec(i)) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree) ' and time basis ' num2str(T_Polydegree) ' basis.mat'])
err(i)=W_1_inf_L2_err ;         
dof(i) = dim_FEM;
end
eOrder=NaN(length(NO_elem_vec)-1,1);
for i=1:length(NO_elem_vec)-1
    eOrder(i)=log(err(i)/err(i+1))/log(2);
end
T1=table(NO_elem_vec',err);
disp([ num2str(Norm_type) 'Error for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short e;
disp(T1);
h1= NO_elem_vec(2:length(NO_elem_vec));
disp([ num2str(Norm_type) 'Rate for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short;
T2=table(h1',eOrder);
disp(T2);  

case 'W_1_inf_H1'
for i=1 :length(NO_elem_vec)
load(['Error of ' num2str(NO_elem_vec(i)) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree) ' and time basis ' num2str(T_Polydegree) ' basis.mat'])
err(i)=W_1_inf_H1_err ;         
dof(i) = dim_FEM;
end
eOrder=NaN(length(NO_elem_vec)-1,1);
for i=1:length(NO_elem_vec)-1
    eOrder(i)=log(err(i)/err(i+1))/log(2);
end
T1=table(NO_elem_vec',err);
disp([ num2str(Norm_type) 'Error for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short e;
disp(T1);
h1= NO_elem_vec(2:length(NO_elem_vec));
disp([ num2str(Norm_type) 'Rate for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short;
T2=table(h1',eOrder);
disp(T2);  

case 'Final_t_L2_L2'
for i=1 :length(NO_elem_vec)
load(['Error of ' num2str(NO_elem_vec(i)) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree) ' and time basis ' num2str(T_Polydegree) ' basis.mat'])
err(i)=Final_t_L2_L2_err ;         
dof(i) = dim_FEM;
end
eOrder=NaN(length(NO_elem_vec)-1,1);
for i=1:length(NO_elem_vec)-1
    eOrder(i)=log(err(i)/err(i+1))/log(2);
end
T1=table(NO_elem_vec',err);
disp([ num2str(Norm_type) 'Error for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short e;
disp(T1);
h1= NO_elem_vec(2:length(NO_elem_vec));
disp([ num2str(Norm_type) 'Rate for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short;
T2=table(h1',eOrder);
disp(T2);  

case 'Final_t_L2_H1'
for i=1 :length(NO_elem_vec)
load(['Error of ' num2str(NO_elem_vec(i)) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree) ' and time basis ' num2str(T_Polydegree) ' basis.mat'])
err(i)=Final_t_L2_H1_err ;         
dof(i) = dim_FEM;
end
eOrder=NaN(length(NO_elem_vec)-1,1);
for i=1:length(NO_elem_vec)-1
    eOrder(i)=log(err(i)/err(i+1))/log(2);
end
T1=table(NO_elem_vec',err);
disp([ num2str(Norm_type) 'Error for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short e;
disp(T1);
h1= NO_elem_vec(2:length(NO_elem_vec));
disp([ num2str(Norm_type) 'Rate for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short;
T2=table(h1',eOrder);
disp(T2); 

case 'Final_t_H1_L2'
for i=1 :length(NO_elem_vec)
load(['Error of ' num2str(NO_elem_vec(i)) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree) ' and time basis ' num2str(T_Polydegree) ' basis.mat'])
err(i)=Final_t_H1_L2_err ;         
dof(i) = dim_FEM;
end
eOrder=NaN(length(NO_elem_vec)-1,1);
for i=1:length(NO_elem_vec)-1
    eOrder(i)=log(err(i)/err(i+1))/log(2);
end
T1=table(NO_elem_vec',err);
disp([ num2str(Norm_type) 'Error for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short e;
disp(T1);
h1= NO_elem_vec(2:length(NO_elem_vec));
disp([ num2str(Norm_type) 'Rate for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short;
T2=table(h1',eOrder);
disp(T2); 

case 'Final_t_H1_H1'
for i=1 :length(NO_elem_vec)
load(['Error of ' num2str(NO_elem_vec(i)) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree) ' and time basis ' num2str(T_Polydegree) ' basis.mat'])
err(i)=Final_t_H1_H1_err ;         
dof(i) = dim_FEM;
end
eOrder=NaN(length(NO_elem_vec)-1,1);
for i=1:length(NO_elem_vec)-1
    eOrder(i)=log(err(i)/err(i+1))/log(2);
end
T1=table(NO_elem_vec',err);
disp([ num2str(Norm_type) 'Error for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short e;
disp(T1);
h1= NO_elem_vec(2:length(NO_elem_vec));
disp([ num2str(Norm_type) 'Rate for space polydegree= ' num2str(S_Polydegree) ' ' Basis_Type ...
      ' rectangle Elements - ' ' and time polydegree=' num2str(T_Polydegree) ' with respect to space h refinement ' ]);
format short;
T2=table(h1',eOrder);
disp(T2); 


saveas(L2_L2,['viscouswave_L2_norm_error_h_refine_' ],'fig');

print(L2_L2,'-depsc',['viscouswave_L2_norm_error_h_refine_' '.eps']);

end
%matlabpool close;

