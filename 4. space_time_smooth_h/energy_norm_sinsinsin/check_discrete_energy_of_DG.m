
%% 
clear all;

time_starting = 0;   
time_end = 1;  


S_Polydegree = 2; 
T_Polydegree = 2;

Mesh_type = 'rect';   % rect, 
NO_elem = 16; 
NO_time_step  = 20;
Basis_Type  = 'Q';  % Q, Se

ene =NaN(  NO_time_step ,1); 
time = NaN( NO_time_step ,1);

load(['Energy of ' num2str(NO_elem ) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree) ' and time basis ' num2str(T_Polydegree) ' basis.mat'])

for i=1 : NO_time_step 
ene(i)=DG_ene_vector(i) ;         
time(i) = i.*( (time_end-time_starting)./NO_time_step+time_starting ) ;
end
    
 
semilogy(time,ene,'r','LineWidth',1,'MarkerSize',8);
 
legend([ num2str(NO_time_step) ' time steps with P' num2str(S_Polydegree) ' time basis '])
%xlabel('Polynomial order of basis','FontSize',18);
xlabel('time','FontSize',18);
ylabel('Energy','FontSize',18);
title('Energy','FontSize',18); set(gca,'FontSize',15); 

set(gca,'FontSize',15)
 
%matlabpool close;

