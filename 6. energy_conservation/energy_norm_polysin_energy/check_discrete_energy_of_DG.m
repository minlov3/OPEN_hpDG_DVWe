
%% 
clear all;

time_starting = 0;   
time_end = 100;  
color = {'b','r','g','k','c','m'};
shape = {'-','-.','-o','-d'};
Basis_Type  = 'Q';  % Q, Se
Mesh_type = 'rect';   % rect, 
NO_elem = 4; 
S_Polydegree = 2; 

NO_time_step_vector = [200]; 
T_Polydegree_vector = [2,3];


for t=1:length(T_Polydegree_vector)
    

for s=1:length(NO_time_step_vector)

NO_time_step  = NO_time_step_vector(s);

ene =NaN(  NO_time_step ,1); 
time = NaN( NO_time_step ,1);

load(['Energy of ' num2str(NO_elem ) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree) ' and time basis ' num2str(T_Polydegree_vector(t))  ' basis.mat'])

for i=1 : NO_time_step 
ene(i)=DG_ene_vector(i) ;         
time(i) = i.*( (time_end-time_starting)./NO_time_step+time_starting ) ;
end

logerr_P1 = abs(log(ene(:,1))); 
slope_P1 = 0.5.*((logerr_P1(2:end)-logerr_P1(1:end-1))./(time(2:end)-time(1:end-1)));
slope_P1 = roundn(slope_P1,-3);

 

plot(time,ene,[ num2str(color{T_Polydegree_vector(t)}), num2str(shape{1}) ],'LineWidth',1,'MarkerSize',8);
hold on;
end





end

format shortE
legend([ '$P_{2}$'  ],...
       [ '$P_{3}$' ],...
       'Location','NorthEast','interpreter','LaTex');
%xlabel('Polynomial order of basis','FontSize',18);
xlabel('time $t_n$','FontSize',18,'interpreter','LaTex');
ylabel('$||{u}^{\prime}_{\tau}(t_n^-)||^2+||\nabla {u}_{\tau}(t_n^-)||^2$ ','FontSize',18,'interpreter','LaTex');
title(['$N_T = $' num2str(NO_time_step_vector(1))],'FontSize',18,'interpreter','LaTex');
set(gca,'FontSize',13.5); 




 
%matlabpool close;

