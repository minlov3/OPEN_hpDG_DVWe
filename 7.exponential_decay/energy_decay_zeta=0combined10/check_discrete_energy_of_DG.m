
%% 
clear all;

time_starting = 0;   
time_end = 1;  
color = {'b','r','g','k','c','m'};
shape = {'-*','-s','-o','-d'};
Basis_Type  = 'Q';  % Q, Se
Mesh_type = 'rect';   % rect, 
NO_elem = 4; 
S_Polydegree = 2; 
NO_time_step_vector = [50]; 
T_Polydegree_vector = [6];
gamma=[0,1,5,10];
for ze = 1:length(gamma)
for t=1:length(T_Polydegree_vector)
    
for s=1:length(NO_time_step_vector)

NO_time_step  = NO_time_step_vector(s);

ene =NaN(  NO_time_step ,1); 
time = NaN( NO_time_step ,1);

load([num2str(gamma(ze)) ' Energy of ' num2str(NO_elem ) ' rectangle space Elements combine ' num2str(NO_time_step) ' time elements with space basis '...
            Basis_Type num2str(S_Polydegree) ' and time basis ' num2str(T_Polydegree_vector(t)) ' basis.mat'])

for i=1 : NO_time_step 
ene(i)=DG_ene_vector(i) ;         
time(i) = i.*( (time_end-time_starting)./NO_time_step+time_starting ) ;
end

logerr_P1 = abs(log(ene(:,1))); 
slope_P1 = 0.5.*abs((logerr_P1(2:end)-logerr_P1(1:end-1))./(time(2:end)-time(1:end-1)));
slope_P1 = roundn(slope_P1,-3);

slope_Poly(ze) = mean(slope_P1(end));

semilogy(time,ene,[ num2str(color{ze}), num2str(shape{ze}) ],'LineWidth',1,'MarkerSize',8);
hold on;
end
end
end
format shortE
legend([ '\gamma=' num2str(gamma(1)) ',  \eta=0.5   with slope ' num2str(-slope_Poly(1))  ],...
       [ '\gamma=' num2str(gamma(2)) ',  \eta=0.45 with slope ' num2str(-slope_Poly(2))  ],...
       [ '\gamma=' num2str(gamma(3)) ',  \eta=0.25 with slope ' num2str(-slope_Poly(3))  ],...
       [ '\gamma=' num2str(gamma(4)) ', \eta=0     with slope ' num2str(-slope_Poly(4))  ],...
       'Location','SouthWest')
%xlabel('Polynomial order of basis','FontSize',18);
xlabel('time $t_n$','FontSize',18,'interpreter','LaTex');
ylabel('$||{u}^{\prime}_{\tau}(t_n^-)||$','FontSize',18,'interpreter','LaTex');
title(['$||{u}^{\prime}_{\tau}(t_n^-)||$  under the coefficient $\gamma+20\eta=10$, $q=$' num2str(T_Polydegree_vector(1)), '$, N_T=50$'],'FontSize',18,'interpreter','LaTex');
set(gca,'FontSize',13.5); 




 
%matlabpool close;

