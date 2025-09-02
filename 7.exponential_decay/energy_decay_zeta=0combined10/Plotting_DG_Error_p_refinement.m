%% plotting errors for the DG for p refinement

clear all;close all;

%parameters

penalty=5; 

NO_elem = [8,64,512,4096];

Polynomial_degree = [15,12,8,6];  color = {'b','r','g','k'};





L2=figure;


for i=1:size(NO_elem,2)

errL_2 =NaN(Polynomial_degree(i),1); dof = NaN(Polynomial_degree(i),1);


for j=1:Polynomial_degree(i)
    
    
load(['Error ' num2str(NO_elem(i)) ' cube Elements penalty ' num2str(penalty) ' P' num2str(j) ' basis.mat'])


errL_2(j)=L2_err;         dof(j) = dim_FEM;

    
    
end


semilogy(dof.^(1/3),errL_2,['-h' num2str(color{i})],'LineWidth',1,'MarkerSize',8);

hold on;

end

hold off;

legend([num2str(NO_elem(1)) ' P basis'],...
       [num2str(NO_elem(2)) ' P basis'],...
       [num2str(NO_elem(3)) ' P basis'],...
       [num2str(NO_elem(4)) ' P basis']);

xlabel('Dof^{1/3}','FontSize',18);

ylabel('||u-u_h ||_{L^2{(\Omega)}}','FontSize',18);

set(gca,'FontSize',15);

title(['L2 norm error under p-refinement penalty ' num2str(penalty) ],'FontSize',20)


saveas(L2,['Poisson_L2_norm_error_p_refine_penalty_' num2str(penalty) ],'fig');

print(L2,'-depsc',['Poisson_L2_norm_error_p_refine_penalty_' num2str(penalty) '.eps']);







H1=figure;


for i=1:size(NO_elem,2)

errL_2 =NaN(Polynomial_degree(i),1); dof = NaN(Polynomial_degree(i),1);


for j=1:Polynomial_degree(i)
    
    
load(['Error ' num2str(NO_elem(i)) ' cube Elements penalty ' num2str(penalty) ' P' num2str(j) ' basis.mat'])


errL_2(j)=H1_err;         dof(j) = dim_FEM;

    
    
end


semilogy(dof.^(1/3),errL_2,['-h' num2str(color{i})],'LineWidth',1,'MarkerSize',8);

hold on;

end

hold off;

legend([num2str(NO_elem(1)) ' P basis'],...
       [num2str(NO_elem(2)) ' P basis'],...
       [num2str(NO_elem(3)) ' P basis'],...
       [num2str(NO_elem(4)) ' P basis']);

xlabel('Dof^{1/3}','FontSize',18);

ylabel('|u-u_h |_{H^1{(\Omega)}}','FontSize',18);


set(gca,'FontSize',15);

title(['H1 semi-norm error under p-refinement penalty ' num2str(penalty) ],'FontSize',20)

saveas(H1,['Poisson_H1_semi_norm_error_p_refine_penalty_' num2str(penalty) ],'fig');

print(H1,'-depsc',['Poisson_H1_semi_norm_error_p_refine_penalty_' num2str(penalty) '.eps']);





DG=figure;


for i=1:size(NO_elem,2)

errL_2 =NaN(Polynomial_degree(i),1); dof = NaN(Polynomial_degree(i),1);


for j=1:Polynomial_degree(i)
    
    
load(['Error ' num2str(NO_elem(i)) ' cube Elements penalty ' num2str(penalty) ' P' num2str(j) ' basis.mat'])


errL_2(j)=DG_err;         dof(j) = dim_FEM;

    
    
end


semilogy(dof.^(1/3),errL_2,['-h' num2str(color{i})],'LineWidth',1,'MarkerSize',8);

hold on;

end

hold off;

legend([num2str(NO_elem(1)) ' P basis'],...
       [num2str(NO_elem(2)) ' P basis'],...
       [num2str(NO_elem(3)) ' P basis'],...
       [num2str(NO_elem(4)) ' P basis']);

xlabel('Dof^{1/3}','FontSize',18);

ylabel('||| u-u_{h}|||_{DG}','FontSize',18);

set(gca,'FontSize',15);

title(['DG norm error under p-refinement penalty ' num2str(penalty) ],'FontSize',20)

saveas(DG,['Poisson_DG_norm_error_p_refine_penalty_' num2str(penalty) ],'fig');

print(DG,'-depsc',['Poisson_DG_norm_error_p_refine_penalty_' num2str(penalty) '.eps']);

