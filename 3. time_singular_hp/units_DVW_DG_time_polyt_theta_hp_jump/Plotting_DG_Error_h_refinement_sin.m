%clear all; close all;

%% different elements we use different order of basis

NO_Elements = [16,64,256,1024,4096];    color = {'b','r','g','k'};



DG_err=figure;

slope_Poly = NaN(1,4);

order= 4;

Basis_Type  = 'Se';  % Q

Norm_type = 'L2_H1';  % L2_L2, L2_H1, L_inf_L2, L_inf_H1



for k = 1:order


%% plotting errors for the L2 norm for p refinement

%%P basis


Polynomial_degree = k;    penalty=10;



NO_elem = NO_Elements;   


err_P = NaN(length(NO_elem),1); dof_P = NaN(length(NO_elem),1);

for i=1 :length(NO_elem)
    

     load(['Error ' num2str(NO_elem(i)) ' rectangle Elements for FEM '  Basis_Type num2str(Polynomial_degree) ' basis.mat'])

    
     switch  Norm_type
         
     case 'L2_L2'
    
        err_P(i,1)=L2_L2_err;      
    
     case 'L2_H1'
    
        err_P(i,1)=L2_H1_err;      

     case 'L_inf_L2'

        err_P(i,1)= L_inf_L2_err;
        
     case 'L_inf_H1'

        err_P(i,1)= L_inf_H1_err;   
        
     end

dof_P(i) = dim_FEM;

       

end




logerr_P1 = abs(log(err_P(:,1))); 



slope_P1 = abs((logerr_P1(2:end)-logerr_P1(1:end-1))./(log(dof_P(2:end).^(1./3))-log(dof_P(1:end-1).^(1./3))));



slope_Poly(k) = mean(slope_P1(end-1:end));

%slope_Poly(k) = max(slope_P1);

%{

switch  Mesh_type
    
     case 'poly'

loglog(dof_P.^(1./3),err_P(:,1),['-h' num2str(color{k})],'LineWidth',1.5,'MarkerSize',10);

     case 'rect'

loglog(dof_P.^(1./3),err_P(:,1),['-s' num2str(color{k})],'LineWidth',1.5,'MarkerSize',10);

end

%}


switch  Basis_Type
         
    
     case 'P'

loglog(dof_P.^(1./3),err_P(:,1),['-.s' num2str(color{k})],'LineWidth',1.5,'MarkerSize',10);

     case 'PQ'

loglog(dof_P.^(1./3),err_P(:,1),['-*' num2str(color{k})],'LineWidth',1.5,'MarkerSize',10);


     case 'Q'

loglog(dof_P.^(1./3),err_P(:,1),['-o' num2str(color{k})],'LineWidth',1.5,'MarkerSize',10);

end

hold on;

end



legend(['FEM rect ' Basis_Type '1 slope ' num2str(slope_Poly(1))]...
      ,['FEM rect ' Basis_Type '2 slope ' num2str(slope_Poly(2))]...
      ,['FEM rect ' Basis_Type '3 slope ' num2str(slope_Poly(3))]...
      ,['FEM rect ' Basis_Type '4 slope ' num2str(slope_Poly(4))]...
      ,'Location','SouthWest')





xlabel('Dof^{1/3}','FontSize',18);



switch  Norm_type
    
    case 'L2_L2'
    
        ylabel('|| u-u_{h}||_{L^2((0,T),L^2(\Omega))}','FontSize',18);
        
        title('L^2((0,T),L^2(\Omega)) norm error under h-refinement','FontSize',18)

    
     case 'L2_H1'
    
        ylabel('|| u-u_{h}||_{L^2((0,T),H^1(\Omega))}','FontSize',18);
        
        title('L^2((0,T),H^1(\Omega)) norm error under h-refinement','FontSize',18)

     case 'L_inf_L2'

        ylabel('|| u-u_{h}||_{L^\infty((0,T),L^2(\Omega))}','FontSize',18);
        
        title('L^\infty((0,T),L^2(\Omega)) norm error under h-refinement','FontSize',18)
        
     case 'L_inf_H1'

        ylabel('|| u-u_{h}||_{L^\infty((0,T),H^1(\Omega))}','FontSize',18);
        
        title('L^\infty((0,T),H^1(\Omega)) norm error under h-refinement','FontSize',18)
        
end



set(gca,'FontSize',15)


