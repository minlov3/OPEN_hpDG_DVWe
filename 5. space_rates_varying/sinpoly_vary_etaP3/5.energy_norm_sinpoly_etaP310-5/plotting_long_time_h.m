%clear all; close all;

%% different elements we use different order of basis

NO_Elements = 1024;    color = {'b','r','g','k','m'};

time_step=[400,800, 1600,3200];


DG_err=figure;

slope_Poly = NaN(1,4);

order=4;

Basis_Type  = 'Q'; 

Norm_type = 'L_inf_L2';  % L2_H1, L_inf_L2, L_inf_H1




for k = 1:order


%% plotting errors for the L2 norm for p refinement

%%P basis


Polynomial_degree = k;    penalty=10;

NO_elem = NO_Elements;   

err_P = NaN(length(time_step),1); dof_P = NaN(length(time_step),1);


for i=1 :length(time_step)
    
     
     load(['Error ' num2str(NO_elem) ' rectangle Elements time step ' num2str(time_step(i)) ' for FEM ' Basis_Type num2str(Polynomial_degree) ' basis.mat'])


     switch  Norm_type
    
     case 'L2_H1'
    
        err_P(i,1)=L2_H1_err;      

     case 'L_inf_L2'

        err_P(i,1)= L_inf_L2_err;
        
     case 'L_inf_H1'

        err_P(i,1)= L_inf_H1_err;   
        
     end

dof_P(i) = dim_FEM;


end




%slope for poly


logerr_P1 = abs(log(err_P(:,1))); 



slope_P1 = abs((logerr_P1(2:end)-logerr_P1(1:end-1))./((log(dof_P(2:end))-log(dof_P(1:end-1)))./1));

slope_P1 = sort(slope_P1);

%slope_Poly(k) = mean(slope_P1(end-1:end));

slope_Poly(k) = max(slope_P1);



  


loglog(dof_P,err_P(:,1),['-o' num2str(color{k})],'LineWidth',1.5,'MarkerSize',10);





hold on;

end



legend(['FEM rect ' Basis_Type '1 slope ' num2str(slope_Poly(1))]...
      ,['FEM rect ' Basis_Type '2 slope ' num2str(slope_Poly(2))]...
      ,['FEM rect ' Basis_Type '3 slope ' num2str(slope_Poly(3))]...
      ,['FEM rect ' Basis_Type '4 slope ' num2str(slope_Poly(4))]...
      ,'Location','SouthWest')
     



xlabel('Dof','FontSize',18);



switch  Norm_type
    
     case 'L2_H1'
    
        ylabel('|| u-u_{h}||_{L^2((0,T),H^1(\Omega))}','FontSize',18);
        
        title('L^2((0,T),H^1(\Omega)) norm error under \tau-refinement','FontSize',18)

     case 'L_inf_L2'

        ylabel('|| u-u_{h}||_{L^\infty((0,T),L^2(\Omega))}','FontSize',18);
        
        title('L^\infty((0,T),L^2(\Omega)) norm error under \tau-refinement','FontSize',18)
        
     case 'L_inf_H1'

        ylabel('|| u-u_{h}||_{L^\infty((0,T),H^1(\Omega))}','FontSize',18);
        
        title('L^\infty((0,T),H^1(\Omega)) norm error under \tau-refinement','FontSize',18)
        
end



set(gca,'FontSize',15)


