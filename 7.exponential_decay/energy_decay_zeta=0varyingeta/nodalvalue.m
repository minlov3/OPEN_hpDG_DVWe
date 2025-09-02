function [u_DG_val, space_m ] = nodalvalue(Space_Node,Time_Node,BDbox , coef ,Lege_ind)


% information about the bounding box

h = (BDbox(2,:)-BDbox(1,:))./2;  

m = 0.5.*sum(BDbox);  space_m = m(1:2);

dim_elem =size(Lege_ind,1); % number of basis for each element.



 P_Qpoints = [m(1:2), Time_Node(2)];
    
 

    
    % construct the matrix for all the local basis function
    
     P = zeros(size(P_Qpoints,1) ,dim_elem);

    
    for i =1:dim_elem
        
        P(:,i)= tensor_leg3(P_Qpoints,m,h,Lege_ind(i,:));
        
        
    end
  
    u_DG_val = P*coef;   %DG solution;
    
  
 
 
 
     
end







