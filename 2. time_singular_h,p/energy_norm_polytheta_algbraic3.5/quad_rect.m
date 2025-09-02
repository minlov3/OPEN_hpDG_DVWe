function [P_Qpoints, weights] = quad_rect(node,order)

[w_x,x] = quad_GL(order);

[w_y,y] = quad_GL(order);


weights = kron(w_x,w_y) ;


 B = [node(1,:)-node(3,:); node(2,:)-node(3,:)];
   
    De_tri = abs(det(B));  De = De_tri./4; 
        
    x_P_Qpoints = x.*(max(node(:,1))-min(node(:,1)))*0.5+ 0.5*sum(max(node(:,1))+min(node(:,1)));
    
    y_P_Qpoints = y.*(max(node(:,2))-min(node(:,2)))*0.5+ 0.5*sum(max(node(:,2))+min(node(:,2)));    
    
    quad_x = kron(x_P_Qpoints,ones(size(w_y,1),1)); quad_y = kron(ones(size(w_x,1),1),y_P_Qpoints);


P_Qpoints = [quad_x,quad_y];


weights = weights.*De;

end