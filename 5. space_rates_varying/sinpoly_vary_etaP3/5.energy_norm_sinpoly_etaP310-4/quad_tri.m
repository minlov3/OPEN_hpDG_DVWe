function [P_Qpoints, weights] = quad_tri(node,order)

[w_x,x] = quad_GL(order);

[w_y,y] = quad_GJ1(order);

quad_x = kron(x,ones(size(w_y,1),1)); quad_y = kron(ones(size(w_x,1),1),y);

weights = kron(w_x,w_y) ;

shiftpoints = [(1+quad_x).*(1-quad_y).*0.5-1, quad_y ];
 
 ref_points = 0.5.*shiftpoints+0.5; 

 B = [node(1,:)-node(3,:); node(2,:)-node(3,:)];
   
    De_tri = abs(det(B));  De = De_tri./4;  Duffy_y = 0.5;
        
    phy = reference_to_physical_t3 ( node', size(ref_points,1), ref_points' );
    
    P_Qpoints  = phy'; 

weights = weights.*De.*Duffy_y;

end