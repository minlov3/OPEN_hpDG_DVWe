function val = nodal_basis(x,node,index,order_derivative)




if index >=2
   
    error('nodal basis is only allowed for index = 0 or 1');
    
end

% index = 0 means the nodal basis = 1|x_min, = 0|x_max    
% index = 1 means the nodal basis = 0|x_min  = 1|x_max 

x_min = min(node);  x_max = max(node);

switch order_derivative
   
    % standard nodal interpolant
    
    case 0
        
        if index ==0
           
            val = (x - x_max)./(x_min - x_max);
            
        else
            
            val = (x - x_min)./(x_max - x_min);
            
        end
        
    
        
        
    % derivative of nodal interpolant
    
    case 1 
    
        
         if index ==0
           
            val = 1./(x_min - x_max).*ones(size(x));
            
        else
            
            val = 1./(x_max - x_min).*ones(size(x));
            
        end
        
        
        
    
end





end



