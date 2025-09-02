function val = modal_basis(x,m,h,order_basis,order_derivative)



if order_basis < 2
   
    error('modal basis is only allowed for order_basis greater than 2');
    
end



switch order_derivative
   
    % standard modal basis order k is the antiderivative of shifted legengre polynomial
    % at order k-1
    
    case 0
          
        coe = -1./(order_basis.*(order_basis-1)).*(1 - ((x - m)./h ).^2).*h^2;
      
        val = coe.*shift_leg_derivative(x,m,h,order_basis-1,1);
        
        
    % derivative of modal basis order k is the shifted legengre polynomial
    % at order k-1
    
    case 1 
    
    
        val = shift_leg_derivative(x,m,h,order_basis-1,0);
        
        
    
end





end



