function val = basis1D (x,node,m,h,order_basis,order_derivative)



if order_basis < 2  % call nodal basis
    
    index = order_basis ;
    
    
    val = nodal_basis(x,node,index,order_derivative);
    
else % call modal basis 
    
    
    val = modal_basis(x,m,h,order_basis,order_derivative);
    
end



end