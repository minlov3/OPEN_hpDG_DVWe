function val = FEM2D_DG_basis(x,BDbox,m,h,order_basis)

% node value is passed by bounding box

t=x(:,3);

val=basis1D(x(:,1),BDbox(:,1),m(1),h(1),order_basis(1),0).*basis1D(x(:,2),BDbox(:,2),m(2),h(2),order_basis(2),0).*shift_leg_derivative(t,m(3),h(3),order_basis(3),0);
    
end