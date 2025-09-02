function val = tensor_leg3(x,m,h,order)

%val=shift_leg(x(:,1),m(1),h(1),order(1)).*shift_leg(x(:,2),m(2),h(2),order(2)).*shift_leg(x(:,3),m(3),h(3),order(3));


val=shift_leg_derivative(x(:,1),m(1),h(1),order(1),0).*shift_leg_derivative(x(:,2),m(2),h(2),order(2),0).*shift_leg_derivative(x(:,3),m(3),h(3),order(3),0);


end