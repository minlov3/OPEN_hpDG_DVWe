function val = gradtensor_leg3(x,m,h,order)

val = NaN(size(x,1),3);


val(:,1)=shift_leg_derivative(x(:,1),m(1),h(1),order(1),1).*shift_leg_derivative(x(:,2),m(2),h(2),order(2),0).*shift_leg_derivative(x(:,3),m(3),h(3),order(3),0);

val(:,2)=shift_leg_derivative(x(:,1),m(1),h(1),order(1),0).*shift_leg_derivative(x(:,2),m(2),h(2),order(2),1).*shift_leg_derivative(x(:,3),m(3),h(3),order(3),0);

val(:,3)=shift_leg_derivative(x(:,1),m(1),h(1),order(1),0).*shift_leg_derivative(x(:,2),m(2),h(2),order(2),0).*shift_leg_derivative(x(:,3),m(3),h(3),order(3),1);


end