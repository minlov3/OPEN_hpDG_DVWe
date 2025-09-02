function u = AA_u_true(x)


if size(x,2)~= 3
    
   error('input should be 2 dimensional points')
    
end


t = x(:,3);


% % case 1 smooth
% 
% 
% z = sin(pi.*x(:,1)).*sin(pi.*x(:,2)).*exp(-2*pi^2.*t);


% % case 2 singular at t =0
% 
% alpha = 3/4; 
% 
% z = t.^(alpha).*x(:,1).*(1-x(:,1)).*x(:,2).*(1-x(:,2));

theta=3.51;
u = t.^theta.*x(:,1).*(  1 - x(:,1)  ).*x(:,2).*( 1 - x(:,2) );

end