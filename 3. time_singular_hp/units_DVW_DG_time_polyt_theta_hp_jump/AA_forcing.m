function z= AA_forcing(x)


if size(x,2)~= 3
    
   error('input should be 2 dimensional points')
    
end




t = x(:,3);


% % case 1 smooth
% 
% z = zeros(size(t,1),1);


% % case 2 singular at t =0
% 
% alpha = 3/4; 
% 
% s1 = 2.*t.^(alpha).*(x(:,1).*(1-x(:,1)) + x(:,2).*(1-x(:,2)));
% 
% s2 = alpha.*t.^(alpha-1).*x(:,1).*(1-x(:,1)).*x(:,2).*(1-x(:,2));
% 
% z = s1+s2;


% case 3 oscillatory

zeta = 1;
gamma = 1;
eta = 1;

theta=3.51;
u = t.^theta.*x(:,1).*(  1 - x(:,1)  ).*x(:,2).*( 1 - x(:,2) );

ux = t.^theta.*(  1 - 2.*x(:,1)  ).*x(:,2).*( 1 - x(:,2) );
uy = t.^theta.*x(:,1).*(  1 - x(:,1)  ).*( 1 - 2.*x(:,2) );
ut = theta.*t.^(theta-1).*x(:,1).*(  1 - x(:,1)  ).*x(:,2).*( 1 - x(:,2) );

uxt = theta.*t.^(theta-1).*(  1 - 2.*x(:,1)  ).*x(:,2).*( 1 - x(:,2) );
uyt = theta.*t.^(theta-1).*x(:,1).*(  1 - x(:,1)  ).*( 1 - 2.*x(:,2) );
utt = theta.*(theta-1).*t.^(theta-2).*x(:,1).*(  1 - x(:,1)  ).*x(:,2).*( 1 - x(:,2) );
 
uxx = t.^theta.*(  -2  ).*x(:,2).*( 1 - x(:,2) );
uyy = t.^theta.*x(:,1).*(  1 - x(:,1)  ).*(  -2  ); 

uxxt = theta.*t.^(theta-1).*(  -2  ).*x(:,2).*( 1 - x(:,2) );
uyyt = theta.*t.^(theta-1).*x(:,1).*(  1 - x(:,1)  ).*(  -2  );  

z = utt + gamma.*ut - zeta.*(uxx+uyy) - eta.*(uxxt+uyyt);

end