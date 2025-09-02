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
zeta=1;
gamma=1;
eta=1;

u = sin( 2.*pi.*x(:,1)  ).*sin( 2.*pi.*x(:,2)  ).*t.^2;

ux = 2.*pi.*cos( 2.*pi.*x(:,1)  ).*sin( 2.*pi.*x(:,2)  ).*t.^2;
uy = 2.*pi.*sin( 2.*pi.*x(:,1)  ).*cos( 2.*pi.*x(:,2)  ).*t.^2;
ut = sin( 2.*pi.*x(:,1)  ).*sin( 2.*pi.*x(:,2)  ).*2.*t ;

uxt = 2.*pi.*cos( 2.*pi.*x(:,1)  ).*sin( 2.*pi.*x(:,2)  ).*2.*t;
uyt = 2.*pi.*sin( 2.*pi.*x(:,1)  ).*cos( 2.*pi.*x(:,2)  ).*2.*t;
utt = sin( 2.*pi.*x(:,1)  ).*sin( 2.*pi.*x(:,2)  ).*2;

uxx = -2.*pi.*2.*pi.*sin( 2.*pi.*x(:,1)  ).*sin( 2.*pi.*x(:,2)  ).*t.^2;
uyy = -2.*pi.*2.*pi.*sin( 2.*pi.*x(:,1)  ).*sin( 2.*pi.*x(:,2)  ).*t.^2;

uxxt = -2.*pi.*2.*pi.*sin( 2.*pi.*x(:,1)  ).*sin( 2.*pi.*x(:,2)  ).*2.*t;
uyyt = -2.*pi.*2.*pi.*sin( 2.*pi.*x(:,1)  ).*sin( 2.*pi.*x(:,2)  ).*2.*t;

z = utt + gamma.*ut - zeta.*(uxx+uyy) - eta.*(uxxt+uyyt);

end