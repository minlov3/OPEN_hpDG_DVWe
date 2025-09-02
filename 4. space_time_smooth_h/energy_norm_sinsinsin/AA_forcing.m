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

u = sin(pi.*x(:,1)).*sin(pi.*x(:,2)).*sin(pi.*t);

ux = pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2)).*sin(pi.*t);
uy = pi.*sin(pi.*x(:,1)).*cos(pi.*x(:,2)).*sin(pi.*t);
ut = pi.*sin(pi.*x(:,1)).*sin(pi.*x(:,2)).*cos(pi.*t);

uxt = pi.^2.*cos(pi.*x(:,1)).*sin(pi.*x(:,2)).*cos(pi.*t);
uyt = pi.^2.*sin(pi.*x(:,1)).*cos(pi.*x(:,2)).*cos(pi.*t);
utt = -pi.^2.*sin(pi.*x(:,1)).*sin(pi.*x(:,2)).*sin(pi.*t);

uxx = -pi.^2.*sin(pi.*x(:,1)).*sin(pi.*x(:,2)).*sin(pi.*t);
uyy = -pi.^2.*sin(pi.*x(:,1)).*sin(pi.*x(:,2)).*sin(pi.*t);

uxxt = -pi.^3.*sin(pi.*x(:,1)).*sin(pi.*x(:,2)).*cos(pi.*t);
uyyt = -pi.^3.*sin(pi.*x(:,1)).*sin(pi.*x(:,2)).*cos(pi.*t);


z = utt + gamma.*ut - zeta.*(uxx+uyy) - eta.*(uxxt+uyyt);

end