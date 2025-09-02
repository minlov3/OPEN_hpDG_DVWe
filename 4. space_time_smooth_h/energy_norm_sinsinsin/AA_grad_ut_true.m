function z = AA_grad_ut_true(x)


if size(x,2)~= 3
    
   error('input should be 2 dimensional points')
    
end

t = x(:,3);

% % case 1 smooth
% 
% z = [pi*cos(pi.*x(:,1)).*sin(pi.*x(:,2)).*exp(-2*pi^2.*t),...
%      pi*sin(pi.*x(:,1)).*cos(pi.*x(:,2)).*exp(-2*pi^2.*t),...
%      -2*pi.^2*sin(pi.*x(:,1)).*sin(pi.*x(:,2)).*exp(-2*pi^2.*t)];


% % case 2 singular at t =0
% 
% alpha = 3/4; 
% 
% s1 = t.^(alpha).*(1-2.*x(:,1)).*x(:,2).*(1-x(:,2));
% 
% s2 = t.^(alpha).*x(:,1).*(1-x(:,1)).*(1-2.*x(:,2));
% 
% s3 = alpha.*t.^(alpha-1).*x(:,1).*(1-x(:,1)).*x(:,2).*(1-x(:,2));
% 
% 
% z = [s1 , s2 , s3];


% case 3 oscillatory

u = sin(pi.*x(:,1)).*sin(pi.*x(:,2)).*sin(pi.*t);

ux = pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2)).*sin(pi.*t);
uy = pi.*sin(pi.*x(:,1)).*cos(pi.*x(:,2)).*sin(pi.*t);
ut = pi.*sin(pi.*x(:,1)).*sin(pi.*x(:,2)).*cos(pi.*t);

uxt = pi.^2.*cos(pi.*x(:,1)).*sin(pi.*x(:,2)).*cos(pi.*t);
uyt = pi.^2.*sin(pi.*x(:,1)).*cos(pi.*x(:,2)).*cos(pi.*t);
utt = -pi.^2.*sin(pi.*x(:,1)).*sin(pi.*x(:,2)).*sin(pi.*t);

z =  [uxt , uyt , utt];

end