function z = A_gamma(x)

if size(x,2)~= 3
    
   error('input should be 3 dimensional points')
    
end

t = x(:,3);

gamma=10;

x11 = zeros(size(x,1),1); 

x12 = zeros(size(x,1),1);

x13 = zeros(size(t,1),1);

x21 = zeros(size(x,1),1);

x22 = zeros(size(x,1),1);

x23 = zeros(size(t,1),1);

x31 = zeros(size(t,1),1);

x32 = zeros(size(t,1),1);

x33 = gamma.*ones(size(t,1),1);


% case 1


z = [x11, x12, x13,...
     x21, x22, x23,...
     x31, x32, x33];

% % case 3
% 
% 
% z = [x11, x12, x13,...
%      x21, x22, x23,...
%      x31, x32, x33];



end