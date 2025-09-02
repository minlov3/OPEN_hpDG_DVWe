function z = adv(x)

if size(x,2)~= 3
    
   error('input should be 3 dimensional points')
    
end

t = x(:,3);


z = [ zeros(size(x,1),1) , zeros(size(x,1),1) , ones(size(t,1),1) ];




end
