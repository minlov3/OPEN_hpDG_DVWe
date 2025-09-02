function z = gradient_b(x)


if size(x,2)~= 3
    
   error('input should be 2 dimensional points')
    
end


t = x(:,3);

z = zeros(size(t,1),1);

end