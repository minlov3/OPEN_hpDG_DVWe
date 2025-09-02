function BDbox = Bounding_box(Node)


if size(Node,2) ==2


bdbox_x = sort(Node(:,1)); bdbox_y = sort(Node(:,2)); 
    
BDbox =[bdbox_x([1,end]) ,  bdbox_y([1,end]) ];

end


if size(Node,2) ==3


bdbox_x = sort(Node(:,1)); bdbox_y = sort(Node(:,2)); bdbox_z = sort(Node(:,3));
    
BDbox =[bdbox_x([1,end]) ,  bdbox_y([1,end]) ,bdbox_z([1,end])];

end

end