%% Show the mesh



figure;
hold on;
for i = 1: NT
 %% plot element
 
b = Node(Elem{i,1},:);

%% 2 face

b1 = b(1:4,:);

b2 = b(5:8,:);

plot3([b1(:,1) ;b1(1,1)],[b1(:,2); b1(1,2)],[b1(:,3); b1(1,3)],'k-','LineWidth',1)    

plot3([b2(:,1) ;b2(1,1)],[b2(:,2); b2(1,2)],[b2(:,3); b2(1,3)],'k-','LineWidth',1)    

%% 4 edges
edge1 = b([1,5],:); edge2 = b([2,6],:);  edge3 = b([3,7],:);  edge4 = b([4,8],:);


plot3(edge1(:,1),edge1(:,2),edge1(:,3),'k-','LineWidth',1);

plot3(edge2(:,1),edge2(:,2),edge2(:,3),'k-','LineWidth',1); 

plot3(edge3(:,1),edge3(:,2),edge3(:,3),'k-','LineWidth',1); 

plot3(edge4(:,1),edge4(:,2),edge4(:,3),'k-','LineWidth',1); 

%% barycenter

bc= sum(b)./size(b,1);

text(bc(:,1),bc(:,2), bc(:,3),num2str(i),'color','r');

end; 


view(3)

axis on;

xlabel('X');ylabel('Y');zlabel('Z');