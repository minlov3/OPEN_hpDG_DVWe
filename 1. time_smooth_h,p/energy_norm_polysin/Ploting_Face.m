% plotting Face


figure(1); hold on;


face = inflowBDFace;

for t = 1:size(face,1) 

order = [1,2,4,3]';    
    
t1 = Node(face(t,order),:);

plot3([t1(:,1) ;t1(1,1)],[t1(:,2); t1(1,2)],[t1(:,3); t1(1,3)],'r*-','LineWidth',1)    

xlim([0 ,1]); ylim([0 ,1]); zlim([0 ,1]);


end

view(3)

axis on;

xlabel('X');ylabel('Y');zlabel('Z');

title('inflow Face','FontSize',18);



figure(2); hold on;

face = DiribdFace;

%face = outflowFace;

for t = 1:size(face,1) 

order = [1,2,4,3]';    
    
t2 = Node(face(t,order),:);

plot3([t2(:,1) ;t2(1,1)],[t2(:,2); t2(1,2)],[t2(:,3); t2(1,3)],'k*-','LineWidth',1)    

xlim([0 ,1]); ylim([0 ,1]); zlim([0 ,1]);


end

view(3)

axis on;

xlabel('X');ylabel('Y');zlabel('Z');

title('Dirichlet boundary Face','FontSize',18);




figure(3); hold on;


face = intFace;



for t = 1:size(face,1) 

order = [1,2,4,3]';    
    
t3 = Node(face(t,order),:);

plot3([t3(:,1) ;t3(1,1)],[t3(:,2); t3(1,2)],[t3(:,3); t3(1,3)],'*-','LineWidth',1)    

xlim([0 ,1]); ylim([0 ,1]); zlim([0 ,1]);


end

view(3)

axis on;

xlabel('X');ylabel('Y');zlabel('Z');


title('Interior Face','FontSize',18);
