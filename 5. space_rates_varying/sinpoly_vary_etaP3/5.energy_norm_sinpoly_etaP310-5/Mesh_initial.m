%%Initialize the Polymesh generator to construct the information about the
%%underlying FEM space and the geometry of the computational domain.



clear all; 


%[Node,Element,Supp,Load,P]= PolyMesher(@MbbDomain,4096,200);


%% Generating the squ are domain [0,1]^2 with ranctangel elements.
  dx =1/2; dy =1/2;
  
  [X,Y] = meshgrid(0+dx/2:dx:1,0+dy/2:dy:1);
  
  
  P = [X(:) Y(:)];
 
  %% grading mesh
  
  % dist_x=(P(:,1)+1)./2;    sigma=0.6;
  
  % P(:,1)= dist_x.^sigma.*(P(:,1)+1)-1;
  
  
  
  
 %% Lshape  
 % [i,~,~] = find(P(:,1)>0);
  
 % [j,~,~] = find(P(:,2)<0);
  
 % index = intersect(i,j,'rows');
  
 % true_P = P;   true_P(index,:)=[];
  
 

  [Node,Element,Supp,Load,P]= PolyMesher(@MbbDomain,25,0,P);
  
  

  
%% Generate bounding box for each polygonal element

% move to the correct place.

n = 8;   Node = round(10^n.*Node)/10^n;



% anisotropic mesh

%[ix] = find(Node(:,1)~=0);  [iy] = find(Node(:,2)~=0);

%scale_x = 2.^(Node(ix,1)./dx);  scale_y = 2.^(Node(iy,2)./dy);

%Node(ix,1) = scale_x./max(scale_x);  Node(iy,2) = scale_y./max(scale_y);


%[iy] = find(Node(:,2)~=0);

%scale_y = 2.^(Node(iy,2)./dy);
 
%Node(iy,2) = scale_y./max(scale_y);







% radial meshes  I

%scale=max(abs(Node),[],2);

%Node= [Node(:,1).*scale, Node(:,2).*scale];


% radial meshes  II

%scale=max(abs(Node),[],2);

%Node= [Node(:,1).*scale.^2, Node(:,2).*scale.^2];


%  and  finding all the edges. Label element to edges 



N = size(Node,1); NT = size(Element,1); 

Elem = cell(NT,2);

totalEdge = NaN(1,2);  elemperedge = NaN(NT,1);

for i =1: NT
    
    % calculate all the bounding box
    
    Elem{i,1}= Element{i}';
   
    % barycenter
    
    P(i,:) = sum(Node(Elem{i,1},:))./size(Node(Elem{i,1},:),1);
    
    
    bdbox_x = sort(Node(Elem{i,1},1));
    
    bdbox_y = sort(Node(Elem{i,1},2));
    
    Elem{i,2} = [bdbox_x([1,end]) ,  bdbox_y([1,end])]; %bounding box is [x_min y_min
                                                        %                 x_max y_max]
                                                 
   % find all the local index of triangle in polygon!
      
    Elem{i,3} = delaunay(Node(Elem{i,1},:)) ; 
    
%     % calculate the measure of polygons based on subtriangulation
%     
%         ve = zeros(size(tri,1),2,3);
%         
%         ve(:,:,1) = Node(tri(:,3),:)-Node(tri(:,2),:);
%         ve(:,:,2) = Node(tri(:,1),:)-Node(tri(:,3),:);
%         ve(:,:,3) = Node(tri(:,2),:)-Node(tri(:,1),:);
%         
%         area = 0.5*abs(-ve(:,1,3).*ve(:,2,2)+ve(:,2,3).*ve(:,1,2));
%              
%         
%     Elem{i,4}=sum(area);
    
    
% %      % local edge and the outward normal vector!      
% %     
% %     localedge = [Elem{i,1}, [Elem{i,1}(2:end); Elem{i,1}(1) ] ];
% %     
% %     localedge = sort(localedge,2);
% %     
% %     Elem{i,5} =  localedge;
% %     
% %     % tangert vector of the localedge%%%%%%
% % localedgetan_vec = Node(localedge(:,1),:) - Node(localedge(:,2),:);
% %  
% % 
% %  
% %   % unit normal vector of the localedge%%%%%%
% % localedgenorvec = [-localedgetan_vec(:,2)./sqrt(localedgetan_vec(:,1).^2 +localedgetan_vec(:,2).^2) , localedgetan_vec(:,1)./sqrt(localedgetan_vec(:,1).^2 +localedgetan_vec(:,2).^2)];
% % 
% % 
% % % outward normal vector of the edge%%%%%%
% %   
% %   internalP =  kron(P(i,:),ones(size(localedge,1),1));
% % 
% %   localoutward = Node(localedge(:,1),:)-internalP ;
% %   
% %  
% %  localindex =  max(sum(localedgenorvec.*localoutward,2),0);
% %  
% % 
% %  [k ,j ,s] = find( localindex==0);  
% %  
% %  localedgenorvec(k,:) = - localedgenorvec(k,:);
% %     
% %  Elem{i,6} = localedgenorvec;          
    
                                                        
    % take out all the edges  
    
    totalEdge = [totalEdge; [Elem{i,1}, [Elem{i,1}(2:end); Elem{i,1}(1) ] ] ] ;                                                       
                    
    % label edges per element
    
    elemperedge(i) = size(Elem{i,1},1);
end

elemtotaledge = cumsum(elemperedge);

totalEdge(1,:)=[];



%ploting the bounding box to cover  element!

if NT <= 300

figure;
hold on;
for i = 1: NT
 %% plot element
 
b = Node(Elem{i,1},:);
plot([b(:,1) ;b(1,1)],[b(:,2); b(1,2)],'k-','LineWidth',1)    
 %% bounding box 
 
a=[Elem{i,2}(1,:);[Elem{i,2}(1,1) Elem{i,2}(2,2)]  ;Elem{i,2}(2,:);  [Elem{i,2}(2,1) Elem{i,2}(1,2)] ] ;
plot([a(:,1) ;a(1,1)],[a(:,2); a(1,2)],'o-','LineWidth',1,'MarkerSize',5)

%text(P(i,1),P(i,2), num2str(i));

end; 


end;
%% Classify all the edges

% totalEdge is all the edges from summation of each element
% edge is the all the edges from triagulation 
% bdEdge is the boundary edge and intEdge is the interior edge

totalEdge = sort(totalEdge, 2);

[i , j ,s ] = find(sparse(totalEdge(:,2),totalEdge(:,1),1));

edge = [j,i]; bdEdge = [j(s==1), i(s==1)]; intEdge = [j(s==2), i(s==2)];

%% The relation between edge to element and

% internal edge is shared by 2 elements and boundary edge is only used by 1
% edge

int_edge2elem = NaN(size(intEdge,1),2); bd_edge2elem = NaN(size(bdEdge,1),1);

for i = 1: size(intEdge,1)
   
    edge2elem = elem_share_edge(intEdge(i,:),totalEdge,elemtotaledge);
    
    int_edge2elem(i,:) =   edge2elem';
    
end

for i = 1: size(bdEdge,1)
   
    edge2elem = elem_share_edge(bdEdge(i,:),totalEdge,elemtotaledge);
    
    bd_edge2elem(i) =   edge2elem;
    
end


%% Plot for bdedge


% figure; hold on;
% 
% for i = 1 : size(bdEdge,1) 
%     
%    
%    
%   t=Node(bdEdge(i,:),:);
%    
%   plot(t(:,1), t(:,2),'o-','LineWidth',1,'MarkerSize',5)
%        
% end



%% Outward Normal vectors for all boundary edges

% For bdedge, only one normal vector
% For internal edge there two normnal vector, but we only use one.

% tangert vector of the edge%%%%%%
 bdtan_vec = Node(bdEdge(:,1),:) - Node(bdEdge(:,2),:);
 
inttan_vec = Node(intEdge(:,1),:) - Node(intEdge(:,2),:);
 
  % normal vector of the edge%%%%%%
 bdnorvec = [bdtan_vec(:,2)./sqrt(bdtan_vec(:,1).^2 +bdtan_vec(:,2).^2) , -bdtan_vec(:,1)./sqrt(bdtan_vec(:,1).^2 +bdtan_vec(:,2).^2)];

intnorvec = [inttan_vec(:,2)./sqrt(inttan_vec(:,1).^2 +inttan_vec(:,2).^2) , -inttan_vec(:,1)./sqrt(inttan_vec(:,1).^2 +inttan_vec(:,2).^2)];
 
% outward normal vector of the edge%%%%%%
  bdoutward = Node(bdEdge(:,1),:)-P(bd_edge2elem(:),:);
  
  intoutward = Node(intEdge(:,1),:)-P(int_edge2elem(:,1),:); % the first element
 
  
 
 bdindex =  max(sum(bdnorvec.*bdoutward,2),0);
 
 intindex =  max(sum(intnorvec.*intoutward,2),0);
 
 [i ,j ,s] = find(bdindex==0);   [m ,n ,k] = find(intindex==0);
 
 bdnorvec(i,:) = - bdnorvec(i,:);
  
 intnorvec(m,:) = - intnorvec(m,:);

 

 %% Data for Oliver!

%elements = Element; vertices = Node; boundaryVertices=unique(bdEdge(:));


%savefile = [ num2str(NT) ' rectangle Elements for Oliver.mat'];

%save(savefile, 'elements','vertices','boundaryVertices','-v7.3');

 


%save all the useful data

%savefile = [ num2str(NT) ' polygonal Elements.mat'];

savefile = [ num2str(NT) ' rectangle Elements.mat'];

%savefile = [ num2str(NT) ' anisotropic rectangle Elements.mat'];
save(savefile, 'Elem','Node','N','NT','bdEdge','bd_edge2elem','intEdge','int_edge2elem','bdnorvec','intnorvec','-v7.3');

%%%%%%%%%%%%%%%%%% The intialization of Mesh is finished %%%%%%%%%%%%%%%%%%%%%

