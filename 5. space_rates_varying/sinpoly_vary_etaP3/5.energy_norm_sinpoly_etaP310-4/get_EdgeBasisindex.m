function  global_EdgeBasis_ind = get_EdgeBasisindex(local_edge,TotalEdge,Polydegree)   

%% find all the global index for all the edge basis function 

  sort_edge = sort(local_edge,2);
  
  Edge_ind = NaN(size(sort_edge,1),1);
  
  global_EdgeBasis_ind = NaN(size(sort_edge,1).*(Polydegree-1),1);
    
  for i=1:size(sort_edge,1)
     
      [~,ia,~]=intersect(TotalEdge,sort_edge(i,:),'rows');
      
      Edge_ind(i) = ia;
      
      ind = (i-1)*(Polydegree-1)+1:i*(Polydegree-1);
      
      global_EdgeBasis_ind(ind) = (Edge_ind(i)-1)*(Polydegree-1)+1:Edge_ind(i)*(Polydegree-1);
      
  end
       
       


end