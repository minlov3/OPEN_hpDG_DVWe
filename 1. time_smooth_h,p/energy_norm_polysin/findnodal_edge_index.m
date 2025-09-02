function [Dirinodal_id, Diriedge_id] = findnodal_edge_index(Diriedge_node,BDbox,FEM_index)


%% we want to evalate the edge function along x, y seperately 
   
   test_x_y = find(Diriedge_node(1,:)-Diriedge_node(2,:)~=0);
   
   switch test_x_y
      
       % along x direction
       
       case 1
           
           if Diriedge_node(1,2) == max(BDbox(:,2))
              
               Dirinodal_id = find(FEM_index(1:4,2) ==1);               
               
               Diriedge_id = find(FEM_index(5:end,2) == 1);
               
               Diriedge_id = Diriedge_id+4;
               
           else
               
                Dirinodal_id = find(FEM_index(1:4,2) ==0);               
               
                Diriedge_id = find(FEM_index(5:end,2) == 0);
               
                Diriedge_id = Diriedge_id+4;
               
           end   
                     
           
       % along y direction    
           
       case 2
           
           
           if Diriedge_node(1,1) == max(BDbox(:,1))
              
               Dirinodal_id = find(FEM_index(1:4,1) ==1);               
               
               Diriedge_id = find(FEM_index(5:end,1) == 1);
               
               Diriedge_id = Diriedge_id+4;
               
           else
               
                Dirinodal_id = find(FEM_index(1:4,1) ==0);               
               
                Diriedge_id = find(FEM_index(5:end,1) == 0);
               
                Diriedge_id = Diriedge_id+4;
               
           end              
       
   end
   
end