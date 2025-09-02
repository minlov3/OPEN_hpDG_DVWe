 function total_index = get_space_tensor_time_index(space_index,No_basis_time)
 
 
 No_basis_space = size(space_index,1);
 
 total_index = NaN(No_basis_space.*No_basis_time,1);
 
 
 for i=1:size(space_index,1)
     
     ind = (i-1)*No_basis_time+1:i*No_basis_time;
     
     total_index(ind) = (space_index(i)-1)*No_basis_time+1: space_index(i)*No_basis_time; 
     
 end
 
 
 
 end