function FEM_index = FEM_Basis_index2D(Polydegree,Basis_Type)


%% Q basis Nodal index  4 in 2D

% 0 is the nodal basis interpolating on the xmin or ymin

Nodal_index = [0 ,0; ...
               1 ,0; ...
               1 ,1; ...
               0 ,1];
if Polydegree == 1
    
    FEM_index = Nodal_index ;
   
    
else
%% edge index  4*(p-1) the orientation is 
% P4(0,1) <----- P3(1,1)
%   |       E3       |
% E4|                |E2
%   |                |
% P1(0,0) -----> P2(1,0)   
%           E1

I = 2:Polydegree; I=I';

E1_index = [I, zeros(size(I))];

E2_index = [ones(size(I)),  I];

E3_index = [I,  ones(size(I))];

E4_index = [zeros(size(I)), I];

Edge_index = [E1_index ;E2_index; E3_index; E4_index];
    
    
%% Modal index (p-1)^2

t=0:Polydegree;  t=t';

temp = ones(Polydegree+1,1);

auxilary_index = [kron(t,temp) , kron(temp,t)];

ind1 = find(auxilary_index(:,1) > 1);

ind2 = find(auxilary_index(ind1,2)> 1);

Modal_index = auxilary_index(ind1(ind2),:);

switch Basis_Type
    
    case 'Q'

%% The FEM index

FEM_index = [Nodal_index ;Edge_index;Modal_index];

    case 'Se'        

%reduce the high moment 

Sind = find(sum(Modal_index,2)<=Polydegree);

%% The FEM index

SModal_index = Modal_index(Sind,:);

FEM_index = [Nodal_index ;Edge_index;SModal_index];
end

end