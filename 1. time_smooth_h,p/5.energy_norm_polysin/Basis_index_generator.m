function  Lege_ind = Basis_index_generator(Polydegree,Basis_Type,Dimension)

if Dimension ==2
   

switch Basis_Type 
    

    case 'P'
               

% P basis

t=0:Polydegree;  t=t';

temp = ones(Polydegree+1,1);

Lege_ind = [kron(t,temp) , kron(temp,t)];

% delete the index 

index = find( sum(Lege_ind,2) > Polydegree );

Lege_ind(index,:) = [];

    
    case 'Q'

% Q basis

t=0:Polydegree;  t=t';

temp = ones(Polydegree+1,1);

Lege_ind = [kron(t,temp) , kron(temp,t)];

end


end





if Dimension ==3
   

switch Basis_Type 
    

    case 'P'
               

% P basis

t=0:Polydegree;  t=t';

temp = ones(Polydegree+1,1);

Lege_ind = [  kron(kron(t,temp),temp) , kron(kron(temp,t),temp), kron(kron(temp,temp),t)];


% delete the index 

index = find( sum(Lege_ind,2) > Polydegree );

Lege_ind(index,:) = [];

    
    case 'Q'

% Q basis

t=0:Polydegree;  t=t';

temp = ones(Polydegree+1,1);

Lege_ind = [  kron(kron(t,temp),temp) , kron(kron(temp,t),temp), kron(kron(temp,temp),t)];



 case 'PQ'
               

% 2D space with P basis tensor 1D time

t=0:Polydegree;  t=t';

temp = ones(Polydegree+1,1);

Lege_ind = [  kron(kron(t,temp),temp) , kron(kron(temp,t),temp), kron(kron(temp,temp),t)];


% delete the index 

index = find( sum(Lege_ind(:,1:2),2) > Polydegree );

Lege_ind(index,:) = [];

end


end


end