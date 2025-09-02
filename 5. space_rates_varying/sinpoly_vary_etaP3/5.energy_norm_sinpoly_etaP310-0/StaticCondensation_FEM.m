function  U = StaticCondensation_FEM(A,F,NT,Size_Modal)


Modal_index = size(A,1)-NT*Size_Modal+1 :size(A,1);

Bd_index = 1:size(A,1)-NT*Size_Modal;

% subdivide the large matrix A into 4 blocks!

MB = A(Bd_index,Bd_index); MBI = A(Bd_index,Modal_index);

MI = A(Modal_index,Modal_index); MIB = A(Modal_index,Bd_index);

FB = F(Bd_index);  FI = F(Modal_index);

U =NaN(size(F));

clear A F;

%% inverse the MI matrix elementwise for each Block on diagonal


i =zeros(Size_Modal.^2,NT ); j =zeros(Size_Modal.^2,NT ); s =zeros(Size_Modal.^2,NT );


parfor t =1:NT
    
    ind = (t-1)*Size_Modal+1:1:t*Size_Modal;
    
    ind = ind';
    
    local_inv_block_diag = inv(MI(ind,ind));
        
      
    i(:,t) = kron(ones(Size_Modal,1),ind )  ;
   
    j(:,t) = kron(ind , ones(Size_Modal,1)) ;
   
    s(:,t) = local_inv_block_diag(:);
    
end

inv_MI = sparse(i,j ,s ,size(MI,1),size(MI,1) );


% Calculating the boundary part
B = MB-MBI*inv_MI*MIB;

L = FB - MBI*inv_MI*FI;

% reordering

P = symrcm(B);

u =B(P,P)\L(P);

%spy(B(P,P))

u(P) = u;


U(Bd_index) =u;

% Calculating the interior part

U(Modal_index) = inv_MI*FI -inv_MI*MIB*U(Bd_index);

end