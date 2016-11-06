function A_block = calc_transition_matrix(A,p,Order)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%
    % Create a block diagonal matrix with i-th block having parameters of
    % the ith ROI which is
    % ai or  [ai(1) ai(2)]   for AR(1) and AR(2) respectively
    %        [1      0   ]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch Order
    case 1
        if isequal(size(A,1),size(A,2))
        A_block = A;
        else
            A_block = diag(A);
        end
        
    case 2
    Elements = A';
    % diagonal elements
    D = sparse(1:2*p,1:2*p,upsample(Elements(1,:),2),2*p,2*p);
    % upper diagonal elements
    r = upsample(Elements(2,:),2);
    U = sparse(1:2*p-1,2:2*p,r(1:end-1),2*p,2*p);
    % lower diagonal elements
    c = upsample(ones(1,p),2);
    L = sparse(2:2*p,1:2*p-1,c(1:end-1),2*p,2*p);
    % tridiagonal 
    A_block = D+L+U;

end

end