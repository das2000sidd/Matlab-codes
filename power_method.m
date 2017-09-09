A=[1 2 2 ; 2 4 5 ; 1 2 3];
x=[0 0 1].';
%%EigVectorMat = zeros(3,50);
eigvectors=zeros(3,50);
%%eigvectors(:,1) = x';
eigvalues = zeros(50,1);
%% doing dominant eigen value and eigen vec calculation
for c = 1:50
    x1=A*x;
    %%EigVectorMat(:,c) = x1;
    eigvectors(:,c) = x1./max(x1);
    %%eigvectors(:,c) = x1./x1(1,1);
    eigvalues(c,:) = max(x1);
    x=eigvectors(:,c);
end
dominant_eigenvector = eigvectors(:,4); %% stabilizes after 4th iteration A^4.x
dominant_eigenvalue = eigvalues(6,:);
%%normalisedeigvec 
normalisedEigVec = dominant_eigenvector./norm(dominant_eigenvector);

[V,D] = eig(A);
eigvec = V(1,:);
%%eigvec_new = eigvec.';
eigval = D(1,1);
distance1 = norm(dominant_eigenvector-eigvec.'); %% same as line 25
distance_no_norm_formula = sqrt(sum((dominant_eigenvector-eigvec.').^2));
negative_dominant_eigenvector = -(dominant_eigenvector);
distance2 = norm(negative_dominant_eigenvector-eigvec);
eigval_diff = dominant_eigenvalue-eigval;

   %% CORRECT SOLUTION 
    
    