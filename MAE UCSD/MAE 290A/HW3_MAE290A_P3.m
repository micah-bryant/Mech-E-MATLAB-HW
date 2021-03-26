number = 5;
uNorm = zeros(number,1); vNorm = zeros(number,1); 
sNorm = zeros(number,1); aNorm = zeros(number,1); 
aCond = zeros(number,1);

for i = 1:number
    %setting up orthogonal and diagonal matrices
    U = rand(50,50); V = rand(50,50); S = rand(50,1);
    S = sort(S,'descend');
    S = (1/S(1)).*S;
    S = diag(S.^6);
    [U,~] = qr(U);
    [V,~] = qr(V);
    A = U*S*V';
    
    [U2,S2,V2] = svd(A);
    
    sNorm(i) = norm(S-S2); aNorm(i) = norm(A-(U2*S2*V2'));
    
    %fixing the signs of the orthogonal matrices
    U2f = diag(U2'*U);
    U2f = U2*diag(U2f);
    
    V2f = diag(V2'*V);
    V2f = V2*diag(V2f);
    
    uNorm(i) = norm(U-U2f); vNorm(i) = norm(V-V2f);
    
    aCond(i) = cond(A);
end