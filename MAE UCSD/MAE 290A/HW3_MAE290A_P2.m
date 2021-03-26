A = rand(500,10000); time = zeros(3,1);

tic
[L,U] = lu(A);
toc
time(1) = toc;

tic
[L,U] = LU(A);
toc
time(2) = toc;

tic
[L,U] = LUm(A); %modified algorithm
toc
time(3) = toc;