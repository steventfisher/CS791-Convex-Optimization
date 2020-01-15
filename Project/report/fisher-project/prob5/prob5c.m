% Compressive Sensing
% Optimization Problem
%    minimize ||x||_1
%        s.t. y = phi * x
% where phi = [k,n] matrix with k << n
%         x = [n,1] vector with S < n non-zero elements
%         e = upper bound of norm of noise vector

% Size of x
size_n = 100;

% Number of non-zero elements in x
nz_S = 30;

% Noise levels
e = [0.0001 0.0005 0.001 0.005 0.01 0.05];

% Vector x with nz_S number of uniformly distributed values
x = sprand(size_n, 1, (nz_S./size_n));
opt_f = norm(x,1);

% Set of selectable frequencies
freq = (0 : size_n - 1);

% Pick k using constraints
size_k = 58;

% Set of frequency randomly chosen from set of selectable frequencies
freq_k = sort(randsample(freq, size_k));

% Form the DFT matrix
phi = zeros( size_k, size_n);
for index1 = 1 : size_k
    for index2 = 1 : size_n
        phi(index1, index2) = size_n ^(-0.5) * exp(-1i * 2 * pi * freq_k(index1) * ((index2 - 1)/size_n));
    end
end

% Randomly generated Noise Vectors
n1 = rand(size_k, 1);
while(norm(n1) > e(1))
    n1 = n1.*rand(size_k, 1);
end

n2 = rand(size_k, 1);
while(norm(n2) > e(2))
    n2 = n2.*rand(size_k, 1);
end

n3 = rand(size_k, 1);
while(norm(n3) > e(3))
    n3 = n3.*rand(size_k, 1);
end

n4 = rand(size_k, 1);
while(norm(n4) > e(4))
    n4 = n4.*rand(size_k, 1);
end

n5 = rand(size_k, 1);
while(norm(n5) > e(5))
    n5 = n5.*rand(size_k, 1);
end

n6 = rand(size_k, 1);
while(norm(n6) > e(6))
    n6 = n6.*rand(size_k, 1);
end

% Form the Sample Matrix
y1 = phi * x + n1;
y2 = phi * x + n2;
y3 = phi * x + n3;
y4 = phi * x + n4;
y5 = phi * x + n5;
y6 = phi * x + n6;

% Solve the system using CVX
cvx_begin
    variable cvx_x1(size_n)
    minimize(norm(cvx_x1, 1))
    subject to
        norm(y1 - phi * cvx_x1) <= e(1);
cvx_end

cvx_begin
    variable cvx_x2(size_n)
    minimize(norm(cvx_x2, 1))
    subject to
        norm(y2 - phi * cvx_x2) <= e(2);
cvx_end

cvx_begin
    variable cvx_x3(size_n)
    minimize(norm(cvx_x3, 1))
    subject to
        norm(y3 - phi * cvx_x3) <= e(3);
cvx_end

cvx_begin
    variable cvx_x4(size_n)
    minimize(norm(cvx_x4, 1))
    subject to
        norm(y4 - phi * cvx_x4) <= e(4);
cvx_end

cvx_begin
    variable cvx_x5(size_n)
    minimize(norm(cvx_x5, 1))
    subject to
        norm(y5 - phi * cvx_x5) <= e(5);
cvx_end

cvx_begin
    variable cvx_x6(size_n)
    minimize(norm(cvx_x6, 1))
    subject to
        norm(y6 - phi * cvx_x6) <= e(6);
cvx_end

% Calculate the norm of differences between noiseless recovery
% and noisy recovery
nrm(1) = norm(x - cvx_x1);
nrm(2) = norm(x - cvx_x2);
nrm(3) = norm(x - cvx_x3);
nrm(4) = norm(x - cvx_x4);
nrm(5) = norm(x - cvx_x5);
nrm(6) = norm(x - cvx_x6);

% Plot the norm vs. error bound
semilogx(e, nrm, 'b-o');
set(gca, 'FontSize', 12);
xlabel('Upper bound on norm of noise vector \epsilon', 'FontSize', 12);
ylabel('||x^* - x||_x', 'FontSize', 12);