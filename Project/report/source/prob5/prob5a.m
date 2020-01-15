% Compressive Sensing
% Optimization Problem
%    minimize ||x||_1
%        s.t. y = phi * x
% where phi = [k,n] matrix with k << n
%         x = [n,1] vector with S < n non-zero elements

% Size of x
size_n = 64;

% Number of non-zero elements in x
nz_S = 8;

% Vector x as a comb shape of size nz_S
x = zeros(size_n, 1);
for index1 = 1 : nz_S;
    x(index1 * (size_n / nz_S)) = 1;
end

% Iteration counter
itr = 0;

% Set of selectable frequencies
freq = (0 : size_n - 1);

%Pick k using constraints
size_k = size_n - nz_S;

% Non uniformly distributed frequency choices
freq_k = freq;
freq_k(~mod(freq_k, (size_n/nz_S))) = [];

% Form the DFT matrix
phi = zeros( size_k, size_n);
for index1 = 1 : size_k
    for index2 = 1 : size_n
        phi(index1, index2) = size_n ^(-0.5) * exp(-1i * 2 * pi * freq_k(index1) * ((index2 - 1)/size_n));
    end
end

% Form the Sample Matrix
y = phi * x;

% Solve the system using CVX
cvx_begin
    variable cvx_x(size_n)
    minimize(norm(cvx_x, 1))
    subject to
        phi * cvx_x == y
cvx_end

