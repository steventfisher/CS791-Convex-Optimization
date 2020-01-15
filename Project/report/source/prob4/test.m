% Log-Barrier Method for Testing Feasiblility of LP (Phase I Optimization)

% Problem size variables
size_m = 100;
size_n = 50;

% Generate random coefficient matrix and RHS 
A = randn(size_m, size_n);
b = randn(size_m, 1);

% Random initial value for x
start_x = randn(size_n, 1);

% Initial value for s based on constraints
start_s = 2 * max(A * start_x - b);

% Barrier method variables
t = 1;
mu = 10;
epsilon = 0.000001;

% Counter for indexing
index = 0;

% Flag to end iterations
stop = 0;

% Counter for number of inequalities satisfied by x
numIneqSat = 0;

% Optimal values from previous iteration
x_1 = start_x;
s_1 = start_s;

while(~stop)
    %Centering Step
    [ x , s , itr ] = newtonLP( t, A, b, x_1, s_1 );
    
    index = index + 1;
    
    % Test the number of inequalities satisfied by x
    numIneqSat(index) = sum(A * x < b);
    
    % Tracking progression of x and s through each iteration
    x_t(:, index) = x;
    s_t(:, index) = s;
    itr_str(index) = itr;
    
    %Set the optimal values of the previous iteration as starting points
    %for the next iteration
    x_1 = x;
    s_1 + s;
    
    % Test the stopping condition
    if( x< 0 || abs(size_m/t) < epsilon)
        stop = 1;
    else
        t - mu * t;
    end
end