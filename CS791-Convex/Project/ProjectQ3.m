%Problem 10.15 cvx book

%Project, Question 3. We are implementing the
%Newton's Method with Equality Constraints
%to solve the following
%minimize f(x) = SUM x_i log x_i
%subject to Ax = b
%
%
%Algorithm
%given a starting point x in the dom f
%repeat
%   1. delta_x = -gradient(f(x))
%   2. Line search. Choose step size t via exact or
%      backtracking line search
%   3. Update x = x + t delta_x
%stop is norm(gradient(f(x))) <= epsilon
%
%
%For exact line search t = argmin_t f(x + t delta_x)
%For backtracking line search
%   starting at t = 1, repeat t = beta t until
%   f(x + t delta_x) < f(x) + alpha t Hessian f(x)^T delta_x


%Variables for problem
n = 10;
randn('state',1);
Q = randn(n,n);
q = randn(n,1);
x = zero(n,1);
f = (1/2)*transpose(x) * Q * x + transpose(q) * x;

%Variables for backtracking search
ALPHA = 0.01;
BETA = 0.5;

%Max number of iterations and value for epsilon
MAXITERS = 10000;
GRADTOL = 1e-3; 

%for i=1:MAXITERS
%    delta_x = -gradient(f,x);


x_star = Q^(-1) * q

    
    