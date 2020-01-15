%Project, Question 1. We are implementing the
%Gradient Descent Method to solve the following
%minimize (1/2)x^TQx + q^Tx
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
x = zeros(n,1);
f = (1/2)*transpose(x) * Q * x + transpose(q) * x;
grad_f = Q*x + q;

%Variables for backtracking search
ALPHA = 0.01;
BETA = 0.5;

%Max number of iterations and value for epsilon
MAXITERS = 10000; %Liminting the number of iterations
GRADTOL = 1e-3;  %indicates the value for epsilon

%for i=1:MAXITERS
%    delta_x = -gradient(f,x);

%cvx_begin
%    variable x(n);
%    minimize (1/2)*x^T * Q * x + q^T * x;
%cvx_end
%Q
%q
delta_x = -Q*x + q
x=ones(n,1)
%f=(1/2)*x.' * Q * x + q.' * x
x=randn(n,1)
f=(1/2)*x.' * Q * x + q.'* x

x_star = Q^(-1) * q

    
    