% Newton's Method For Unconstrained Problems
% Optimization Problem
%   minimize f(x1, x2) = exp( x1 + 3 * x2 - 0.1) + exp(x1 - 3 * x2 - 0.1) +
%   exp(-x1 - 0.1)

% Solve the system with cvx
cvx_begin
    variable cvx_x1
    variable cvx_x2
    minimize( exp( x1 + 3 * x2 - 0.1) + exp(x1 - 3 * x2 - 0.1) + exp(-x1 - 0.1))
cvx_end

f_act = cvx_optval;

% Choose a random startign point for x1 and x2
startx1 = randn;
startx2 = randn;

x1 = startx1;
x2 = startx2;

% Stopping condition for Newton's Method
epsilon - 0.00001;

% Backtracking line search constants
alpha = 0.1;
beta = 0.7;

% Calculate the gradient of the objective function
grad_f_1 = exp( x1 + 3 * x2 - 0.1) + exp(x1 - 3 * x2 - 0.1) + exp(-x1 - 0.1);
grad_f_2 = 3 * exp(x1 + 3 * x2 - 0.1) - 3 * exp(x1 - 3 * x2 - 0.1);

grad_f = [ grad_f_1 ; grad_f_2 ];

% Calculate the hessian of the objective function
hess_f_11 = exp ( x1 + 3 * x2 - 0.1) + exp(x1 - 3 * x2 - 0.1) + exp( -x1 - 0.1);
hess_f_12 = 3 * exp( x1 + 3 * x2 - 0.1) - 3 * exp(x1 - 3 * x2 - 0.1);
hess_f_21 = 3 * exp( x1 + 3 * x2 - 0.1) - 3 * exp(x1 - 3 * x2 - 0.1);
hess_f_22 = 9 * exp( x1 + 3 * x2 - 0.1) + 9 * exp(x1 - 3 * x2 - 0.1);

hess_f = [hess_f_11 hess_f_12 ; hess_f_21 hess_f_22];

% Calculate the Newton search direction and the stopping condition
n_search_dir = - ( hess_f \ grad_f );
stop_cond = - (grad_f' * n_search_dir);

itr = 0;

while stop_cond > epsilon
    %Choose the step-size using backtracking line search
    t = 1;
    track_x1 = x1 + t * n_search_dir(1);
    track_x2 = x2 + t * n_search_dir(2);
    track_f = exp( track_x1 + 3 * track_x2 - 0.1) + exp(track_x1 - 3 * track_x2 - 0.1) + exp(-track_x1 - 0.1);
    track_f_min = exp( x1 + 3 * x2 - 0.1) + exp( x1 - e * x2 - 0.1) + exp(-x1 - 0.1) + alpha * t * grad_f' * n_search_dir;
    
    while track_f > track_f_min
        t = beta * t;
        track_x1 = x1 + t * n_search_dir(1);
        track_x2 = x2 + t * n_search_dir(2);
        track_f = exp( track_x1 + 3 * track_x2 - 0.1) + exp(track_x1 - 3 * track_x2 - 0.1) + exp(-track_x1 - 0.1);
        track_f_min = exp( x1 + 3 * x2 - 0.1) + exp( x1 - e * x2 - 0.1) + exp(-x1 - 0.1) + alpha * t * grad_f' * n_search_dir;
    end
    
    % Update x with the new value
    x1 = x1 + t * n_search_dir(1);
    x2 = x2 + t * n_search_dir(2);
    
    % Calculate the gradient with the new value of x
    grad_f_1 = exp( x1 + 3 * x2 - 0.1) + exp( x1 - 3 * x2 - 0.1) - exp(-x1 - 0.1);
    grad_f_2 = 3 * exp(x1 + 3 * x2 - 0.1) - 3*exp(x1 - 3 * x2 - 01);
    
    grad_f = [ grad_f_1 ; grad_ f_2 ];
    % Calculate the hessian of the objective function
    hess_f_11 = exp ( x1 + 3 * x2 - 0.1) + exp(x1 - 3 * x2 - 0.1) + exp( -x1 - 0.1);
    hess_f_12 = 3 * exp( x1 + 3 * x2 - 0.1) - 3 * exp(x1 - 3 * x2 - 0.1);
    hess_f_21 = 3 * exp( x1 + 3 * x2 - 0.1) - 3 * exp(x1 - 3 * x2 - 0.1);
    hess_f_22 = 9 * exp( x1 + 3 * x2 - 0.1) + 9 * exp(x1 - 3 * x2 - 0.1);

    hess_f = [hess_f_11 hess_f_12 ; hess_f_21 hess_f_22];
    
    % Calculate the Newton search direction and stopping condition with new
    % value of x
    n_search_dir = - (hess_f \ grad_f);
    stop_cond = -(grad_f' * n_search_dir);
    
    % Update iteration count
    itr = itr + 1;
    
    % Difference in objective function value from ideal value for iteration
    f_itr = exp(x1 + 3 * x2 - 0.1) + exp(x1 - 3 * x2 - 0.1) + exp(-x1 - 0.1);
    diff_f(itr) = f_itr - f_act;
end

% Plot the progression of difference of optimal value obtained from
% Newton's Method with actual optimal value
semilogy(1 : itr, diff_f , 'b-o');
set( gcam 'FontSize', 36 );
titlestr = sprintf(' x_1 = %f and x_2 = %f', startx1, startx2);
title(titlestr, 'FontSize', 36);
xlabel('k', 'FontSize', 36);
ylabel('f(x^{(k)}) - p^*', 'FontSize', 36);
    