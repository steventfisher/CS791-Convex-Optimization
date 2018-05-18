% Gradient Descent FUnction with Exact Line Search
% Input
%   param_Q : nxn positive definite matrix
%   param_q : nx1 column vector
% Output
%   ret_x : optimal minimized x
%   ret_itr : number of iterations
%   ret_diff_f : difference in approximated objective
%                function value from actual value for
%                each iteration
% Optimation Problem
%    minimize f(x) = (1/2)x'Qx + q'x

function [ret_x, ret_itr, ret_diff_f] = gradientDescent2( param_Q, param_q, x )

% Correct solution for x
x_act = -param_Q \ param_q;
f_act = (1/2) * x_act' * param_Q * x_act + param_q' * x_act;

% Stopping Condition Limit
inta = 0.000001;

% Exact Line Search Limit
epsilon = 0.000001;

%Size of the input parameters
size_Q = size( param_Q);
disp( 'Size of matrix Q')
disp(size_Q)


size_q = size( param_q);
disp( 'Size of vector q')
disp(size_q)

% Make sure param_Q is a positive definite matrix
eigen_Q = eig( param_Q );
if ~( all(eigen_Q) > 0)
    disp('Parameter Q is not a positive definite matrix')
end

% Condition number of positive definite matrix param_Q
eigen_max_Q = max( eigen_Q );
eigen_min_Q = min( eigen_Q );
cond_num_ub = eigen_max_Q / eigen_min_Q;

disp('Condition number of the positive definite matrix Q is ')
disp(cond_num_ub);

% Make sure param_Q and param_q have correct dimensions
if ( size_Q(1) ~= size_Q(2) || size_Q(1) ~= size_q(1) )
    disp('Parameter Q and q have incorrect dimensions')
end

% Pick a random starting poitn x in dom(f)
%x = rand( size_q(1), 1);

% Gradient of objective function
grad_f = param_Q * x + param_q;

itr = 0;

while norm( grad_f ) > inta
    
    % Gradient of objective function
    grad_f = param_Q * x + param_q;
    
    % Determine the search direction
    search_dir = -grad_f;
    
    % Choose step-size using exact line search
    t_low = 0;
    t_up = 1;
    
    track_f = param_Q * ( x + t_up * search_dir ) + param_q;
        
    while( ( search_dir' * track_f ) < 0)
        t_low = t_up;
        t_up = t_up * 2;
    end
    
    while( ( t_up - t_low) > epsilon )
        t_mid = (1/2) * (t_low + t_up);
        track_f = param_Q * (x + t_up * search_dir) + param_q;
        if((search_dir' * track_f) < 0)
            t_low = t_mid;
        else
            t_up = t_mid;
        end
        t = (1/2) * (t_low + t_up);
    end
        
    %Update x with new value
    x = x + t * search_dir;
    
    %Update iteration count
    itr = itr + 1;
    
    % Difference in the objective funtion value fro the ideal value
    f_itr = (1/2) * x' * param_Q * x + param_q' * x;
    diff_f( itr ) = f_itr - f_act;
end

%Return from function
ret_x = x;
ret_itr = itr;
ret_diff_f = diff_f;