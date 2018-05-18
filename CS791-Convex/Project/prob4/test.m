% define the Hesssian of the objective
function h = Hessian(x1,x2)
h = [exp(x1+3.*x2-0.1)+exp(x1-3.*x2-0.1)+exp(-x1-0.1) 3.*exp(x1+3.*x2-0.1)-3.*exp(x1-3.*x2-0.1)
    3.*exp(x1+3.*x2-0.1)-3.*exp(x1-3.*x2-0.1) 9.*exp(x1+3.*x2-0.1)+9.*exp(x1-3.*x2-0.1)];
end

% define the gradient of the objective
function g = grad(x1,x2)
g = [exp(x1+3.*x2-0.1)+exp(x1-3.*x2-0.1)-exp(-x1-0.1);3.*exp(x1+3.*x2-0.1)-3.*exp(x1-3.*x2-0.1)];
end

function t = backtrack_linesearch(f,dx,x,beta,alpha)
t = 1;
fk = feval(f,x);
xx = x;
x = x + t*dx;
fk1 = feval(f,x);
while fk1 > fk + alpha*t*(dx'*dx),
  ff=fk + alpha*t*(dx'*dx);
  t = t*beta;
  x = xx + t*dx;
  fk1 = feval(f,x);
end

clear all
clc
 
% define starting point
% x0 = round(3*rand(2,1));
x0=[-1;2];
% termination tolerance
tol = 1e-6;
 
% maximum number of allowed iterations
maxiter = 1000;
 
% minimum allowed perturbation
dxmin = 1e-6;
 
% step size ( beta=0.8 causes instability, beta=0.6 accurate with lowest iteration number)
alpha = 0.01;
beta=0.5;
% initialize gradient norm, optimization vector, iteration counter, perturbation
gnorm = inf;
x = x0;
niter = 0;
dx = inf;
 
% define the objective function:
f=@(x1,x2) exp(x1+3*x2-0.1)+exp(x1-3*x2-0.1)+exp(-x1-0.1);
% plot objective function contours for visualization for 2D:
figure(1); clf; fcontour(f); %axis equal; 
hold on
 
% redefine objective function syntax for use with optimization:
f2 = @(x) f(x(1),x(2));
 
% gradient descent algorithm:
while and(gnorm>=tol, and(niter <= maxiter, dx >= dxmin))
    % calculate gradient:
    g = grad(x(1),x(2));
    gnorm = norm(g);
    d=-g;
    % take step:
    t=backtrack_linesearch(f2,d,x,beta,alpha);
    xnew = x + t*d;
    % check step
    if ~isfinite(xnew)
        display(['Number of iterations: ' num2str(niter)])
        error('x is inf or NaN')
    end
    
    % plot current point for 2D
    plot([x(1) xnew(1)],[x(2) xnew(2)],'ko-')
    refresh
 
    % update termination metrics
    niter = niter + 1;
    dx = norm(xnew-x);
    x = xnew;
    
end
xopt = x;
fopt = f2(xopt);
niter = niter - 1;
end
...............................................

clear all
clc
 
% define starting point
 x0 = randn(2,1);
%x0=[-2;-1];
% termination tolerance
tol = 1e-6;
 
% maximum number of allowed iterations
maxiter = 1000;
 
% minimum allowed perturbation
dxmin = 1e-6;
 
% step size 
alpha = 0.1;
beta=0.7;
% initialize gradient norm, optimization vector, iteration counter, perturbation
landa = inf;
x = x0;
niter = 0;
dx = inf;
 
% define the objective function:
f=@(x1,x2) exp(x1+3*x2-0.1)+exp(x1-3*x2-0.1)+exp(-x1-0.1);
 
% redefine objective function syntax for use with optimization:
f2 = @(x) f(x(1),x(2));
 
% plot objective function contours for visualization for 2D:
% figure(1); clf; fcontour(f,[-5 5;-5 5]); axis equal; hold on
 
[X,Y]=meshgrid(-5:0.1:5,-5:0.1:5);
Z=exp(X+3.*Y-0.1)+exp(X-3.*Y-0.1)+exp(-X-0.1);
figure(1); clf; contour(X,Y,Z); axis equal; hold on
 
% gradient descent algorithm:
while and((landa/2)>=tol,(niter <= maxiter))%, dx >= dxmin
    
    % calculate gradient:
    g = grad(x(1),x(2));
    % calculate hessian 
    h=Hessian(x(1),x(2));
    landa = g'*inv(h)*g;
    % direction value
    d=-inv(h)*g;
    % take step:
    t=backtrack_linesearch(f2,d,x,beta,alpha);
    xnew = x + t*d;
 
    % check step
    if ~isfinite(xnew)
        display(['Number of iterations: ' num2str(niter)])
        error('x is inf or NaN')
    end
    
    % plot current point for 2D
    plot([x(1) xnew(1)],[x(2) xnew(2)],'ko-')
    refresh
 
    % update termination metrics
    niter = niter + 1;
 
    % calculate error
    xx=[-0.3465;1.612e-07];
    p_star=feval(f2,xx);
    fiter=f2(x);  
    error(niter)=fiter-p_star;
 
    % dx = norm(xnew-x);
    x = xnew;
    
end
xopt = x;
fopt = f2(xopt);
niter = niter - 1;
figure(2)
%%
loglog(error,'ko-')
ylabel('f(x)-p')
xlabel('k')
title('x1=-0.0117,x2=-.79')
