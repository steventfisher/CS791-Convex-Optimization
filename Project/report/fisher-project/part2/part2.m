%% maximizing algebraic connectivity of a graph
max_alg_conn_data;
Q = null(ones(1,n)); % columns of Q are orthonormal basis for 1^\perp
cvx_begin
    variable w(m)
    L = A*diag(w)*A';
    maximize (lambda_min(Q'*L*Q))
    subject to
        w >= 0;
        F*w <= g;
cvx_end
w(abs(w) < 1e-4) = 0;

% compare algebraic connectivities
L_unif = (1/m)*A*A';
dunif = eig(L_unif);
dopt = eig(L);
fprintf(1, 'Algebraic connectivity of L_unif: %f\n', dunif(2));
fprintf(1, 'Algebraic connectivity of L_opt: %f\n', dopt(2));

% plot topology of constant-weight graph
figure(1), clf
gplot(L_unif,xy);
hold on;
plot(xy(:,1), xy(:,2), 'ko','LineWidth',4, 'MarkerSize',4);
axis([0.05 1.1 -0.1 0.95]);
title('Constant-weight graph')
hold off;
print -deps graph_plot1.eps;

% plot topology of optimal weight graph
figure(2), clf
gplot(L,xy);
hold on;
plot(xy(:,1), xy(:,2), 'ko','LineWidth',4, 'MarkerSize',4);
axis([0.05 1.1 -0.1 0.95]);
title('Optimal weight graph')
hold off;
print -deps graph_plot2.eps;