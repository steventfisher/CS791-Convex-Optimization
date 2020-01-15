cx_begin
    variable x(1)
    minimize( x^2 + 1 )
    subject to
        (x - 2)*(x - 4) <= 0
cvx_end