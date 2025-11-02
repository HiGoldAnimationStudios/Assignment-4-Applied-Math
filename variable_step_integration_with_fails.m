function [t_list, X_list, h_avg, num_fails, num_evals, h_rec] = variable_step_integration_with_fails(rate_func_in, step_func, tspan, X0, h_ref, BT_struct, p, error_desired)
    ti = tspan(1); tf = tspan(2);
    N  = ceil((tf - ti)/h_ref); h = (tf - ti)/N;
    nx = numel(X0); 
    t_list = 1; t_list(1) = ti;
    X_list = zeros(nx, 1); X_list(:,1) = X0;
    num_evals = 0; XA = X0;
    h_rec = []; 
    num_fails = 0;
    t = ti;
    while t <= tf
        redo = true; 
        while redo == true
            h_prev = h;
            h = min(h,tf-t+1e-15);
            [XB, add_evals, h, redo] = step_func(rate_func_in, t, XA, h, BT_struct, p, error_desired);
            if redo == true
                num_fails = num_fails+1;
            end
            num_evals = num_evals + add_evals;
        end
        h_rec(end+1,:) = h_prev;
        t = t+h_prev;
        t_list(end+1) = t;
        X_list(:,end+1) = XB; XA = XB;
    end
    h_avg = mean(h_rec);
end
