%Runs numerical integration arbitrary RK method using variable time steps
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%p: how error scales with step size (error = k*hˆp)
%error_desired: the desired local truncation error at each step
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, num_evals, fail_rate, h_rec] = explicit_RK_variable_step_integration(rate_func_in,tspan,X0,h_ref,BT_struct, p, error_desired)
    ti = tspan(1); tf = tspan(2);
    N  = ceil((tf - ti)/h_ref); h = (tf - ti)/N;
    nx = numel(X0); 
    t_list = 1; t_list(1) = ti;
    X_list = zeros(nx, 1); X_list(:,1) = X0;
    num_evals = 0; XA = X0;
    h_rec = []; 
    num_fails = 0;
    tnow = ti;
    while tnow <= tf
        redo = true; 
        while redo == true
            h_prev = h;
            h = min(h,tf-tnow+1e-15);
            [XB, add_evals, h, redo] = explicit_RK_variable_step(rate_func_in, tnow, XA, h, BT_struct, p, error_desired);
            if redo == true
                num_fails = num_fails+1;
            end
            num_evals = num_evals + add_evals;
        end
        h_rec(end+1,:) = h_prev;
        tnow = tnow+h_prev;
        t_list(end+1) = tnow;
        X_list(:,end+1) = XB; XA = XB;
    end
    h_avg = mean(h_rec);
    fail_rate = num_fails/num_evals;
end