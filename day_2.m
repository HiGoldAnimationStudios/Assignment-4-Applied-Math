function day_2()

    DormandPrince = struct();
    DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;...
    5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    DormandPrince.A = [0,0,0,0,0,0,0;
    1/5, 0, 0, 0,0,0,0;...
    3/40, 9/40, 0, 0, 0, 0,0;...
    44/45, -56/15, 32/9, 0, 0, 0,0;...
    19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
    9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
    35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];
    
    t = 5;
    XA = 1;
    h = 0.5;
    my_rate = @(t, XA) rate_func01(t,XA);


    [XB1, XB2, num_evals] = RK_step_embedded(my_rate,t,XA,h,DormandPrince);

    % disp(['XB1 = ', num2str(XB1)])
    % disp(['XB2 = ', num2str(XB2)])
    % disp(['num_evals = ', num2str(num_evals)])

    X = solution01(t+h);
    %disp(['X(',num2str(t+h),') = ',num2str(X)])
    
    XB_diff = norm(XB1 - XB2);

    %disp(XB_diff)
 
    %plotting XB_diff vs h_ref_list at one specific time constant
    n_samples=160;
    h_ref_list=logspace(-6,1,n_samples);
    XB_diff_list = zeros(1,n_samples);

    for i = length(h_ref_list)
        
        h = h_ref_list(i); 
        [XB1, XB2, ~] = RK_step_embedded(my_rate,t,XA,h,DormandPrince);
        XB_diff = norm(XB1 - XB2);
        XB_diff_list(i) = XB_diff;
        
    end 
    
    disp(XB_diff_list)
    figure; 
    plot(h_ref_list, XB_diff_list)
end 