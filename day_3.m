function day_3()
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1/330000;
    orbit_params.G = 4*pi^2/orbit_params.m_sun;
    ti = 0;
    % tf = 10;
    tf = 100;
    % V0 = [1; 0; 0; 6.28;];
    V0 = [1.8; 0; 0; 6.28;];

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

    tspan = [ti,tf];
    my_rate_func = @(t, V) gravity_rate_func(t, V, orbit_params);

    h0 = 0.1;
    %h0 = 1;
    %test working
    des_err = 1e-5;
    [t_list,X_list,h_avg, num_evals, num_fails, h_rec] = variable_step_integration_with_fails(my_rate_func, tspan, V0, h0, DormandPrince, 5, des_err);

    figure
    hold on
    plot(X_list(1,:), X_list(2,:), "r", "DisplayName", "Dormand-Prince Aprroximation of Planet Path",'LineWidth',1)
    plot(0,0, 'ko','MarkerFaceColor','y','MarkerSize',10, "DisplayName","Sun Position")
    axis equal
    xlabel("x position")
    ylabel("y position")
    legend("Location","best")
    title('Position of the Planet')

    h_ref_list = logspace(-5,-1, 30);
    num_evals_list = [];
    h_avg_list = [];
    tr_error_list = [];


    fail_rate_list=[];

    %for your chosen method, do the following... (variable runge kutta)
    n_samples = 60;
    h_ref_list = logspace(-3, 1, n_samples);

    abs_diff_list = zeros(1, n_samples);
    approx_diff_list = zeros(1, n_samples);
    V_real = gravity_rate_func(ti,V0,orbit_params);
    XB1s = zeros(4, n_samples);
    XB2s = zeros(4, n_samples);
    for i = 1:length(h_ref_list)
        h_ref = h_ref_list(i);
        V_list = compute_planetary_motion(ti+h_ref,V0,orbit_params);
        % V_next = gravity_rate_func(t,V,orbit_params)
        [XB1, XB2, ~] = RK_step_embedded(my_rate_func, ti, V0, h_ref, DormandPrince);
        % XB1s(:, i) = XB1';
        % XB2s(:, i) = XB2';
        abs_diff_list(i) = norm(V_list-V0);
        approx_diff_list(i) = norm(XB1-XB2);
        tr_error_list1(i) = norm(XB1-V_list);
        tr_error_list2(i) = norm(XB2-V_list);

        [t_list,X_list,h_avg, num_evals, num_fails, h_rec] = variable_step_integration_with_fails(my_rate_func, tspan, V0, h_ref, DormandPrince, 5, des_err);
        fail_rate_list(end+1)=num_fails/num_evals;
        h_avg_list(end+1)=h_avg;
    end

    figure
    %fix - global error
    loglog(h_ref_list, abs_diff_list, DisplayName='f(t+h)-f(t)')
    hold on
    loglog(h_ref_list, approx_diff_list, DisplayName='XB2-XB1')
    loglog(h_ref_list, tr_error_list1, DisplayName=['XB1 ' char(949) '_{local}'])
    loglog(h_ref_list, tr_error_list2, DisplayName=['XB2 ' char(949) '_{local}'])
    legend(Location='southeast')

    
    % Plot the local truncation errors of XB1 and XB2 as a function of their difference, |XB1 − XB2|
    figure
    loglog(approx_diff_list, tr_error_list1, DisplayName=['XB1 ' char(949) '_{local}'])
    hold on
    loglog(approx_diff_list, tr_error_list2, DisplayName=['XB2 ' char(949) '_{local}'])
    loglog(approx_diff_list, approx_diff_list, 'k--', DisplayName='(XB1-XB2)')
    title(['XB1 and XB2 ', char(949), '_{local}' ' vs (XB1-XB2)'])


    figure()
    plot(h_avg_list, fail_rate_list)
    xlabel("h_a_v_g")
    ylabel("fail_rate")
    title("Failure rate vs Average step size")
    legend("Location","best")
end