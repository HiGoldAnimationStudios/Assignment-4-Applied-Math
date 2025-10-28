function day_1()
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    
    V0 = [x0;y0;dxdt0;dydt0];
    t_range = linspace(0,10,100);
    V_list = compute_planetary_motion(t_range,V0,orbit_params);
    
    my_rate_func = @(t_in, V_in) gravity_rate_func(t_in,V_in,orbit_params);

    %axis equal; axis square;
    %axis([-20,20,-20,20])
    %hold on
    %plot(0,0,'ro','markerfacecolor','r','markersize',5);
    %plot(V_list(:,1),V_list(:,2),'k');

    tspan=[0,10];
    h_ref=0.02;
    
    %Forward Euler
    Forward_Euler=struct();
    Forward_Euler.A=[0];
    Forward_Euler.B=[1];
    Forward_Euler.C=[0];

    %Explict Midpoint
    Explict_Midpoint=struct();
    Explict_Midpoint.A=[0, 0; 0.5, 0];
    Explict_Midpoint.B=[0, 0.5];
    Explict_Midpoint.C=[0, 1];

    %Heun's Third Order Method
    Heun_Third=struct();
    Heun_Third.A=[0,0,0;1/3,0,0;0,2/3,0];
    Heun_Third.B=[1/4,0,3/4];
    Heun_Third.C=[0,1/3,2/3];

    %DormandPrince (because orion wanted to test it)
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

    BT_struct = Forward_Euler;
    
    [t_list_FE,X_list_FE,h_avg_FE, num_evals_FE] = explicit_RK_fixed_step_integration(my_rate_func,tspan,V0,h_ref,BT_struct);

    figure(); 
    subplot(2,1,1); hold on;
    plot(t_range,V_list(:,1),"r")
    plot(t_range,V_list(:,2),"b")
    plot(t_list_FE,X_list_FE(:,1),"r--")
    plot(t_list_FE,X_list_FE(:,2),"b--")
    xlabel("time"); ylabel("pos"); title("Forward Euler")

    subplot(2,1,2);hold on;
    plot(t_range,V_list(:,3),"r")
    plot(t_range,V_list(:,4),"b")
    plot(t_list_FE,X_list_FE(:,3),"r--")
    plot(t_list_FE,X_list_FE(:,4),"b--")
    xlabel("time"); ylabel("velocity")

    

    [t_list_EM,X_list_EM,h_avg_EM, num_evals_EM] = explicit_RK_fixed_step_integration(my_rate_func,tspan,V0,h_ref,BT_struct);

    figure(); 
    subplot(2,1,1); hold on;
    plot(t_range,V_list(:,1),"r")
    plot(t_range,V_list(:,2),"b")
    plot(t_list_EM,X_list_EM(:,1),"r--")
    plot(t_list_EM,X_list_EM(:,2),"b--")
    xlabel("time"); ylabel("pos"); title("Explicit Midpoint")

    subplot(2,1,2);hold on;
    plot(t_range,V_list(:,3),"r")
    plot(t_range,V_list(:,4),"b")
    plot(t_list_EM,X_list_EM(:,3),"r--")
    plot(t_list_EM,X_list_EM(:,4),"b--")
    xlabel("time"); ylabel("velocity")

    %Local Truncation Error
    n_samples=30;
    h_ref_list=logspace(-6,1,n_samples);
    tr_error_list_FE=zeros(1,n_samples);
    tr_error_list_EM=zeros(1,n_samples);

    for n=1:length(h_ref_list)
        h_ref=h_ref_list(n);
        V_list=compute_planetary_motion(tspan(1)+h_ref,V0,orbit_params);
        %Forward Euler
        BT_struct = Forward_Euler;
        [XB,~]=explicit_RK_step(my_rate_func,tspan(1),V0,h_ref,BT_struct);
        tr_error_list_FE(n)=norm(XB-V_list);
        %Explicit Midpoint
        BT_struct = Explict_Midpoint;
        [XB,~]=explicit_RK_step(my_rate_func,tspan(1),V0,h_ref,BT_struct);
        tr_error_list_EM(n)=norm(XB-V_list);
    end

    filter_params=struct();
    filter_params.min_y_val=1e-13;
    filter_params.max_y_val=1e-6;

    [p_FE,k_FE]=loglog_fit(h_ref_list,tr_error_list_FE,filter_params)
    [p_EM,k_EM]=loglog_fit(h_ref_list,tr_error_list_EM,filter_params)

    hh = logspace(log10(min(h_ref_list)), log10(max(h_ref_list)), 200);


    figure();
    loglog(h_ref_list,tr_error_list_FE,'o', 'MarkerSize',5, 'Color',[1 0 0], 'DisplayName','Forward Euler Error');
    hold on;
    loglog(hh, k_FE .*hh.^p_FE,'r', 'LineWidth',1, 'DisplayName', sprintf('Forward Euler fit: p=%.2f', p_FE));
    loglog(h_ref_list,tr_error_list_EM,'o', 'MarkerSize',5, 'Color',[0 1 0], 'DisplayName','Explicit Midpoint Error');
    loglog(hh, k_EM .*hh.^p_EM,'g', 'LineWidth',1, 'DisplayName', sprintf('Explicit Midpoint fit: p=%.2f', p_EM));
    title("Local Truncation Error"); xlabel("step size h"); ylabel("local error")
    legend("Location","best"); hold off
    
    %Global Truncation Error
    n_samples=20;
    h_ref_list=logspace(-4,-1,n_samples);

    h_avg_list_FE=zeros(1,n_samples);
    num_evals_list_FE=zeros(1,n_samples);
    tr_error_list_FE=zeros(1,n_samples);

    h_avg_list_EM=zeros(1,n_samples);
    num_evals_list_EM=zeros(1,n_samples);
    tr_error_list_EM=zeros(1,n_samples);

    
    for n=1:length(h_ref_list)
        h_ref=h_ref_list(n);

        BT_struct = Forward_Euler;

        [t_list_FE,X_list_FE,h_avg_FE, num_evals_FE] = explicit_RK_fixed_step_integration(my_rate_func,tspan,V0,h_ref,BT_struct);
        V_list=compute_planetary_motion(t_list_FE-tspan(1),V0,orbit_params);
        tr_error=norm(X_list_FE(end,:)-V_list(end,:));
        tr_error_list_FE(n)=tr_error;
        num_evals_list_FE(n)=num_evals_FE;
        h_avg_list_FE(n)=h_avg_FE;

        BT_struct = Explict_Midpoint;
        [t_list_EM,X_list_EM,h_avg_EM, num_evals_EM] = explicit_RK_fixed_step_integration(my_rate_func,tspan,V0,h_ref,BT_struct);
        V_list=compute_planetary_motion(t_list_FE-tspan(1),V0,orbit_params);
        tr_error=norm(X_list_EM(end,:)-V_list(end,:));
        tr_error_list_EM(n)=tr_error;
        num_evals_list_EM(n)=num_evals_EM;
        h_avg_list_EM(n)=h_avg_EM;
    end

    filter_params=struct();

    % filter_params.min_yval=1e-10;
    filter_params.max_yval=10;
    filter_params.max_xval=10^(-2.5);


    [p_FE_h,k_FE_h]=loglog_fit(h_avg_list_FE,tr_error_list_FE,filter_params)   
    [p_EM_h,k_EM_h]=loglog_fit(h_avg_list_EM,tr_error_list_EM,filter_params)   

    filter_params=struct();
    filter_params.min_xval=1e4;

    
    [p_EM_ne,k_EM_ne]=loglog_fit(num_evals_list_EM,tr_error_list_EM,filter_params) 
    [p_FE_ne,k_FE_ne]=loglog_fit(num_evals_list_FE,tr_error_list_FE,filter_params)  

    figure();
    loglog(h_avg_list_FE,tr_error_list_FE,'o', 'MarkerSize',5, 'Color','r', 'DisplayName','Forward Euler Error')
    hold on;
    loglog(h_avg_list_FE,k_FE_h*h_avg_list_FE.^p_FE_h,"b")

    title("Global Truncation Error vs Number of Function Calls"); xlabel("step size h"); ylabel("global error")
    legend("Location","best"); hold off
    
    figure();
    loglog(num_evals_list_FE,tr_error_list_FE,"r")
    hold on;
    loglog(num_evals_list_FE,k_FE_ne*num_evals_list_FE.^p_FE_ne,"b")

    title("Global Truncation Error vs Number of Function Calls"); xlabel("number of function calls"); ylabel("global error")
    legend("Location","best"); hold off
end