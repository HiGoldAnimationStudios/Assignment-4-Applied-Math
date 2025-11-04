function day_1()
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1/330000;
    orbit_params.G = 4*pi^2/orbit_params.m_sun;
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
    Explict_Midpoint.A=[0 0; 0.5 0];
    Explict_Midpoint.B=[0 1];
    Explict_Midpoint.C=[0 0.5];

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

    %{
    figure(); 
    subplot(2,1,1); hold on;
    plot(t_range,V_list(:,1),"r", "DisplayName","x exact position")
    plot(t_range,V_list(:,2),"b", "DisplayName","y exact position")
    plot(t_list_FE,X_list_FE(:,1),"r--", "DisplayName","x predicted position")
    plot(t_list_FE,X_list_FE(:,2),"b--", "DisplayName","y predicted position")
    xlabel("time"); ylabel("pos"); title("Forward Euler Method")
    legend("Location", "best")

    subplot(2,1,2);hold on;
    plot(t_range,V_list(:,3),"r", "DisplayName","x exact velocity")
    plot(t_range,V_list(:,4),"b", "DisplayName","y exact velocity")
    plot(t_list_FE,X_list_FE(:,3),"r--", "DisplayName","x predicted velocity")
    plot(t_list_FE,X_list_FE(:,4),"b--", "DisplayName","y predicted velocity")
    xlabel("time"); ylabel("velocity")
    legend("Location", "best")
    %}

    BT_struct=Explict_Midpoint;

    [t_list_EM,X_list_EM,h_avg_EM, num_evals_EM] = explicit_RK_fixed_step_integration(my_rate_func,tspan,V0,h_ref,BT_struct);
    
    %{
    figure(); 
    subplot(2,1,1); hold on;
    plot(t_range,V_list(:,1),"r", "DisplayName","x exact position")
    plot(t_range,V_list(:,2),"b", "DisplayName","y exact position")
    plot(t_list_EM,X_list_EM(:,1),"r--", "DisplayName","x predicted position")
    plot(t_list_EM,X_list_EM(:,2),"b--", "DisplayName","y predicted position")
    xlabel("time"); ylabel("pos"); title("Explicit Midpoint Method")
    legend("Location", "best")

    subplot(2,1,2);hold on;
    plot(t_range,V_list(:,3),"r", "DisplayName","x exact velocity")
    plot(t_range,V_list(:,4),"b", "DisplayName","y exact velocity")
    plot(t_list_EM,X_list_EM(:,3),"r--", "DisplayName","x predicted velocity")
    plot(t_list_EM,X_list_EM(:,4),"b--", "DisplayName","y predicted velocity")
    xlabel("time"); ylabel("velocity")
    legend("Location", "best")
    %}

    BT_struct=Heun_Third;
    
    [t_list_HM,X_list_HM,h_avg_HM, num_evals_HM] = explicit_RK_fixed_step_integration(my_rate_func,tspan,V0,h_ref,BT_struct);

    %{
    figure(); 
    subplot(2,1,1); hold on;
    plot(t_range,V_list(:,1),"r", "DisplayName","x exact position")
    plot(t_range,V_list(:,2),"b", "DisplayName","y exact position")
    plot(t_list_HM,X_list_HM(:,1),"r--", "DisplayName","x predicted position")
    plot(t_list_HM,X_list_HM(:,2),"b--", "DisplayName","y predicted position")
    xlabel("time"); ylabel("pos"); title("Heun's Third Order Method")
    legend()

    subplot(2,1,2);hold on;
    plot(t_range,V_list(:,3),"r", "DisplayName","x exact velocity")
    plot(t_range,V_list(:,4),"b", "DisplayName","y exact velocity")
    plot(t_list_HM,X_list_HM(:,3),"r--", "DisplayName","x predicted velocity")
    plot(t_list_HM,X_list_HM(:,4),"b--", "DisplayName","y predicted velocity")
    xlabel("time"); ylabel("velocity")
    legend()
    %}

    %Local Truncation Error
    n_samples=30;
    h_ref_list=logspace(-6,1,n_samples);
    tr_error_list_FE=zeros(1,n_samples);
    tr_error_list_EM=zeros(1,n_samples);
    tr_error_list_HM=zeros(1,n_samples);

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
        %Heun Third
        BT_struct = Heun_Third;
        [XB,~]=explicit_RK_step(my_rate_func,tspan(1),V0,h_ref,BT_struct);
        tr_error_list_HM(n)=norm(XB-V_list);
    end

    filter_params=struct();
    filter_params.min_y_val=1e-13;
    filter_params.max_y_val=1e-6;

    filter_params.min_xval=10^(-3);

    [p_FE,k_FE]=loglog_fit(h_ref_list,tr_error_list_FE,filter_params)
    [p_EM,k_EM]=loglog_fit(h_ref_list,tr_error_list_EM,filter_params)
    [p_HM,k_HM]=loglog_fit(h_ref_list,tr_error_list_HM,filter_params)

    hh = logspace(log10(min(h_ref_list)), log10(max(h_ref_list)), 200);


    figure();
    loglog(h_ref_list,tr_error_list_FE,'o', 'MarkerSize',5, 'Color',[1 0 0], 'DisplayName','Forward Euler Error');
    hold on;
    loglog(hh, k_FE .*hh.^p_FE,'r', 'LineWidth',1, 'DisplayName', sprintf('Forward Euler fit: p=%.2f', p_FE));
    loglog(h_ref_list,tr_error_list_EM,'o', 'MarkerSize',5, 'Color',[0 0 1], 'DisplayName','Explicit Midpoint Error');
    loglog(hh, k_EM .*hh.^p_EM,'b', 'LineWidth',1, 'DisplayName', sprintf('Explicit Midpoint fit: p=%.2f', p_EM));
    loglog(h_ref_list,tr_error_list_HM,'o', 'MarkerSize',5, 'Color',[0 1 0], 'DisplayName','Heun Third Error');
    loglog(hh, k_HM .*hh.^p_HM,'g', 'LineWidth',1, 'DisplayName', sprintf('Heun Third fit: p=%.2f', p_HM));
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

    h_avg_list_HM=zeros(1,n_samples);
    num_evals_list_HM=zeros(1,n_samples);
    tr_error_list_HM=zeros(1,n_samples);

    
    for n=1:length(h_ref_list)
        h_ref=h_ref_list(n);

        BT_struct = Forward_Euler;
        [t_list_FE,X_list_FE,h_avg_FE, num_evals_FE] = explicit_RK_fixed_step_integration(my_rate_func,tspan,V0,h_ref,BT_struct);
        V_list_FE=compute_planetary_motion(t_list_FE-tspan(1),V0,orbit_params);
        tr_error=norm(X_list_FE(end,:)-V_list_FE(end,:));
        tr_error_list_FE(n)=tr_error;
        num_evals_list_FE(n)=num_evals_FE;
        h_avg_list_FE(n)=h_avg_FE;

        BT_struct = Explict_Midpoint;
        [t_list_EM,X_list_EM,h_avg_EM, num_evals_EM] = explicit_RK_fixed_step_integration(my_rate_func,tspan,V0,h_ref,BT_struct);
        V_list_EM=compute_planetary_motion(t_list_EM-tspan(1),V0,orbit_params);
        tr_error=norm(X_list_EM(end,:)-V_list_EM(end,:));
        tr_error_list_EM(n)=tr_error;
        num_evals_list_EM(n)=num_evals_EM;
        h_avg_list_EM(n)=h_avg_EM;

        BT_struct = Heun_Third;
        [t_list_HM,X_list_HM,h_avg_HM, num_evals_HM] = explicit_RK_fixed_step_integration(my_rate_func,tspan,V0,h_ref,BT_struct);
        V_list_HM=compute_planetary_motion(t_list_HM-tspan(1),V0,orbit_params);
        tr_error=norm(X_list_HM(end,:)-V_list_HM(end,:));
        tr_error_list_HM(n)=tr_error;
        num_evals_list_HM(n)=num_evals_HM;
        h_avg_list_HM(n)=h_avg_HM;
    end

    filter_params=struct();

    % filter_params.min_yval=1e-10;
    filter_params.max_yval=10;
    filter_params.max_xval=10^(-2.5);


    [p_FE_h,k_FE_h]=loglog_fit(h_avg_list_FE,tr_error_list_FE,filter_params)   
    [p_EM_h,k_EM_h]=loglog_fit(h_avg_list_EM,tr_error_list_EM,filter_params)   
    [p_HM_h,k_HM_h]=loglog_fit(h_avg_list_HM,tr_error_list_HM,filter_params)   

    filter_params=struct();
    filter_params.min_xval=1e4;

    
    [p_EM_ne,k_EM_ne]=loglog_fit(num_evals_list_EM,tr_error_list_EM,filter_params) 
    [p_FE_ne,k_FE_ne]=loglog_fit(num_evals_list_FE,tr_error_list_FE,filter_params)  
    [p_HM_ne,k_HM_ne]=loglog_fit(num_evals_list_HM,tr_error_list_HM,filter_params)  

    figure();
    loglog(h_avg_list_FE,tr_error_list_FE,'o', 'MarkerSize',5, 'Color','r', 'DisplayName','Forward Euler Error')
    hold on;
    loglog(h_avg_list_FE,k_FE_h*h_avg_list_FE.^p_FE_h,"r", 'DisplayName','Forward Euler Fit Line')
    loglog(h_avg_list_EM,tr_error_list_EM,'o', 'MarkerSize',5, 'Color','b', 'DisplayName','Explicit Midpoint Error')
    loglog(h_avg_list_EM,k_EM_h*h_avg_list_EM.^p_EM_h,"b", 'DisplayName','Explicit Midpoint Fit Line')
    loglog(h_avg_list_HM,tr_error_list_HM,'o', 'MarkerSize',5, 'Color','g', 'DisplayName','Heun Third Error')
    loglog(h_avg_list_HM,k_HM_h*h_avg_list_HM.^p_HM_h,"g", 'DisplayName','Heun Third Fit Line')

    title("Global Truncation Error vs Step Size h"); xlabel("step size h"); ylabel("global error")
    legend("Location","best"); hold off
    
    figure();
    loglog(num_evals_list_FE,tr_error_list_FE,'o', 'MarkerSize',5, 'Color','r', 'DisplayName','Forward Euler Error')
    hold on;
    loglog(num_evals_list_FE,k_FE_ne*num_evals_list_FE.^p_FE_ne,"r", 'DisplayName','Forward Euler Fit Line')
    loglog(num_evals_list_EM,tr_error_list_EM, 'o', 'MarkerSize',5, 'Color','b', 'DisplayName','Explicit Midpoint Error')
    loglog(num_evals_list_EM,k_EM_ne*num_evals_list_EM.^p_EM_ne,"b", 'DisplayName','Explicit Midpoint Fit Line')
    loglog(num_evals_list_HM,tr_error_list_HM, 'o', 'MarkerSize',5, 'Color','g', 'DisplayName','Heun Third Error')
    loglog(num_evals_list_HM,k_HM_ne*num_evals_list_HM.^p_HM_ne,"g", 'DisplayName','Heun Third Fit Line')

    title("Global Truncation Error vs Number of Function Calls"); xlabel("number of function calls"); ylabel("global error")
    legend("Location","best"); hold off

    %Conservation of Physical Quantities
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1/330000;
    orbit_params.G = 4*pi^2/orbit_params.m_sun;
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    
    V0 = [x0;y0;dxdt0;dydt0];
    tspan=[0,10];
    h=0.1;
    my_rate_func = @(t_in, V_in) gravity_rate_func(t_in,V_in,orbit_params);
    
    [t_list_FE, X_list_FE] = explicit_RK_fixed_step_integration(my_rate_func, tspan, V0, h, Forward_Euler);
    [E_FE,H_FE] = physical_quantities(X_list_FE, orbit_params); 

    [t_list_EM, X_list_EM] = explicit_RK_fixed_step_integration(my_rate_func, tspan, V0, h, Explict_Midpoint);
    [E_EM,H_EM] = physical_quantities(X_list_EM, orbit_params); 

    [t_list_HM, X_list_HM] = explicit_RK_fixed_step_integration(my_rate_func, tspan, V0, h, Heun_Third);
    [E_HM,H_HM] = physical_quantities(X_list_HM, orbit_params); 

    figure()
    plot(t_list_FE, E_FE, "r", 'DisplayName',"Forward Euler")
    hold on
    plot(t_list_EM, E_EM, "b--", 'DisplayName',"Explicit Midpoint")
    plot(t_list_HM, E_HM, "go", 'DisplayName',"Forward Euler")
    title("Mechanical Energy vs time")
    xlabel("time (t)"); ylabel("mechanical energy (E)");
    legend();

    figure()
    plot(t_list_FE, H_FE, "r", 'DisplayName',"Forward Euler")
    hold on
    plot(t_list_EM, H_EM, "b--", 'DisplayName',"Explicit Midpoint")
    plot(t_list_HM, H_HM, "go", 'DisplayName',"Forward Euler")
    title("Angular Momentum vs time")
    xlabel("time (t)"); ylabel("angular momentum (H)");
    legend();
end