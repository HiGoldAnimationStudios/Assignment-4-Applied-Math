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
    t_range = linspace(0,30,100);
    V_list = compute_planetary_motion(t_range,V0,orbit_params);
    
    
    %axis equal; axis square;
    %axis([-20,20,-20,20])
    %hold on
    %plot(0,0,'ro','markerfacecolor','r','markersize',5);
    %plot(V_list(:,1),V_list(:,2),'k');


    tspan=[0,10];
    h_ref=0.2;
    
    %Forward Euler
    BT_struct=struct();
    BT_struct.A=[0];
    BT_struct.B=[1];
    BT_struct.C=[0];
    
    my_rate_func = @(t_in, V_in) gravity_rate_func(t_in,V_in,orbit_params);

    [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(@my_rate_func,tspan,V0,h_ref,BT_struct);

    
end



function orion_stuff()
    
end