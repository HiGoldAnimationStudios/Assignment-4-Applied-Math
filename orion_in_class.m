function orion_in_class()
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

    my_rate=@(t_in,V_in)gravity_rate_func(t_in,V_in,orbit_params);

    DormandPrince = struct();
    DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    DormandPrince.B = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    DormandPrince.A = [0,0,0,0,0,0,0;
1/5, 0, 0, 0,0,0,0;...
3/40, 9/40, 0, 0, 0, 0,0;...
44/45, -56/15, 32/9, 0, 0, 0,0;...
19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];
   
    h_ref=0.01;
    tspan=[0,30];
    [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(myrate,tspan,V0,h_ref,DormandPrince);


figure()
    subplot(2,1,1);
    hold on;
    plot(t_range,V_list(:,1))
    plot(t_range,V_list(:,2))

    plot(t_list,X_list(:,1))
    plot(t_list,X_list(:,2))
    

    xlabel("time")
    ylabel("pos")

    subplot(2,1,3);
    hold on;
    plot(t_range,V_list(:,3))
    plot(t_range,V_list(:,4))

%local
    n_samples=60;
    h_ref_list=logspace(-6,1,n_samples)
    tr_error_list1=zeros(1,n_samples)
    tr_error_list2=zeros(1,n_samples)

    abs_diff_list=zeros(1,n_samples)

    for n=1:length_h_ref_list
        h_ref=h_ref_list(n)
        V_list=compute_planetary_motion(tspan(1)+h_ref,V0,orbit_params)
        [XB1,XB2,~]=explicit_RK_step_embedded(my_rate,tspan(1),V0,h_ref,BT_struct);
        abs_diff_list(n)=norm(V_list-V0);
        tr_error_list1(n)=norm(XB1-V_list);
        tr_error_list2(n)=norm(XB2-V_list);
    end

    filter_parmas=struct();
    filter_parmas.min_y_val=1e-13;
    filter_parmas.max_y_val=1e-6;

    [p1,k1]=loglog_fit(h_ref_list,tr_error_list1,filter_parmas);
    [p2,k2]=loglog_fit(h_ref_list,tr_error_list2,filter_params);

    figure();
    loglog(h_ref_list,tr_error_list1);
    hold on
    loglog(h_ref_list,tr_error_list2);
    loglog(h_ref_list,abs_diff_list);



%global
    n_samples=30;
    h_ref_list=logspace(-5,1,30);

    h_avg_list=zeros(1,n_samples);
    
    num_evals_list=zeros(1,n_samples);
    tr_error_list=zeros(1,n_samples);
    for n=1:length(h_ref_list)
        h_ref=h_ref_list(n);

        [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(myrate,tspan,V0,h_ref,DormandPrince);

        tr_error=norm(X_list(end,:)-V_list(end,:));
        tr_error_list(n)=tr_error;
        num_evals_list(n)=num_evals;
        h_avg_list(n)=h_avg;
    end

    filter_params=struct();

    filter_params.min_y_val=1e-10;
    filter_params.max_y_val=1;

    [p1,k1]=loglog_fit(h_avg_list,tr_error_list,varagin);
    [p2,k2]=loglog_fit(num_evals_list,tr_error_list,varagin);
    abs(p1)
    abs(p2)

    figure();
    loglog(h_avg_list,tr_error_list)
    loglog(h_avg_list,k1*h_avg_list.^p1)
    
    figure();
    loglog(num_evals_list,tr_error_list)
    loglog(num_evals_list,k2*num_evals_list.^p2)
end

function [XB, num_evals] = explicit_RK_step(rate_func_in,t,XA,h,BT_struct)
    
    K=zeros(length(XA),length(BT_struct.C)); 
    num_evals=0;

    for n=1:length(BT_struct.B)
        t_temp=t+BT_struct.C(n)*h;
        x_temp=XA+h*(K*BT_struct.A(n,:)');
        K(:,n)=rate_func_in(t_temp,x_temp);
        num_evals=num_evals+1;
    end

    XB=XA+h*(K*BT_struct.B');

end

function [XB1, XB2, num_evals] = explicit_RK_step_embedded(rate_func_in,t,XA,h,BT_struct)
    
    K=zeros(length(XA),length(BT_struct.C)); 
    num_evals=0;

    for n=1:length(BT_struct.B)
        t_temp=t+BT_struct.C(n)*h;
        x_temp=XA+h*(K*BT_struct.A(n,:)');
        K(:,n)=rate_func_in(t_temp,x_temp);
        num_evals=num_evals+1;
    end

    XB1=XA+h*(K*BT_struct.B(1,:)');

end

function [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration ...
(rate_func_in,tspan,X0,h_ref,BT_struct)
    %your code here
    N=ceil((tspan(2)-tspan(1))/h_ref);
    h_avg=(tspan(2)-tspan(1))/N;

    t_list=linspace(tspan(1),tspan(2),N+1);
    X_list=zeros(N+1,length(X0));

    X_list(1,:)=X0';
    num_evals=0;

    for m=1:N
        [XB, num_evals_temp] = explicit_RK_step(rate_func_in,t_list(m),X_list(m,:)',h_avg,BT_struct);
        num_evals=num_evals+num_evals_temp;
        X_list(m+1,:)=XB';
    end


end


function dVdt = gravity_rate_func(t,V,orbit_params)
    x=V(1);
    y=V(2);
    Vx=V(3);
    Vy=V(4);



    m_sun=orbit_params.m_sun;
    G=orbit_params.G;


    r=[x;y];
    r_mag=norm(r);


    a=-(m_sun*G)/(r_mag^3)*r;
    ax=a(1);
    ay=a(2);

    dVdt=[Vx,Vy,ax,ay];
    
end

%this function computes the orbit of a planet about a sun
%the sun is assumed to located at the origin (0,0)
%and motion is restricted to the x-y plane
%INPUTS:
%t_list: a list of times to compute the position & velocity of the planet
%V0: initial condition.  V0 is a column vector consisting
%   of the initial position and velocity of the planet:
%       V0 = [x(0); y(0); dx/dt(0); dy/dt(0)]
%orbit_params: a struct describing the system parameters
%   orbit_params.m_sun: mass of the sun
%   orbit_params.m_planet: mass of the planet
%   orbit_params.G: gravitational constant
%       Force = -m_planet*m_sun*G/r^2
%OUTPUTS:
%V_list: the state of the planet at each time in t_list
%   if t_list is a list, then V_list is a Nx4 MATRIX
%   where each ROW has the form [x_i,y_i,dx/dt_i,dy/dt_i]
%   corresponding to the values at the ith time
%   if t_list is a SCALAR (i.e. t_list = t),
%   then V_list is a COLUMN VECTOR of the form:
%   [x(t); y(t); dx/dt(t); dy/dt(t)]
%NOTES:
%This function needs all the other functions in this file to run
%I HIGHLY RECOMMEND JUST SAVING THIS FUNCTION IN ITS OWN FILE
%DON'T COPY AND PASTE INTO YOUR CODE! IT'S NOT WORTH IT!
%
%USAGE EXAMPLE:
%At the bottom of this file is a function usage_example()
%which shows how to use compute_planetary_motion(...)
%You can start from there and then tweak it.
function V_list = compute_planetary_motion(t_list,V0,orbit_params)
    m_sun = orbit_params.m_sun;
    m_planet = orbit_params.m_planet;
    G = orbit_params.G;

    i = sqrt(-1);

    x0 = V0(1);
    y0 = V0(2);
    vx0 = V0(3);
    vy0 = V0(4);

    Q1 = abs(x0*vy0-y0*vx0);
    Q2 = abs(1/2*(vx0^2+vy0^2)-m_sun*G/sqrt(x0^2+y0^2));

    tau = Q1/Q2;
    l_scale = sqrt(Q1^2/Q2);
    velocity_scale = sqrt(Q2);

    x0=x0/l_scale;
    y0=y0/l_scale;
    vx0 = vx0/velocity_scale;
    vy0 = vy0/velocity_scale;

    t_list = t_list/tau;


    G = m_sun*G/(Q1*sqrt(Q2));
    m_sun = 1;
    m_planet = 1;


    H0 = m_planet*(x0*vy0-y0*vx0);

    dArea_dt = .5*H0/m_planet;


    r0 = abs(x0+y0*i);
    theta0 = angle(x0+y0*i);

    rdot0 = (vx0*x0+vy0*y0)/r0;
    theta_dot0 = H0/(m_planet*r0^2);

    alpha = m_sun*m_planet^2*G/H0^2;

    u0=1/r0;

    M = [-cos(theta0),sin(theta0);sin(theta0),cos(theta0)];
    B = [u0-alpha;-rdot0*u0^2/theta_dot0];
    Q = M\B;

    A_polar = abs(Q(1)+Q(2)*i);
    phi = angle(Q(1)+Q(2)*i);

    num_pts = length(t_list);
    V_list = zeros(num_pts,4);
    

    if 0<=alpha-A_polar && alpha-A_polar<1e-4
        d = 2/(alpha+A_polar);
        h0 = sin(theta0+phi)/(alpha+A_polar*cos(theta0+phi));
        area0 = kepler_func_parabola(h0,d,0);
    else
        a = alpha/abs(alpha^2-A_polar^2);
        b = sqrt(1/abs(alpha^2-A_polar^2));
        c = A_polar/abs(A_polar^2-alpha^2);

        if alpha>A_polar
            ellipse_area = pi*b*a;
            area0 = kepler_func_ellipse(theta0+phi,a,b,0);
        else
            h0 = (b^2/(a-c*cos(theta0+phi)))*sin(theta0+phi);
            area0 = kepler_func_hyperbola(h0,a,b,0);
        end
    end

    for n = 1:num_pts
        t = t_list(n);

        if 0<=alpha-A_polar && alpha-A_polar<1e-4
            target_area = -t*dArea_dt+area0;

            error_func = @(h_in) kepler_func_parabola(h_in,d,target_area);
            h_out = mini_secant_method(error_func,0);

            l_out = .5*h_out.^2/d -d/2;

            x_out = cos(-phi)*l_out-sin(-phi)*h_out;
            y_out = sin(-phi)*l_out+cos(-phi)*h_out;

            theta_out = angle(x_out+i*y_out);
            r_out = abs(x_out+i*y_out);
        else
            if alpha>A_polar
                target_area = t*dArea_dt+area0;
                target_area = mod(target_area,ellipse_area);
    
                if target_area>ellipse_area/2
                    target_area = target_area-ellipse_area;
                end
    
                error_func = @(theta_in) kepler_func_ellipse(theta_in,a,b,target_area);
                theta_out = mini_secant_method(error_func,0)-phi;
    
                r_out =  1/(alpha-A_polar*cos(theta_out+phi));

                x_out = r_out*cos(theta_out);
                y_out = r_out*sin(theta_out);
            else
                target_area = -t*dArea_dt+area0;

                error_func = @(h_in) kepler_func_hyperbola(h_in,a,b,target_area);
                h_out = mini_secant_method(error_func,0);
    
                l_out = a*sqrt(1+h_out^2/b^2)-c;
                
    
                x_out = cos(-phi)*l_out-sin(-phi)*h_out;
                y_out = sin(-phi)*l_out+cos(-phi)*h_out;
    
                theta_out = angle(x_out+i*y_out);
                r_out = abs(x_out+i*y_out);
            end
        end

        dtheta_dt_out = H0/(m_planet*r_out^2);
        drdt_out = -A_polar*sin(theta_out+phi)*dtheta_dt_out/(alpha-A_polar*cos(theta_out+phi))^2;

        vx_out = drdt_out*cos(theta_out)-dtheta_dt_out*r_out*sin(theta_out);
        vy_out = drdt_out*sin(theta_out)+dtheta_dt_out*r_out*cos(theta_out);

        V_list(n,:)=[x_out,y_out,vx_out,vy_out];
    end


    V_list = [V_list(:,1)*l_scale, V_list(:,2)*l_scale,...
            V_list(:,3)*velocity_scale, V_list(:,4)*velocity_scale];

    if isscalar(t_list)
        V_list = V_list';
    end

end

function error_val = kepler_func_parabola(y,d,target_area)
    error_val = y^3/(12*d)+y*d/4 - target_area;
end

function error_val = kepler_func_hyperbola(y,a,b,target_area)
    c = sqrt(a^2+b^2);
    q = y/b;
    hyperbolic_angle = log(q+sqrt(1+q^2));

    error_val = y*c/2-a*b*hyperbolic_angle/2 - target_area;
end

function error_val = kepler_func_ellipse(theta,a,b,target_area)
    i = sqrt(-1);

    c = sqrt(a^2-b^2);
    r = b^2/(a-c*cos(theta));

    x = (r*cos(theta)-c)/a;
    y = (r*sin(theta))/b;
    psi = angle(x+i*y);

    error_val = a*b*psi/2 + (c/2)*r*sin(theta) - target_area;
end

function x_out = mini_secant_method(func_in,xa)
    my_tol = 1e-13;

    xb = xa+1e-2;

    ya = func_in(xa);
    yb = func_in(xb);

    count = 0;
    while count<20 && abs(yb)>my_tol
        count = count+1;

        xc = (xa*yb-xb*ya)/(yb-ya);
        yc = func_in(xc);

        xa = xb;
        ya = yb;

        xb = xc;
        yb = yc;
    end

    x_out = xb;
end

%example for how to use compute_planetary_motion(...)
function usage_example()
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
    
    
    axis equal; axis square;
    axis([-20,20,-20,20])
    hold on
    plot(0,0,'ro','markerfacecolor','r','markersize',5);
    plot(V_list(:,1),V_list(:,2),'k');
end

function stuff()
    N=ceil((tspan(2)-tspan(1))/h_ref);
    h_avg=(tspan(2)-tspan(1))/N;

    t_list=linspace(tspan(1),tspan(2),N+1);
    X_list=zeros(N+1,length(X0));

    X_list(1,:)=X0';
    num_evals=0;

    for m=1:N
        [XB, num_evals_temp] = explicit_RK_step(rate_func_in,t_list(m),X_list(m,:)',h_avg,BT_struct);
        num_evals=num_evals+num_evals_temp;
        X_list(m+1,:)=XB';
    end
end