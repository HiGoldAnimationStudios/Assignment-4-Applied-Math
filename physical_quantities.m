function[E, H] = physical_quantities(V, orbit_params)
    x = V(:,1); y = V(:,2); vx = V(:,3); vy = V(:,4);
    mp = orbit_params.m_planet;
    ms = orbit_params.m_sun;
    G  = orbit_params.G;
    r  = hypot(x,y);
    v2 = vx.^2 + vy.^2;
    E = 0.5*mp*v2 - G*ms*mp ./ r;
    H = mp*(x.*vy - y.*vx);
end