function [oe] = rv2oe(r0,v0,mu)
    
    %all angles calculated in degrees
    
    %i,j,k vectors
    I = [1 0 0];
    J = [0 1 0];
    K = [0 0 1];
    
    %Norm of r0 and v0
    r = norm(r0);
    v = norm(v0);
    
    %Calculating semi-major axis using vis-viva eqn
    a = 1/((2/r) - ((v^2)/mu));
    
    %Calculating eccentricity vector then eccentricity
    evec = (((v^2)/mu) - (1/r))*r0 - (1/mu)*dot(r0,v0)*v0;
    e = norm(evec);
    
    %Solving for inclination
    hvec = cross(r0,v0);
    h = norm(hvec);
    i = acosd(dot((hvec/norm(hvec)),K));
    
    %Solving for RAN
    nvec = cross(K,hvec);
    n = norm(nvec);

    if dot(nvec,J) < 0
        W = acosd(dot((nvec/norm(nvec)),I));
        W = 360 - W;
    else
        W = acosd(dot((nvec/norm(nvec)),I));
    end
    
    %Solving for argument of periapsis
    if dot(evec,K) < 0
        w = acosd(dot((nvec/norm(nvec)),(evec/norm(evec))));
        w = 360 - w;
    else
        w = acosd(dot((nvec/norm(nvec)),(evec/norm(evec))));
    end
    
    %Solving for true anomaly
    if dot(r0,v0) < 0
        f = acosd(dot((r0/norm(r0)),(evec/norm(evec))));
        f = 360 - f;
    else
        f = acosd(dot((r0/norm(r0)),(evec/norm(evec))));
    end
    
    oe = [a e i W w f]';
    
end