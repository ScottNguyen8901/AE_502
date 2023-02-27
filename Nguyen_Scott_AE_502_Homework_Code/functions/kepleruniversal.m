function [x, z]  = kepleruniversal(r1vec, v1vec, dt, mu)
    
    r1 = norm(r1vec);
    v1 = norm(v1vec);
    vr1 = dot(r1vec,v1vec)/r1;
    alph = (2/r1) - ((v1^2)/mu);
    x0 = sqrt(mu)*abs(alph)*dt;
    
    z0 = alph*x0^2;
    x = x0; z = z0;
    rat = 1; tol = 10^-2;
    
    while abs(rat) > tol
        [S,C] = stumpff(z);
        f = ((r1*vr1)/sqrt(mu))*(x^2)*C + (1 - alph*r1)*(x^3)*S + r1*x - sqrt(mu)*dt;
        df = (r1*vr1)/(sqrt(mu))*x*(1 - alph*(x^2)*S) + (1 - alph*r1)*(x^2)*C + r1;
        rat = f/df;
        x = x - rat;
        z = alph*x^2;
    end

end

