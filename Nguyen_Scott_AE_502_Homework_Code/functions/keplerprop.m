function [r2vec, v2vec] = keplerprop(r1vec, v1vec, dt, mu)
    
    r1 = norm(r1vec);
    v1 = norm(v1vec);
    alph = (2/r1) - ((v1^2)/mu);
    
    [x, z] = kepleruniversal(r1vec, v1vec, dt, mu);
    [S,C] = stumpff(z);
    
    f = 1 - ((x^2)/r1)*C;
    g = dt - (1/sqrt(mu))*(x^3)*S;
    r2vec = f*r1vec + g*v1vec;
    r2 = norm(r2vec);
    
    fdot = (sqrt(mu)/(r2*r1))*(alph*(x^3)*S - x);
    gdot = 1 - ((x^2)/r2)*C;
    v2vec = fdot*r1vec + gdot*v1vec;

end


