function [v1vec, v2vec] = lambertuniversal(r1vec, r2vec, dt, mu, orbit)

    r1 = norm(r1vec); r2 = norm(r2vec); k = cross(r1vec,r2vec); k= k(3);

    if orbit == 0
        if k >= 0
            theta = acos(dot(r1vec,r2vec)/(r1*r2)); 
        elseif k < 0
            theta = 2*pi - acos(dot(r1vec,r2vec)/(r1*r2));
        end
    elseif orbit == 1
        if k < 0
            theta = acos(dot(r1vec,r2vec)/(r1*r2));
        elseif k >= 0
            theta = 2*pi - acos(dot(r1vec,r2vec)/(r1*r2));
        end
    
    end
    
    A = sin(theta)*sqrt((r1*r2)/(1 - cos(theta)));
    
%     z_range = -0:.1:3; l = length(z_range);
%     
%     z = z_range(1);
%     [S, C] = stumpff(z);
%     y = r1 + r2 + A*((z*S - 1)/sqrt(C));
%     
%     dS = (1/(2*z))*(C - 3*S);
%     dC = (1/(2*z))*(1 - z*S - 2*C);
%     dy = (A/(2*C^(3/2)))*((1 - z*S)*dC + 2*(S + z*dS)*C);
%     
%     F = ((y/C)^(3/2))*S + A*sqrt(y) - sqrt(mu)*dt;
%     s1 = sign(F);
%     
%     for i = 1:l
%     
%         z = z_range(i);
%         [S, C] = stumpff(z);
%         y = r1 + r2 + A*((z*S - 1)/sqrt(C));
%         
%         dS = (1/(2*z))*(C - 3*S);
%         dC = (1/(2*z))*(1 - z*S - 2*C);
%         dy = (A/(2*C^(3/2)))*((1 - z*S)*dC + 2*(S + z*dS)*C);
%     
%         F = ((y/C)^(3/2))*S + A*sqrt(y) - sqrt(mu)*dt;
%         s2 = sign(F);
%     
%         if s1 ~= s2
%             index = i;
%             break
%         end
%     
%         s1 = s2;
%         s1
%         s2
%     end
%     
%     z0 = z_range(index); z = z0; 
    z0 = 1; z = z0;
    tol = 10^-8; rat = 1;
    
    n_max = 100; n = 1;

    while abs(rat) > tol
    
        [S, C] = stumpff(z);
        y = r1 + r2 + A*((z*S - 1)/sqrt(C));
        
        dS = (1/(2*z))*(C - 3*S);
        dC = (1/(2*z))*(1 - z*S - 2*C);
        dy = (A/(2*C^(3/2)))*((1 - z*S)*dC + 2*(S + z*dS)*C);
    
        F = ((y/C)^(3/2))*S + A*sqrt(y) - sqrt(mu)*dt;
        dF = 1/(2*sqrt(y*C^5))*((2*C*dS - 3*dC*S)*y^2 + (A*C^(5/2) +3*C*S*y)*dy);
    
        rat = F/dF;
    
        z = z - rat;
        
        n = n + 1;

        if n > n_max
            break
        end
    
    end
    
    y = r1 + r2 + A*(z*S - 1)/sqrt(C);
    f = 1 - y/r1;
    g = A*sqrt(y/mu);
    gdot = 1 - y/r2;
    
    v1vec = (1/g)*(r2vec - f*r1vec);
    v2vec = (1/g)*(gdot*r2vec - r1vec);

end