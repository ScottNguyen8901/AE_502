clear all; close all; clc

% UNITS
% Distance: AU
% Time: Days

AUtokm = 1.496*10^8; daystosec = 86400; 
mu_sun = ((daystosec^2)/(AUtokm^3))*1.327*10^11;

r1Evec = [-1.796136509111975*10^-1 9.667949206859814*10^-1 -3.668681017942158*10^-5];
v1Evec = [-1.720038360888334*10^-2 -3.211186197806460*10^-3 7.927736735960840*10^-7];

r1Ovec = [3.515868886595499*10^-2 -3.162046390773074 4.493983111703389]
v1Ovec = [-2.317577766980901*10^-3 9.843360903693031*10^-3 -1.541856855538041*10^-2]

r1Bvec = [7.249472033259724 14.61063037906177 14.24274452216359];
v1Bvec = [-8.241709369476881*10^-3 -1.156219024581502*10^-2 -1.317135977481448*10^-2];

dep = 1:365;
arr = (7*30 + 1):(365*2);
r2Evec = zeros(length(dep),3);
r2Ovec = zeros(length(arr),3);

for i = 1:length(dep)
    [r2Evec, v2Evec] = keplerprop(r1Evec, v1Evec, dep(i), mu_sun);
    for j = 1:length(arr)
        [r2Ovec, v2Ovec] = keplerprop(r1Ovec, v1Ovec, arr(j), mu_sun);
        dt = j - i;
        [v1vec, v2vec] = lambertuniversal(r2Evec, r2Ovec, dt, mu_sun, 0)
        dv(j,i) = norm(v2vec - v2Ovec) + norm(v1vec - v2Evec);
        if dv(j,i) > 50
            dv(j,i) = NaN
        end
    end
end

figure(1)
surf(dep,arr,dv)


for i = 1:length(dep)
    [r2Evec, v2Evec] = keplerprop(r1Evec, v1Evec, dep(i), mu_sun);
    for j = 1:length(arr)
        [r2Ovec, v2Ovec] = keplerprop(r1Ovec, v1Ovec, arr(j), mu_sun);
        dt = j - i;
        [v1vec, v2vec] = lambertuniversal(r2Evec, r2Ovec, dt, mu_sun, 0)
        dv(j,i) = norm(v1vec - v2Evec);
        if dv(j,i) > 20
            dv(j,i) = NaN
        end
    end
end

figure(2)
surf(dep,arr,dv)

for i = 1:length(dep)
    [r2Evec, v2Evec] = keplerprop(r1Evec, v1Evec, dep(i), mu_sun);
    for j = 1:length(arr)
        [r2Bvec, v2Bvec] = keplerprop(r1Bvec, v1Bvec, arr(j), mu_sun);
        dt = j - i;
        [v1vec, v2vec] = lambertuniversal(r2Evec, r2Bvec, dt, mu_sun, 0)
        dv(j,i) = norm(v2vec - v2Bvec) + norm(v1vec - v2Evec);
        if dv(j,i) > 50
            dv(j,i) = NaN
        end
    end
end

figure(3)
surf(dep,arr,dv)


for i = 1:length(dep)
    [r2Evec, v2Evec] = keplerprop(r1Evec, v1Evec, dep(i), mu_sun);
    for j = 1:length(arr)
        [r2Bvec, v2Bvec] = keplerprop(r1Bvec, v1Bvec, arr(j), mu_sun);
        dt = j - i;
        [v1vec, v2vec] = lambertuniversal(r2Evec, r2Bvec, dt, mu_sun, 0)
        dv(j,i) = norm(v1vec - v2Evec);
        if dv(j,i) > 20
            dv(j,i) = NaN
        end
    end
end

figure(4)
surf(dep,arr,dv)
