clear all; close all; clc

% UNITS
% Distance: AU
% Time: Days

AUtokm = 1.496*10^8; daystosec = 86400; 
mu_sun = 1.327*10^11;

r1Evec = AUtokm*[-1.796136509111975*10^-1 9.667949206859814*10^-1 -3.668681017942158*10^-5];
v1Evec = (AUtokm/daystosec)*[-1.720038360888334*10^-2 -3.211186197806460*10^-3 7.927736735960840*10^-7];

r1Ovec = AUtokm*[3.515868886595499*10^-2 -3.162046390773074 4.493983111703389]
v1Ovec = (AUtokm/daystosec)*[-2.317577766980901*10^-3 9.843360903693031*10^-3 -1.541856855538041*10^-2]

r1Bvec = AUtokm*[7.249472033259724 14.61063037906177 14.24274452216359];
v1Bvec = (AUtokm/daystosec)*[-8.241709369476881*10^-3 -1.156219024581502*10^-2 -1.317135977481448*10^-2];

dep = 1:365;
arr = (7*30 + 1):(365*2);
r2Evec = zeros(length(dep),3);
r2Ovec = zeros(length(arr),3);

for i = 1:length(dep)
    [r2Evec(i,:), v2Evec] = keplerprop(r1Evec, v1Evec, daystosec*dep(i), mu_sun);
end

for j = 1:length(arr)
    [r2Ovec(j,:), v2Evec] = keplerprop(r1Ovec, v1Ovec, daystosec*arr(j), mu_sun);
end
% 
% plot3(r2Evec(:,1), r2Evec(:,2), r2Evec(:,3))
% hold on
% plot3(r2Ovec(:,1), r2Ovec(:,2), r2Ovec(:,3))


[r2O, v2Evec] = keplerprop(r1Ovec, v1Ovec, 240*daystosec, mu_sun)





