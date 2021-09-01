function [En] = elasten(theta,L,r,E,n,I)

dL = L/n;
% % I = (0.8e-3^3)*2e-3/12;
% I = 0.25*pi*r^4;
k = E*I/dL;

En = sum(0.5*k*(diff(theta)).^2);
