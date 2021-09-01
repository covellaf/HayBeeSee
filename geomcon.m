function [c,ceq] = geomcon(theta,L,n,delx,dely,nrigid)

dL = L/n;

c = [];
ceq = [dL*sum(cos(theta))-dely;
       dL*sum(sin(theta))-delx];
%        zeros(1,n-nrigid),[ones(1,nrigid)].*theta];