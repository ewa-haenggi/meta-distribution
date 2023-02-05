% Algorithm to find the number of moments needed for bounding the cdf
clear all;

% input: known moment, point of interest, and tolerance 
n = 100;
moment = 1./(1+(1:n)').^(1);
eps = 1e-1;

for i = 10:n % starts from 10 for simplicity
truncMoment = moment(1:i);
pt = linspace(0,1,51);
methodVar = CMClass;
methodVar = methodVar.init(length(truncMoment));

infs = 0*pt;
sups = 0*pt;

for j = 1:length(pt)
    bounds = methodVar.CMBounds(truncMoment, pt(j));
    infs(j) = bounds(1);
    sups(j) = bounds(2);
end
[x1,y1] = revisePtv3(pt, infs);
[x2,y2] = revisePtv3(pt, sups);
[maxdv, avgdv, maxdh, avgdh] = dist2d(x1,y1,x2,y2);
if avgdv < eps
    break
end
end

i
