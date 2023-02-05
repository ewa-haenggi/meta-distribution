% Algorithm for finding the infima and suprema of all possible cdfs from a moment sequence 
clear all;
tic;
% input: moment sequence and the point(s) of interest
% The sequence starts with m1
n = 5;
moment = 1./(1+(1:n));
pt = linspace(0,1,51);

 
methodVar = CMClass;
methodVar = methodVar.init(length(moment));

infs = 0*pt;
sups = 0*pt;

for i = 1:length(pt)
    bounds = methodVar.CMBounds(moment, pt(i));
    infs(i) = bounds(1);
    sups(i) = bounds(2);
end

% output: infs (infima) and sups (suprema)
plot(pt, infs,'-x');
hold on;
plot(pt, sups,'-o');
toc;
