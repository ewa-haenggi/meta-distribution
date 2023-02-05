function [maxdv, avgdv, maxdh, avgdh] = dist2d(x1,y1,x2,y2)
% DIST2D calculate the difference of two sets of points according to 100000
% sampled points 
%
%   [MAXDV,AVGDV,MAXDH,AVGDH] = DIST2D(X1,Y1,X2,Y2) calculate the
%   difference 
n = 100000;
y1 = 1-y1;
y2 = 1-y2;
figure(3);
plot(x1, y1,'-or', x2, y2,'-xb');
if ~isempty(x2) 
    [dv,dh] = ccdfDiff(x1,y1, x2,y2, n);
    avgdv = sum(abs(dv))/(n-1);
    maxdv = max(abs(dv));
    avgdh = sum(abs(dh))/(n-1);
    maxdh = max(abs(dh));
else
    fprintf("No solutions.\n");
    avgdv = 1;
    maxdv = 1;
    avgdh = 1;
    maxdh = 1; 
end