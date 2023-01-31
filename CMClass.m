classdef CMClass
    % CMClass Use the average of the Chebyshev-Markov (CM)
    % bounds to reconstruct the distribution (cdf) from moment sequence
    % (m_1, ..., m_n)
    %
    % CMClass Properties:
    %   order - The highest order of moment that are accepted for
    %   reconstruction 
    %   dig - The variable precision
    %
    % CMClass Methods:
    %   init - Initialize the order and digits for precision
    %   CMBounds - Evaluate the CM bounds at the point of interest
    %   based on the moment sequence
    %   value - Evaluate approximate cdf at the point of interest

    properties
         order 
         dig = 16
    end
    methods
        function obj = init(obj,varargin) 
            if nargin == 3
                obj.order = varargin{1};
                obj.dig = varargin{2};
            elseif nargin == 2
                obj.order = varargin{1};
            end  
        end
        function [bounds] = CMBounds(obj,varargin)
            if nargin == 3
                moment = varargin{1};
                pt = varargin{2};
            end
            digits(obj.dig);
            moment = reshape(moment, size(moment,1)*size(moment,2),1);
            if length(moment) < obj.order
                fprintf('Not enough moment input\n');
                return
            end
            moment = moment(1:obj.order);
            eps = 1e-5;
            n = obj.order;
            pointOfInterest = pt;
            if mod(n,2) == 0
                m = n/2;
                % u_m wrt y dF shifted moment = (m_k)_{k = 1}^{2m}
                % Hankel matrix m-by-(m+1)
                % m_1, ..., m_{m+1}
                % ...
                % m_{m}, ..., m_{2m}
                shiftM = moment(1:2*m);
                c = shiftM(1:m);
                r = shiftM(m:2*m);
                hankelU = hankel(c,r);
                um = @(x) det([hankelU;power(x,0:m)]);
                syms x;
                symum = det([hankelU;power(x,0:m)]);  
                % v_m wrt (1-y) dF shifted moment = (m_k - m_{k+1})_{k = 0}^{2m-1}
                % Hankel matrix m-by-(m+1) (wrt to shifted moment)
                % m_1, ..., m_{m+1}
                % ...
                % m_{m}, ..., m_{2m}
                shiftM = [1; moment(1:2*m-1)] - moment(1:2*m);
                c = shiftM(1:m);
                r = shiftM(m:2*m);
                hankelV = hankel(c,r);
                vm = @(x) det([hankelV;power(x,0:m)]);
                %             syms x;
                symvm = det([hankelV;power(x,0:m)]);
                if sum(double(abs(pointOfInterest - double(solve(symum)))) <eps)>=1
                    %         ismember(pointOfInterest, solve(symum)) == 1
                    root = [0,solve(symum)'];  % actually 0 is enough
                elseif sum(abs(pointOfInterest - double(solve(symvm))) <eps)>=1
                    %         ismember(pointOfInterest, solve(symvm)) == 1
                    root = [solve(symvm)',1];
                elseif sign(um(pointOfInterest))*sign(vm(pointOfInterest)) == 1 
                    % omega of degree m wrt (y-poi) dF shifted moment =
                    % (m_k)_{k = 1}^{2m} - poi * (m_k)_{k = 0}^{2m-1}
                    % Hankel matrix m-by-(m+1) (wrt to shifted moment)
                    % m_1, ..., m_{m+1}
                    % ...
                    % m_{m}, ..., m_{2m}
                    shiftM = moment(1:2*m) - pointOfInterest*[1;moment(1:2*m-1)];
                    c = shiftM(1:m);
                    r = shiftM(m:2*m);
                    hankelOmega = hankel(c,r);
                    %                 syms x;
                    omega = det([hankelOmega;power(x,0:m)]);
                    root = [solve(omega,x)', pointOfInterest];
                elseif sign(um(pointOfInterest))*sign(vm(pointOfInterest)) == -1
                    % omega of degree m-1 wrt y(y-1)(y-poi) dF shifted moment =
                    % (m_k)_{k = 3}^{2m} - (1+poi)*(m_k)_{k = 2}^{2m-1} +
                    % poi*(m_k)_{k = 1}^{2m-2}
                    if m == 1
                        omega = 0*x+1;
                    else
                        % Hankel matrix (m-1)-by-m (wrt to shifted moment)
                        % m_1, ..., m_{m}
                        % ...
                        % m_{m-1}, ..., m_{2m-2}
                        shiftM = moment(3:2*m) - (1+pointOfInterest)*moment(2:2*m-1) + pointOfInterest*moment(1:2*m-2);
                        c = shiftM(1:m-1);
                        r = shiftM(m-1:2*m-2);
                        hankelOmega = hankel(c,r);
                        omega = det([hankelOmega;power(x,0:m-1)]);
                    end
                    root = [solve(omega,x)', pointOfInterest, 0, 1];
                end
            else
                m = (n+1)/2;
                if m > 1
                    % t_m wrt dF shifted moment = (m_k)_{k = 0}^{2m-1}
                    % Hankel matrix m-by-(m+1)
                    % m_1, ..., m_{m+1}
                    % ...
                    % m_{m}, ..., m_{2m}
                    shiftM = [1;moment(1:2*m-1)];
                    c = shiftM(1:m);
                    r = shiftM(m:2*m);
                    hankelT = hankel(c,r);
                    tm = @(x) det([hankelT;power(x,0:m)]);
                    syms x;
                    symtm = det([hankelT;power(x,0:m)]);
                    % w_{m-1} wrt y(1-y) dF shifted moment = (m_k)_{k = 1}^{2m-2} -
                    % (m_k)_{k = 2}^{2m-1}
                    % Hankel matrix (m-1)-by-m
                    % m_1, ..., m_{m}
                    % ...
                    % m_{m-1}, ..., m_{2m-2}
                    shiftM = moment(1:2*m-2) - moment(2:2*m-1);
                    c = shiftM(1:m-1);
                    r = shiftM(m-1:2*m-2);
                    hankelW = hankel(c,r);
                    wm = @(x) det([hankelW;power(x,0:m-1)]);
                    symwm = det([hankelW;power(x,0:m-1)]);
                    root_symwm = double(solve(symwm));
                    root_symtm = double(solve(symtm));
                    if sum(abs(pointOfInterest - root_symtm) <eps)>=1
                        root = [0,solve(symtm)',1];
                    elseif sum(abs(pointOfInterest - root_symwm) <eps)>=1
                        root = [0,solve(symwm)',1];
                    elseif sign(tm(pointOfInterest))*sign(wm(pointOfInterest)) == 1
                        % omega of degree m-1 wrt y(y-poi) dF shifted moment = (m_k)_{k = 2}^{2m-1} -
                        % poi*(m_k)_{k = 1}^{2m-2}
                        % Hankel matrix (m-1)-by-m
                        % m_1, ..., m_{m}
                        % ...
                        % m_{m-1}, ..., m_{2m-2}
                        shiftM = moment(2 :2*m-1) - pointOfInterest*moment(1:2*m-2);
                        c = shiftM(1:m-1);
                        r = shiftM(m-1:2*m-2);
                        hankelOmega = hankel(c,r);
                        %                 syms x;
                        omega = det([hankelOmega;power(x,0:m-1)]);
                        root = [solve(omega)', pointOfInterest, 0];
                    elseif sign(tm(pointOfInterest))*sign(wm(pointOfInterest)) == -1
                        % omega of degree m-1 wrt (y-1)(y-poi) dF shifted moment = (m_k)_{k = 2}^{2m-1} -
                        % (1+poi)*(m_k)_{k = 1}^{2m-2} + poi*(m_k)_{k = 0}^{2m-3}
                        % Hankel matrix (m-1)-by-m
                        % m_1, ..., m_{m}
                        % ...
                        % m_{m-1}, ..., m_{2m-2}
                        shiftM = moment(2:2*m-1) - (1+pointOfInterest)*moment(1:2*m-2) + pointOfInterest*[1;moment(1:2*m-3)];
                        c = shiftM(1:m-1);
                        r = shiftM(m-1:2*m-2);
                        hankelOmega = hankel(c,r);
                        %                 syms x;
                        omega = det([hankelOmega;power(x,0:m-1)]);
                        root = [solve(omega)', pointOfInterest, 1];
                    end
                else
                    root = (pointOfInterest<= moment(1))*[pointOfInterest, 1] + ...
                        (pointOfInterest > moment(1)) * [0, pointOfInterest];
                end
            end
            root = double(root);
%             theoretically roots are all real, but (numerically) solving polynomials may
%             return complex values
            root = real(root);
            root = root(logical((0<=root).*(root<= 1))); 
            root = unique(root);
            root_matrix = power(root,(0:length(root)-1)');
            %             p = root_matrix\[1;moment(1:length(root)-1)];
            p = linsolve(root_matrix, [1;moment(1:length(root)-1)]);
            if ~isreal(p)
                p = real(p);
            end
%            find the index for the point of interest (it is ok to do so
%            b/c the point of interest mush be contained in the roots)
            idxarray = abs(root - pt); 
            idx = find(idxarray == min(idxarray),1);
            if idx>1
                minValue = sum(p(1:idx-1));
                maxValue = sum(p(1:idx));
            else
                minValue = 0;
                maxValue = p(1);
            end
            bounds = [minValue; maxValue];
        end
        function value = value(obj, varargin) %initialize
            if nargin == 3
                moment = varargin{1};
                pt = varargin{2};
            end
            bounds = obj.CMBounds(moment, pt);
            value = 1/2*sum(bounds);
            value = real(value);
        end 
    end
end