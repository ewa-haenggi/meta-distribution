classdef MEClass
    % MEClass maximum entropy (ME) method to reconstruct distribution (cdf) 
    % from moment sequence (m_1, ..., m_n)
    %
    % MEClass Properties:
    %   order - The highest order of moment that are accepted for
    %   reconstruction
    %   lambda - The coefficient vector
    %   dig - The variable precision
    %
    % MEClass Methods:
    %   init - Initialize the order and Coefficients from moment sequence (if no
    %   order is specified, order will be the highest order of the moment
    %   sequence). Coefficients are obtained by solving the convex
    %   minimization problem
    %   value - Evaluate approximate cdf at the point of interest
    
    properties 
        order
        lambda
        dig = 16
    end
    methods
        function obj = init(obj, varargin) %initialize
            if nargin == 4
                moment = varargin{1};
                obj.order = varargin{2};
                obj.dig = varargin{3};
            elseif nargin == 3
                moment = varargin{1};
                obj.order  = varargin{2};  
            elseif nargin == 2
                moment = varargin{1};
                obj.order  = length(moment);
            end
            moment = reshape(moment, size(moment,1)*size(moment,2),1);
            moment = double(moment); % FMINUNC requires all values returned
            % by functions to be of data type double
            if length(moment) - obj.order < 0
                fprintf('Not enough moment input\n');
                return
            end  
            digits(obj.dig);
            n = length(moment);
            fun = @(lambda) lambda*moment(1:n) + log(integral(@(x) exp(-lambda*power(x,(1:n)')), 0, 1));
            x0 = ones(1,n);
            lambda1 = fminunc(fun,x0);
            lambda0 = log(integral(@(x) exp(-lambda1*power(x,(1:n)')), 0, 1));  
            obj.lambda = [lambda0, lambda1]; 
        end 
        function v  = value(obj,pt) 
            v = pt*0;
            pdf = @(x) exp(- obj.lambda*power(x,[0:length(obj.lambda)-1]'));
            for j = 1:length(pt)
                v(j) = integral(pdf, 0,pt(j) ); 
            end
        end
    end
end