classdef GPClass
    % GPClass Gil-Peleaz (GP) method to reconstruct distribution (cdf)
    % from moment sequence m_{js}
    %
    % GPClass Properties:
    %   stepsz - The step size of the integral 
    %   lowlim - The lower integration limit 
    %   uplim - The upper integration limit
    %   imagMoment - The imaginary moments given the step size, lower limit
    %   and upper limit
    %
    % GPClass Methods:
    %   init - Initialize all properties
    %   calImagMoment - Calculate the imaginary moments given the step
    %   size, lower limit and upper limit
    %   value - Evaluate approximate cdf based on moment sequence and the point of
    %   interest
    properties
        stpsz
        lowlim
        uplim
        imagMoment
    end
    methods
        function obj = init(obj, varargin) %initialize
            if nargin == 5
                obj.stpsz = varargin{1};
                obj.lowlim = varargin{2};
                obj.uplim = varargin{3};
                momentClass = varargin{4};
                tmp = obj.calImagMoment(momentClass);
                obj.imagMoment = reshape(tmp,1,size(tmp,1)*size(tmp,2));
            end
        end 
        function imagMoment = calImagMoment(obj, momentClass) % its a class type or function handle, 
            % since there are some other parameters combined with the moments
            tmin = obj.lowlim;
            tmax = obj.uplim;
            step_size = obj.stpsz;
            num_points = floor((tmax - tmin)/step_size) + 1;
            pt = linspace(tmin, tmax, num_points); 
            tic;
            if ~strcmpi(class(momentClass), 'function_handle')
                imagMoment = momentClass.gen(1j*pt);
            else
                imagMoment = momentClass(1j*pt);
            end
            toc;
        end
        function value = value(obj, pt) % using imaginary momoents
            if pt == 0
                value = 0;
            else
            tmin = obj.lowlim;
            tmax = obj.uplim;
            step_size = obj.stpsz;
            num_points = floor((tmax - tmin)/step_size) + 1;
            s_set = linspace(tmin, tmax, num_points);
            v_set = imag(exp(-1j*log(pt).*s_set).*obj.imagMoment./s_set);
            s_set1 = s_set;
            if tmin == 1e-15 % limit
                s_set1(1) = 0;
            end
            area1 = trapz(s_set1, v_set);
            value = 1/2 - 1/pi*area1; 
            end
        end
    end
end
