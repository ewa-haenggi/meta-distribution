classdef FLClass
    % FLClass Linear Transform version of Fourier Legendre (FL) method 
    % to reconstruct distribution (cdf) from moment sequence (m_1, ..., m_n)
    %
    % FLClass Properties:
    %   order - The highest order of moments that are accepted for
    %   reconstruction (highest order = order + 1 ), and the fl expansion
    %   is 0 - order. This is because we use the FL series to approximate
    %   the cdf instead of the pdf
    %   FourierC - The Fourier Coefficients, (order+1)-by-1 column vector
    %   dig - The variable precision
    %
    % FLClass Methods:
    %   init - Initialize the order and Fourier Coefficients from moment sequence (if no
    %   order is specified, order will be the highest order of the moment
    %   sequence minus 1)  
    %   value - Evaluate approximate cdf at the point of interest
    
    properties
        order 
        FourierC
        dig = 16
    end
    methods
        function obj = init(obj, varargin) % to obtain the Fourier Coefficients, 
            % which depend on the polynomial coefficients (0,0) here =
            % Legendre poly
            if nargin == 4
                moment = varargin{1};
                obj.order = varargin{2};
                obj.dig = varargin{3};
            elseif nargin == 3
                moment = varargin{1};
                obj.order = varargin{2};
            elseif nargin == 2
                moment = varargin{1};
                obj.order = length(moment) - 1;
            end
            moment = reshape(moment, size(moment,1)*size(moment,2),1); 
            if length(moment) - obj.order < 1
                fprintf('Not enough moment input. Automatically change the order of FL.\n');
                obj.order = length(moment)-1;
            end   
            
            n = obj.order+1; 
            order = obj.order;
            filename = sprintf("FLA_%d.mat", n); 
            if isfile(filename)
                loadfile = load(filename);
                coeffMatrix = loadfile.coeffMatrix;
            else
                digits(ceil(0.9*order));
                coeffMatrix = vpa(zeros(order+1, order+1)); 
                for k = 0:order % do it row by row 
                    for l = 0:k
                        j = vpa(0:(l));
                        k = vpa(k);
                        l = vpa(l);
                        coeffMatrix(k+1, l+1) = (2*k+1)*sum(gamma(k+1)./gamma(j+1)./gamma(k-j+1)...
                            .*gamma(k+1)./gamma(k-j+1)./gamma(j+1)...
                            .*gamma(k-j+1)./gamma(k-l+1)./gamma(l-j+1))*power(-1,k-l)./(l+1);
                    end
                end
                save(filename, 'coeffMatrix');
            end
            filename = sprintf("FLA1_%d.mat", n);
            if isfile(filename) 
                loadfile = load(filename);
                coeffMatrix1 = loadfile.coeffMatrix1;
            else
                digits(max(order,16));
                coeffMatrix1  = coeffMatrix*vpa(ones(order+1,1));
                save(filename, 'coeffMatrix1');
            end 
            digits(obj.dig);  
            xi = (coeffMatrix1 -coeffMatrix*moment);
            obj.FourierC = xi;  
        end 
        function v = value(obj, pt)  % pt can be scalar or vector
            digits(obj.dig);  
            v = R(obj.order,pt)*obj.FourierC; 
            function v = R(n,x)
                x = reshape(x, size(x,1)*size(x,2),1);
                v = sym(zeros(length(x), n+1));   
                for j = 0:n
                    if j == 0
                        v(:,1) = 1;
                    elseif j == 1
                        v(:,2) = 2*x-1;
                    else
                        v(:,j+1) = (2*j-1)/j*(2*x-1).*v(:,j) - (j-1)/j*v(:,j-1); % recurrence formula of Legendre poly
                    end 
                end 
            end
        end
         
    end
end
