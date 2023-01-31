classdef FCClass
    % FCClass Linear Transform version of Fourier Chebyshev (FC) method to reconstruct distribution (cdf) 
    % from moment sequence (m_1, ..., m_n)
    %
    % FCClass Properties:
    %   order - The highest order of moments that are accepted for
    %   reconstruction (highest order = order + 1 ), and the fl expansion
    %   is 0 - order. This is because we use the FL series to approximate
    %   the cdf instead of the pdf
    %   FourierC - The Fourier Coefficients
    %   dig - The variable precision
    %
    % FCClass Methods:
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
            % which depend on the polynomial coefficients (-1/2,-1/2) here
            % = Chebyshev poly 
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
            digits(obj.dig); 
            alpha = vpa(-1/2);
            beta = vpa(-1/2); 
            n = obj.order+1;
            order = obj.order; 
            coeffm = 2*power(gamma(vpa((0:order)')+1)./gamma(vpa((0:order)')+1/2),2);
            coeffm(1) = 1/vpa(pi);
            filename = sprintf("FCA_%d.mat", n); 
            if isfile(filename)
                loadfile = load(filename);
                coeffMatrix = loadfile.coeffMatrix;
            else
                digits(max(ceil(0.9*order),16));
                coeffMatrix = vpa(zeros(order+1, order+1)); 
                for k = 0:order % do it row by row 
                    for l = 0:k
                        j = vpa(0:(l));
                        k = vpa(k);
                        l = vpa(l);
                        coeffMatrix(k+1, l+1) = coeffm(k+1)*sum(gamma(k+alpha+1)./gamma(j+1)./gamma(k+alpha-j+1)...
                            .*gamma(k+beta+1)./gamma(k-j+1)./gamma(j+beta+1)...
                            .*gamma(k-j+1)./gamma(k-l+1)./gamma(l-j+1))*power(-1,k-l)./(l+1);
                    end
                end
                save(filename, 'coeffMatrix');
            end 
            filename = sprintf("FCA1_%d.mat", n);
            if isfile(filename) 
                loadfile = load(filename);
                coeffMatrix1 = loadfile.coeffMatrix1;
            else
                digits(max(order,16));
                coeffMatrix1  = coeffMatrix*vpa(ones(order+1,1));
                save(filename, 'coeffMatrix1');
            end 
            digits(obj.dig);  
            xi = vpa(coeffMatrix1 - coeffMatrix*moment);
            obj.FourierC = xi;   
        end 
        function v = value(obj, pt)
            digits(obj.dig);
            alpha = sym(-1/2);
            beta = sym(-1/2); 
            v = R(alpha, beta, obj.order, pt)*obj.FourierC; 
            weightFunc = 0*pt';% column vector
            for idx = 1:length(pt)
                if pt(idx)>0 && pt(idx)<1
                    weightFunc(idx) = power(pt(idx),alpha).*power(1-pt(idx), beta);
                end
            end 
            v = weightFunc.*v;
            v(abs(pt)<1e-16) = 0;
            v(abs(pt-1)<1e-16) = 1; 

            function v = R(a, b, n, x)
                x = reshape(x, size(x,1)*size(x,2),1);
                v = sym(zeros(length(x), n+1));   
                for j = 0:n
                    if j == 0
                        v(:, 1) = 1;
                    elseif j == 1
                        v(:, 2) = (a+1) + (a+b+2).*(x-1);
                    else
                        v(:,j+1) = 1/(2*j*(j+a+b)*(2*j+a+b-2))...
                            *((2*j+a+b-1)*((2*j+a+b)*(2*j+a+b-2)*(2.*x-1)+a^2-b^2)...
                            .*v(:,j) - 2*(j+a-1)*(j+b-1)*(2*j+a+b).*v(:,j-1));
                    end
                end
            end
        end

    end
end
