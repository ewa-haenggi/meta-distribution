classdef FJClass
	% FJClass Fourier Jacobi (FJ) method to reconstruct distribution (cdf) 
    % from moment sequence (m_1, ..., m_n)
    %
    % FJClass Properties:
    %   order - The highest order of moments that are accepted for
    %   reconstruction (highest order = order), and FL expansion is
    %   0-order. This is because we use the FJ series to approximate the
    %   pdf
    %   alpha - The alpha coefficient for the Jacobi polynomials
    %   beta - The beta coefficient for the Jacobi polynomials
    %   FourierC - The a coefficients (vector) for the Fourier coefficients of
    %   the Jacobi polynomials
    %   dig - The variable precision 
    %
    % FJClass Methods:
    %   init - Initialize all the properties from moment sequence (if no
    %   order is specified, order will be the highest order of the moment
    %   sequence)
    %   initFourierC - Calculate the a coefficients (Fourier coefficients) of
    %   the Jacobi polynomials based on the moment sequence (from c_1 to c_n, no c_0)
    %   shiftedJacobi - Evaluate shifted Jacobi polunomial R_k at point pt
    %   value - Evaluate approximate cdf at the point of interest
    
    properties 
        order
        alpha
        be
        FourierC 
        dig = 16;
	end
    methods 
        function obj = init(obj, varargin)
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
            if length(moment) - obj.order < 0
                fprintf('Not enough moment input. Automatically change the order of FJ.\n');
                obj.order = length(moment);
            end  
            moment = moment(1:obj.order);  
            digits(obj.dig);
            m1 = moment(1);
            m2 = moment(2);
            a  =  double((m1-m2)*(1-m1)/(m2 - m1*m1)-1);
            b =  double((a +1)*m1/(1-m1)-1);   
            obj.alpha = a ;
            obj.be = b;
            h0 = 1/(a+b+1)*gamma(a+1)*gamma(b+1)/gamma(a+b+1);
            h1 = obj.initFourierC(moment);
            h1 = reshape(h1,1,size(h1,1)*size(h1,2)); 
            obj.FourierC = [1/h0, h1];
        end 
        function vector = initFourierC(obj, moment)
            a = obj.alpha;
            b = obj.be;
            ord = obj.order; 
            digits(obj.dig);
            vector = vpa(zeros(ord, 1));
            moment = reshape(moment, size(moment,1)*size(moment,2),1);
            moment = [1; moment(1:obj.order)];
            for n = 1: ord
                an = 0;
                L = 0:n;
                tmpGamma = vpa(gamma(sym(n+1))./...
                    (gamma(sym(L+1)).*gamma(sym(n - L +1))).*gamma(sym(n + a  + b  +1))./...
                    (gamma(sym(n + a  -L + 1)).*gamma(sym(b  +L +1))));
                for L  = 0:n
                    k = 0:n-L;
                    shifted_mom = vpa(sum(gamma(sym(n-L+1))./(gamma(sym(k+1)).*gamma(sym(n-L-k+1))).*power(-1,k).*moment(n-k+1)'));
                    an = an + tmpGamma(L+1)*shifted_mom;
                end
                vector(n) =  ((2*n+1+a  +b )*an);
            end
            
        end 
        function value = value(obj, pt)  
            digits(obj.dig);
            g = @(x) obj.FourierC*R( obj.alpha,  obj.be, obj.order, x);
            value = g(pt);
            function v = R(alpha, be, n, x)
                tm = zeros(n+1,1);
                for j= 0:n
                    if j == 0
                        tm(1) = betacdf(x,be+1,alpha+1)*beta(be+1,alpha+1);
                    else
                        tm(j+1) =  -1/j.*power(1-x,alpha+1).*power(x,be+1).*r(alpha+1, be+1, j-1, x);
                    end
                end
                v = tm;
                function v = r(alpha, be, n, x)
                    tmp = exp(gammaln((n+alpha +1)) - gammaln(((0:n)+1)) - gammaln((n+alpha ...
                        - (0:n)+1)) + gammaln((n+be+1)) - gammaln((n - (0:n)+1)) ...
                        - gammaln((be+ (0:n)+1))).*power(x-1, n - (0:n)).*power(x,(0:n));
                    v = sum(tmp);
                end
            end 

        end
        function value = valuepdf(obj, pt)
            digits(obj.dig);
            g = @(x) obj.FourierC*R( obj.alpha,  obj.be, obj.order, x)+...
                gamma(obj.alpha+obj.be+2)/gamma(obj.alpha+1)/gamma(obj.be+1).*power(1-x,obj.alpha).*power(x,obj.be);
            value = g(pt);
            function v = R(alpha, be, n, x)
                tm = zeros(n,1);
                for j= 1:n
                    tm(j) =  power(1-x,alpha).*power(x,be).*r(alpha, be, j, x); 
                end
                v = tm;
                function v = r(alpha, be, n, x)
                    tmp = exp(gammaln((n+alpha +1)) - gammaln(((0:n)+1)) - gammaln((n+alpha ...
                        - (0:n)+1)) + gammaln((n+be+1)) - gammaln((n - (0:n)+1)) ...
                        - gammaln((be+ (0:n)+1))).*power(x-1, n - (0:n)).*power(x,(0:n));
                    v = sum(tmp);
                end
            end
        end
        
    end
end
