classdef BMClass
    % BMClass Binomial Mixture method to reconstruct distribution (cdf)
    % from moment sequence (m_1, ..., m_n)
    %
    % BMClass Properties:
    %   AMatrix - The (n+1)-by-(n+1) A matrix 
    %   order - The highest order of moment that are accepted for
    %   reconstruction
    %   dig - The variable precision
    %
    % BMClass Methods:
    %   init - Initialize the A matrix and order
    %   bmTransform - Calculate the A matrix from a given order
    %   value - Evaluate approximate cdf based on moment sequence and the point of
    %   interest
    
    properties
        AMatrix 
        order 
        h
        dig = 16
    end
    methods
        function obj = init(obj,varargin) %initialize
            if nargin == 4 % the number of arguments including obj
                n = varargin{1};
                b = varargin{2}; 
                moment = varargin{3};
            elseif nargin == 3
                n = varargin{1};
                moment = varargin{2};
                b = 16;
            end 
                obj.order= n;
                obj.dig = b;  

                filename = sprintf("BMA_%d.mat", n);
                if isfile(filename)
                    loadfile = load(filename);
                    obj.AMatrix = loadfile.am;
                else
                    digits(ceil(n/2));
                    G=vpa(zeros(1,n+1));
                    for j=0:n
                        G(j+1)=vpa(gammaln(sym(j+1)));
                    end
                    am = vpa(zeros(n+1,n+1));
                    for i = 0:n
                        j = i:n;
                        am(i+1,j+1) = round(exp(G(n+1) - G(n-j+1) - G(i+1) - G(j-i+1))).* vpa(sym(-1).^(j-i));
                    end
                    obj.AMatrix = am;
                    save(filename, 'am');
                end
                digits(b)%digit requirement  
                moment = reshape(moment, size(moment,1)*size(moment,2),1);
                moment = [sym(1); moment(1:obj.order)];
                obj.h = obj.AMatrix * moment;  
        end 
        function v = value(obj, pt) %obtain value at a signle point or vector 
            digits(obj.dig);  
            n = obj.order; 
            cumsumH = cumsum([ obj.h]);
            v = 0*pt;
            int_j = floor(n*pt); 
            v = cumsumH(int_j +1);
            v(pt == 0) = 0; 
        end
    end
end


