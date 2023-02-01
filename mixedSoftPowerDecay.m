classdef mixedSoftPowerDecay 
	% mixedSoftPowerDecay Generate different soft power decay with the same rate of
    % decay (s)
    %
    % mixedSoftPowerDecay Properties:
    %   s - rate of decay
    %   a,b - other decay parameter
    %   c - weight
    %
    % mixedSoftPowerDecay Methods:
    %   init - initialize all properties
    %   gen - generate moments (m_1, ..., m_k) or m_k
    properties 
        a
        b
        c
        s
    end
    methods
        function obj = init(obj,s,a,b,c)
            a = reshape(a,1,length(a));
            b = reshape(b,1,length(b));
            c = reshape(c,1,length(c));
            obj.a = a;
            obj.b = b;
            obj.c = c;
            obj.s = s;
        end
        function moment = gen(obj,varargin) 
            if nargin == 2
                idx = varargin{1};
                dig = 16;
            elseif nargin == 3
                idx = varargin{1};
                dig = varargin{2};
            end 
            digits(dig);
            obj.c = obj.c/sum(obj.c);
            if dig <= 16
                if length(idx) == 1
                    moment = sum(obj.c.* power(obj.a,obj.s).*power(idx + obj.a, -obj.s).*obj.b./(log(idx + 1) + obj.b));
                else
                    idx = reshape(idx, size(idx,1)*size(idx,2),1);
                    C = repmat(obj.c,length(idx), 1);
                    moment = sum(C.*repmat(power(obj.a,obj.s).*obj.b, length(idx),1) .* power(idx + obj.a, -obj.s) ./(log(idx+1)+obj.b),2);
                end
            else
                obj.a = sym(obj.a);
                obj.b = sym(obj.b);
                obj.c = sym(obj.c);
                obj.s = sym(obj.s);
                if length(idx) == 1
                    moment = vpa(sum(obj.c.* power(obj.a,obj.s).*power(sym(idx) + obj.a, -obj.s).*obj.b./(log(sym(idx) + 1) + obj.b)));
                else
                    idx = reshape(idx, size(idx,1)*size(idx,2),1);
                    C = repmat(obj.c,length(idx), 1);
                    moment = vpa(sum(C.*repmat(power(obj.a,obj.s).*obj.b, length(idx),1) .* power(sym(idx) + obj.a, -obj.s) ./(log(sym(idx)+1)+obj.b),2));
                end
            end
        end
    end
end