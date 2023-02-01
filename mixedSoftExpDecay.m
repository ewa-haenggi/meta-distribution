classdef mixedSoftExpDecay 
    % mixedSoftExpDecay Generate different soft exp decay with the same rate of
    % decay (s)
    %
    % mixedSoftExpDecay Properties:
    %   s - rate of decay
    %   a - other decay parameter
    %   c - weight
    %
    % mixedSoftExpDecay Methods:
    %   init - initialize all properties
    %   gen - generate moments (m_1, ..., m_k) or m_k
    properties 
        a
        c
        s
    end
    methods
        function obj = init(obj,s,a,c)
            a = reshape(a,1,length(a));
            c = reshape(c,1,length(c));
            obj.a = a;
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
                    moment = sum(obj.c.* exp(-obj.s * idx) .* obj.a ./(idx +obj.a));
                else
                    idx = reshape(idx, size(idx,1)*size(idx,2),1);
                    C = repmat(obj.c,length(idx), 1);
                    moment = sum(C.*  exp(-obj.s.*(idx)) .* repmat(obj.a,length(idx),1)./((idx) + obj.a) ,2);
                end
            else
                obj.a = sym(obj.a);
                obj.c = sym(obj.c);
                obj.s = sym(obj.s);
                if length(idx) == 1
                    moment = vpa(sum(obj.c.* exp(-obj.s * sym(idx)) .* obj.a ./(sym(idx) +obj.a)));
                else
                    idx = reshape(idx, size(idx,1)*size(idx,2),1);
                    C = repmat(obj.c,length(idx), 1);
                    moment = vpa(sum(C.* exp(-obj.s.*sym(idx)) .* repmat(obj.a,length(idx),1)./(sym(idx) + obj.a) ,2));
                end
            end
        end
    end
end