classdef mixedExpDecay 
	% mixedExpDecay Generate different exp decay with the rate of
    % decay (<=s)
    %
    % mixedExpDecay Properties:
    %   s - rate of decay 
    %   c - weight
    %
    % mixedExpDecay Methods:
    %   init - initialize all properties
    %   gen - generate moments (m_1, ..., m_k) or m_k
    properties  
        c
        s
    end
    methods
        function obj = init(obj,s,c ) 
            c = reshape(c,1,length(c)); 
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
                if length(obj.c) == 1
                    idx = reshape(idx, size(idx,1)*size(idx,2),1);
                    moment = exp(-obj.s .*idx);
                else
                    sn = obj.s;
                    if length(idx) == 1
                        moment = obj.c*exp(-sn'*idx);
                    else
                        idx = reshape(idx, size(idx,1)*size(idx,2),1);
                        C = repmat(obj.c,length(idx), 1);
                        moment = sum(C.*exp(-sn.*idx),2);
                    end
                end
            else
                obj.s = sym(obj.s);
                obj.c = sym(obj.c);
                if length(obj.c) == 1
                    idx = reshape(idx, size(idx,1)*size(idx,2),1);
                    moment = vpa(exp(-obj.s .*sym(idx)));
                else
                    sn = (obj.s);
                    if length(idx) == 1
                        moment = vpa(obj.c*exp(-sn'*sym(idx)));
                    else
                        idx = reshape(idx, size(idx,1)*size(idx,2),1);
                        C = repmat(obj.c,length(idx), 1);
                        moment = vpa(sum(C.*exp(-sn.*sym(idx)),2));
                    end
                end
            end
        end
    end
end