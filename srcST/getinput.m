function r = getinput(r,args)
% To pass in variables directly, just use Name,Value pairing:
% simdeceltrap('initvz',500,'loadvelz',6)
%
% If preferred, just pass a struct containing the variables:
% initstruct.initvz = 500;
% initstruct.loadvelz = 6;
% simdeceltrap(initstruct)
    if ~isempty(args)
        inputs = length(args);
        if inputs==1 && isstruct(args{1})
            rinit = args{1};
            initfields = fields(rinit);
            for i=1:length(initfields)
                r.(initfields{i}) = rinit.(initfields{i});
            end
        elseif mod(inputs,2)
            error('Input a single struct or Name,Value pairs')
        else
            for i=1:inputs/2
                r.(args{2*i-1}) = args{2*i};
            end
        end
    end
end