%% Unpacker converts a struct of arraylists into an array of structs
%
% Example:
%
% a.b = 3;
% a.c = {4,5};
% a.d = {'asdf','aoeu'};
% aa = unpacker(a,'product');
%
% >> aa
% >> 1x4 struct array with fields: b,c,d
%
% >> aa(1)
%   struct with fields:
%     b: 3
%     c: 4
%     d: 'asdf'
%
% >> aa(2)
%   struct with fields:
%     b: 3
%     c: 5
%     d: 'asdf'
%
% >> aa(3)
%   struct with fields:
%     b: 3
%     c: 4
%     d: 'aoeu'
%
% >> aa(4)
%   struct with fields:
%     b: 3
%     c: 5
%     d: 'aoeu'

%% loop through all fields and unpack
function rs = unpacker(r,type)
    switch(type)
        case 'linear'
            fs = fields(r);
            l = 1;
            for i=1:length(fs)
                if iscell(r.(fs{i}))
                    l = length(r.(fs{i}));
                    break;
                end
            end
            
            rs = repmat(r,1,l);
            
            for i=1:length(fs)
                for j=1:l
                    tf = r.(fs{i});
                    if iscell(tf)
                        rs(j).(fs{i}) = tf{j};
                    else
                        rs(j).(fs{i}) = tf;
                    end
                end
            end
        case 'product'
            fs = fields(r);
            for i=1:length(fs)
                if ~iscell(r.(fs{i}))
                    r.(fs{i}) = {r.(fs{i})};
                    if i==1
                        continue
                    end
                end
                these = r.(fs{i});
                len = length(these);
                height = length(r.(fs{i-1}));
                for j=1:i-1
                    r.(fs{j}) = repmat(r.(fs{j}),1,len);
                end
                r.(fs{i}) = reshape(repmat(these,height,1),1,len*height);
            end
            structargs = [fs struct2cell(r)];
            structargs = structargs';
            rs = struct(structargs{:});
    end
end