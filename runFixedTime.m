%% Low Velocity Breakdown
%
% In this script, we study low velocity breakdown of different operating
% modes in a Stark Decelerator
%

%% First for S=1, SF Mode:
ivz = {1000,950,900,850,801,781,761,741,721,701,681,660.7,650.7,640.7,630.7,620.7};
dcs1 = struct('a','longdecel','b','longdecel');
dcsf = struct('a','longdecel','b','singlerod');

for n = [423,393,363,333,304,292,280,268,256,244,232,220,214,208,202,196]
    ct{i} = repmat('ba',1,n);
    rot = [0 0 90 90 180 180 270 270];
    rot = repmat(rot,1,ceil(n/4));
    rot{i} = rot(1:2*n);

    trans = [1 1 0 0];
    trans = repmat(trans,1,ceil(n/2));
    trans{i} = trans(1:2*n);

    ep{i} = repmat([125 235 305 55],1,n); % S=1
    %ri.endphases{i} = [repmat([125 235 305 55],1,n)]; % SF
    %ri.endphases{i} = [repmat([145 234.3 325 54.3],1,n)]; % VSF
    %ri.endphases{i} = [repmat([150 229.35 330 49.35],1,n)]; % XSF
end

rss1 = simdecel('initvz',ivz,'decels',dcs1,'chargetype',ct,'rot',rot,...
    'trans',trans,'endphases',ep,'calctype',repmat('p',1,1000));
rssf = simdecel('initvz',ivz,'decels',dcsf,'chargetype',ct,'rot',rot,...
    'trans',trans,'endphases',ep,'calctype',repmat('p',1,1000));

%% Next for VSF: