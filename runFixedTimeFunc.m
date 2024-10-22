%% Low Velocity Breakdown
%
% In this script, we study low velocity breakdown of different operating
% modes in a Stark Decelerator
%


function rsall = runFixedTimeFunc()

%% For each mode, get corrected initial velocities to fix runtime.
modes = {'s3','s1','sf','vsf','xsf'};
decels = struct('a','longdecel','b','longdecel');
decels(2:length(modes)) = decels(1);
decels(3).b = 'singlerod';
decels(4).b = 'ppgg';
decels(5).b = 'ppmm_2mm';
phases = [125 235 305 55 ; 125 235 305 55 ; 125 235 305 55 ; 145 234.3 325 54.3 ; 150 229.35 330 49.35];
nums = fliplr([190:2:230 235:5:290 300:10:400]);
times = repmat(nums,length(modes),1);
ivels = times;

% make a dummy call to initialize struct array:
demo = callSim(20,1000,decels(1),phases(1,:),false,'num',1);
demo(1:length(modes),1:length(nums)) = demo;
rsall = demo;
demo2 = demo(1,:);

parfor i=1:length(modes)
    % First get the acceleration
    iss3 = strcmp('s3',modes{i});
    out = callSim(400,1000,decels(i),phases(i,:),iss3,'num',1);
    A = (1000 - out.vels(end))/out.time;
    thesetimes = nums;
    thesevels = nums;
    thesers = demo2;
    % Now for each decelerator length...
    for n = nums
        
        % Guess vi assuming constant acceleration
        D = n*5e-3;
        T = 3e-3;
        vavg = D/T;
        vdelta = A*T;
        vi = vavg + vdelta/2;
        finalT = 0;
        
        % Now run with that vi and get the runtime.
        try
            out = callSim(n,vi,decels(i),phases(i,:),iss3,'num',1);
            finalT = out.time;
            vf = out.vels(end);
            thesetimes(nums==n) = finalT*1e3;
        catch E
            E
            vf = 0;
            finalT = 3.1e-3;
            thesetimes(nums==n) = 0;
        end
        
        % Then correct the initial velocity using classical mechanics. This
        % can be derived by using:
        % vi  + vf  = 2 D t, 
        % vi  - vf  = a t, 
        % vi' + vf' = 2 D t'
        % vi' - vf' = a t
        % Eliminate a and D and vf' and get vi' in terms of vi, vf, t/t'.
        dt = finalT/3e-3;
        vip = (vi*(dt+1/dt)+vf*(dt-1/dt))/2;
        thesevels(i,nums==n) = vip;
        % Now do a serious full run.
        out = callSim(n,vip,decels(i),phases(i,:),iss3,'num',1e5);
        fprintf('Mode: %s, N: %d, V: %.1f, T: %.3f, num: %d, vf: %.1f\n',...
            modes{i},n,vip,out.time*1000,out.numleft,out.vels(end));
        thesers(nums==n) = out;
    end
    times(i,:) = thesetimes;
    ivels(i,:) = thesevels;
    rsall(i,:) = thesers;
end
end

%%
function rs = callSim(N,ivz,D,M,iss3,varargin)
    ct = repmat('ba',1,N);
    rot = [0 0 90 90 180 180 270 270];
    rot = repmat(rot,1,ceil(N/4));
    rot = rot(1:2*N);

    tran = [1 1 0 0];
    tran = repmat(tran,1,ceil(N/2));
    tran = tran(1:2*N);
    
    if iss3
        rot = repmat(rot',1,3); 
        rot = rot';
        rot = rot(:)';
        rot = rot(1:2*N);
        tran = repmat(tran',1,3);
        tran = tran';
        tran = tran(:)';
        tran = tran(1:2*N);
    end
    
    % if the user hasn't asked for verbosity, turn it off
    if ~any(strcmp('verbose',varargin))
        varargin = [varargin {'verbose' false}];
    end

    ep = repmat(M,1,N); % S=1
    rs = simdecel('initvz',ivz,...
        'decels',D,...
        'chargetype',ct,...
        'rot',rot,...
        'trans',tran,...
        'endphases',ep,...
        'calctype',repmat('p',1,1000),...
        varargin{:});
end
