%% Low Velocity Breakdown
%
% In this script, we study low velocity breakdown of different operating
% modes in a Stark Decelerator
%


function rsall = runFixedTimeFuncAllPar()

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
demo(1:length(modes)*length(nums)) = demo;
rsall = demo;

parfor i=1:length(rsall)
    
    j = floor(i/length(nums)-1e-9) + 1;
    m = modes{j};
    n = nums(1+mod(i-1,length(nums)));
    
    % First get the acceleration
    iss3 = strcmp('s3',m);
    out = callSim(400,1000,decels(j),phases(j,:),iss3,'num',1);
    A = (1000 - out.vels(end))/out.time;

    
    % Guess vi assuming constant acceleration
    D = n*5e-3;
    T = 3e-3;
    vavg = D/T;
    vdelta = A*T;
    vi = vavg + vdelta/2;
    finalT = 0;

    % Now run with that vi and get the runtime.
    try
        out = callSim(n,vi,decels(j),phases(j,:),iss3,'num',1);
        finalT = out.time;
        vf = out.vels(end);
    catch E
        E
        vf = 0;
        finalT = 3.1e-3;
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
    % Now do a serious full run.
    out = callSim(n,vip,decels(j),phases(j,:),iss3,'num',1e5);
    fprintf('Mode: %s, N: %d, V: %.1f, T: %.3f, num: %d, vf: %.1f\n',...
        m,n,vip,out.time*1000,out.numleft,out.vels(end));

    rsall(i) = out;
end

rsall = reshape(rsall,length(modes),length(nums));
%% Plotting
%
% Even for pathetically low speeds, I don't find much loss. I need to
% include some flight time out of the decelerator, which is universally
% required for use of the beam either for collisions or trapping.
%
if false
figure; hold on
for m=modes
    i = find(strcmp(m,modes));
    l = length(nums);
    s1 = [rsall(i,:).numleft];
    s1v = s1;
    dd = 5e-3;
    for j=1:l
        s = rsall(i,j);
        s1v(j) = s.vels(end);
        tt = dd/s1v(j);
        pf = s.pos + s.vel*tt;
        pf(:,3) = (pf(:,3) - pf(1,3))/2;
        pf = pf/1e-3;
        pf = sum(pf.^2,2);
        s1(j) = sum(pf<1.5);
        svz = sqrt(s.tempz*s.k/s.mOH);
        svx = sqrt(s.tempxy*s.k/s.mOH);
        psv1(j) = s1(j) / s.num * (svz*svx^2*s.spreadz*s.spreadxy^2);
    end
    plot(s1v,psv1,'DisplayName',m{:},'LineWidth',2)
end

xlabel('Final Speed','FontSize',13)
ylabel('Phase Space Volume (m^6/s^3)','FontSize',13)
title('Breakdown of Effective Trap','FontSize',14)
set(gca,'FontSize',13)
set(gca,'YScale','log')
h = legend('show');
set(h,'FontSize',13)
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
