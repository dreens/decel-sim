%% 
% Show Phase Space, various Op Modes
%
%function modeComparePhaseSpace()



%% For each mode, get corrected initial velocities to fix runtime.
modes = {'s=1','s=3','sf','vsf','xsf'};
decels = struct('a','longdecel','b','longdecel');
decels(2:length(modes)) = decels(1);
decels(3).b = 'singlerod';
decels(4).b = 'ppgg';
decels(5).b = 'ppmm_2mm';
phases = [125 235 305 55 ; 125 235 305 55 ; 125 235 305 55 ; 145 234.3 325 54.3 ; 150 229.35 330 49.35];

if ~exist('rsall')

    % make a dummy call to initialize struct array:
    demo = callSim(20,1000,decels(1),phases(1,:),false,'num',1);
    demo(1:length(modes)) = demo;
    rsall = demo;

    for j=1:length(modes)

        m = modes{j};
        n = 333;
        vip = 850;

        % First get the acceleration
        iss3 = strcmp('s=3',m);

        % Now do a serious full run.
        out = callSim(n,vip,decels(j),phases(j,:),iss3,'num',1e6);
        fprintf('Mode: %s, N: %d, V: %.1f, T: %.3f, num: %d, vf: %.1f\n',...
            m,n,vip,out.time*1000,out.numleft,out.vels(end));

        rsall(j) = out;
    end

    
end
    
%% Plotting
%
% Compared to work I did studying various final speeds, the simulation part
% of this script is very straightforward. This plotting section is where
% the real action happens.
%
f = figure; %hold on
f.Position = [100 50 1000 400];
f.PaperPosition = [0 0 10 4];
[ha, pos] = tight_subplot(2,l,0,[.13 .01],[.08 .01]);
l = length(modes);
nn = 100;
big = rsall(strcmp('xsf',modes)).numleft;
for m=modes
    i = find(strcmp(m,modes));
    s = rsall(i);
    axes(ha(i));
    x = (s.pos(:,3)-s.pos(1,3))*1e3;
    y = s.vel(:,3)-s.vel(1,3);
    c = ksdensity([x(1:nn),y(1:nn)],[x,y]);
    c = c * s.numleft/big;
    c = log(c);
    scatter(x,y,2,c,'filled');
    ca = caxis;
    caxis([-7 ca(2)+log(big/s.numleft)])
    %cb = colorbar;
    %cb.Ruler.Scale = 'log';
    xlim([-1.7 1.7])
    ylim([-17 17])
    text(1.6,15,upper(m),'FontSize',13,'FontWeight','bold','HorizontalAlignment','right')
    set(gca,'FontSize',13)
    set(gca,'LineWidth',2)
    set(gca,'TickLength',[.02 .05])
    colormap(jet)
    
    axes(ha(i+l))
    x = s.pos(:,1)*1e3;
    y = s.vel(:,1);
    c = ksdensity([x(1:nn),y(1:nn)],[x,y]);
    c = c * s.numleft/big;
    c = log(c);
    scatter(x,y,2,c,'filled');
    ca = caxis;
    caxis([-7 ca(2)+log(big/s.numleft)])
    %cb = colorbar;
    %cb.Ruler.Scale = 'log';
    xlim([-1.7 1.7])
    ylim([-17 17])
    set(gca,'FontSize',13)
    set(gca,'LineWidth',2)
    set(gca,'TickLength',[.02 .05])
    colormap(jet)
end

set(ha(1:l),'XTickLabel','');
set(ha([2:l l+2:end]),'YTickLabel','');
axes(ha(l+ceil(l/2)))
xlabel('Position (mm)','FontSize',13)
axes(ha(1))
ylabel('Velocity (m/s)','FontSize',13)
axes(ha(1+l))
ylabel('Velocity (m/s)','FontSize',13)
%axes(ha(ceil(l/2)))
%title('Breakdown of Effective Trap','FontSize',14)
%end
print(gcf,'Figures/5x2-PSD-Compare','-dpng','-r300')
system('open Figures/5x2-PSD-Compare.png')
%close(gcf)

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
        'smallt',4e-7,...
        'calctype',repmat('p',1,1000),...
        varargin{:});



end




