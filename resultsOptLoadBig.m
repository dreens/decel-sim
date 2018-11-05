%% Optimize Trapping
% Plot first optimization in tricycle trap with new, sleek load and trap sim
function resultsOptLoadBig(rs)
    
    phase = zeros(size(rs));
    loadV = zeros(size(rs));
    decel = zeros(size(rs));
    trap = zeros(size(rs));
    tt = zeros(size(rs));
    ff = zeros(size(rs));
    for i=1:length(rs)
        phase(i) = str2double([rs(i).contname(end-1) rs(i).contname(end)]);
        loadV(i) = rs(i).endphases(end-1);
        switch(rs(i).decels.l{1})
            case 'tricycleload13broad'
                trap(i) = 1;
            case 'tricycleload25broad'
                trap(i) = 2;
            case 'ringloadbroad'
                trap(i) = 3;
            case 'cryoloadbroad'
                trap(i) = 4;
            case 'cryo2loadbroad'
                trap(i) = 5;
            case 'clover1load'
                trap(i) = 6;
            case 'clover2loadLF'
                trap(i) = 7;
            case 'clover2load'
                trap(i) = 8;
            case 'clover3load'
                trap(i) = 9;
            case 'ringloadbroadCLV'
                trap(i) = 10;
        end
        switch(rs(i).contfillstub)
            case 'XSF'
                decel(i) = 1;
            case 'VSF'
                decel(i) = 2;
            case 'SF'
                decel(i) = 3;
            case 'S1'
                decel(i) = 4;
        end
        rv = rs(i).vel;
        rv = rv(~isnan(rv));
        tt(i) = sum(rv(:).^2)*rs(i).mOH/(3*rs(i).k*rs(i).numleft);
        ff(i) = rs(i).molnum(2)/rs(i).molnum(1);
    end
    phase = sort(unique(phase));
    phase = 56+phase/100;
    loadV = sort(unique(loadV));
    %loadV = loadV(loadV<25);
    %loadV = loadV(loadV~=5e-3);
    [pp ll] = meshgrid(phase,loadV);
    trapnames = {'tri13','tri25','ring','cryo','cryo2','clv1','cl2L','clv2','clv3','mclv'};
    decelnames = {'XSF','VSF','SF','S1'};
    for i=1:length(trapnames)
        fprintf('\n')
        for j=1:length(decelnames)
            rr = rs(decel==j&trap==i);
            [m l] = max([rr.numleft]);
            t = tt(decel==j&trap==i);
            t = t(l);
            f = ff(decel==j&trap==i);
            f = f(l);
            fprintf('Max %d, Temp %.1f mK, Frac %.2f, %s, %s\n',...
                m,1000*t,f,trapnames{i},decelnames{j})
            if false
                figure;
                %surf(pp,ll,rr)
                N = reshape([rr.numleft],fliplr(size(pp)))';
                bar3(loadV,N)
                set(gca,'XTickLabel',phase)
                %set(gca,'XTickLabel',p)
                titlestr = sprintf('%s, %s, %s',...
                    'Optimize Speed and Loading v_F',...
                    trapnames{i},decelnames{j});
                title(titlestr)
                xlabel('Slowing Phase Angle (deg)')
                ylabel('Loading Speed (m/s)')
                zlabel('Final Population')

                figure;
                %surf(pp,ll,rr)
                bar3(loadV,1e3*reshape(tt(decel==j&trap==i),fliplr(size(pp)))')
                set(gca,'XTickLabel',phase)
                %set(gca,'XTickLabel',p)
                titlestr = sprintf('%s, %s, %s',...
                    'Optimize Temperature',...
                    trapnames{i},decelnames{j});
                title(titlestr)
                xlabel('Slowing Phase Angle (deg)')
                ylabel('Loading Speed (m/s)')
                zlabel('Temperature (mK)')
            end
        end
    end
    %{
    tt = reshape(tt,length(phase),length(loadV))';
    figure;
    %surf(pp,ll,rr)
    bar3(loadV,rr)
    set(gca,'XTickLabel',phase)
    %set(gca,'XTickLabel',p)
    title('Optimize Speed and Loading On-Time, Cryocycle VSF 90')
    xlabel('Slowing Phase Angle (deg)')
    ylabel('Loading Speed (m/s)')
    zlabel('Final Population')

    figure;
    %surf(pp,ll,rr)
    bar3(loadV,1e3*tt)
    set(gca,'XTickLabel',phase)
    %set(gca,'XTickLabel',p)
    title('Optimize Loading for Temperature, Cryocycle VSF 90')
    xlabel('Slowing Phase Angle (deg)')
    ylabel('Loading Speed (m/s)')
    zlabel('Temperature (mK)')

    %}


end