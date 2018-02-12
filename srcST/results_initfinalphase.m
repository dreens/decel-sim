%% Plots results
function results_initfinalphase(rs)
    figure('position',[50,50,1300,600])
    subplot(2,2,1)
    colors = get(gca,'ColorOrder');
    colors = [colors ; colors ; colors ; colors];
    ms = 1;
    speedcap = 600;
    
    hold on
    for i=1:length(rs)
        r = rs(i);
%         if i==1
%             plot(r.posiall(:,3)*1e3,r.veliall(:,3),'b.',...
%                 'Color',colors(1,:),'MarkerSize',2);
%         end
        keep = r.vel(:,3) < speedcap;
        
        plot(r.posi(keep,3)*1e3,r.veli(keep,3),'b.','Color',...
            colors(i+1,:),'MarkerSize',ms);
    end
    
    xlabel('Z Position Initial (mm)','FontSize',12)
    ylabel('Z Velocity Initial (m/s)','FontSize',12)
    title('Phase Space Z','FontSize',14)
    set(gca,'FontSize',12)
    grid on
    
    subplot(2,2,2)
    hold on
    for i=1:length(rs)
        r = rs(i);
%         if i==1
%             plot(r.posiall(:,1)*1e3,r.veliall(:,1),'b.',...
%                 'Color',colors(1,:),'MarkerSize',2);
%         end
        keep = r.vel(:,3) < speedcap;
        
        plot(r.posi(keep,1)*1e3,r.veli(keep,1),'b.','Color',...
            colors(i+1,:),'MarkerSize',ms);
    end
    
    xlabel('X Position Initial (mm)','FontSize',12)
    ylabel('X Velocity Initial (m/s)','FontSize',12)
    title('Phase Space X','FontSize',14)
    set(gca,'FontSize',12)
    grid on
    
    subplot(2,2,3)
    hold on
    finalnums = zeros(size(rs));
    for i=1:length(rs)
        r = rs(i);
        keep = r.vel(:,3) < speedcap;
        finalnums(i) = sum(keep)
        plot(r.pos(keep,3)*1e3,r.vel(keep,3),'b.','Color',...
            colors(i+1,:),'MarkerSize',ms);
    end
    
    xlabel('Z Position Final (mm)','FontSize',12)
    ylabel('Z Velocity Final (m/s)','FontSize',12)
    title('Phase Space Z','FontSize',14)
    set(gca,'FontSize',12)
    grid on
    
    subplot(2,2,4)
    hold on
    for i=1:length(rs)
        r = rs(i);
        keep = r.vel(:,3) < speedcap;
        
        plot(r.pos(keep,1)*1e3,r.vel(keep,1),'b.','Color',...
            colors(i+1,:),'MarkerSize',ms);
    end
    
    xlabel('X Position Final (mm)','FontSize',12)
    ylabel('X Velocity Final (m/s)','FontSize',12)
    title('Phase Space X','FontSize',14)
    set(gca,'FontSize',12)
    grid on
    
    
    % Now Lets calculate the peak density
    % We need to define a region in which to count.
    % x,y between +/-4m/s and +/-100um.
    % z within 250um of 818mm, vz within 5m/s of 34m/s
%     dist = [.7745 .774];
%     for i=1:length(rs)
%         pos = rs(i).pos; vel = rs(i).vel;
%         xc = (pos(:,1)/2.5e-4).^2 + (vel(:,1)/4).^2 < 1;
%         yc = (pos(:,2)/2.5e-4).^2 + (vel(:,2)/4).^2 < 1;
%         zc = ((pos(:,3)-dist(i))/5e-4).^2 + ((vel(:,3)-35)/5).^2 < 1;
%         number = sum(xc & yc & zc)
%         volume = pi^3 * (2.5e-4 * 2.5e-4 * 5e-4) * (4 * 4 * 5)
%     end
end