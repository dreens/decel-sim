%% Make 1D plots that pertain to PSD
function results_1DPSD(rs)
    figure; hold on
    for i=1:length(rs)
        r = rs(i);
        levels = 0.5:0.01:2;
        numbers = zeros(size(levels));
        N = max(size(r.vel));
        for j=1:length(levels)
            veld = r.vel-repmat(r.vel(1,:),N,1);
            posd = r.pos-repmat(r.pos(1,:),N,1);
            veld = veld/10;
            posd = posd/1e-3;
            veld = sum(veld.^2,2);
            posd = sum(posd.^2,2);
            numbers(j) = sum(veld + posd < levels(j)^2);
        end
        psd = numbers./(1000*pi^3*levels.^6/6);
        subplot(1,2,1); hold on
        plot(levels,numbers)
        subplot(1,2,2); hold on
        plot(levels,psd)
    end




end