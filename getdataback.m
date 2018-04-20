%% Get data back from a figure.
listC = get(gca,'Children');
for i=1:length(listC)
    tc = listC(i);
    name = tc.DisplayName;
    alldata(i).new = contains(name,'new');
    zlocs = strfind(name,'0');
    mult = 10^length(zlocs);
    alldata(i).speed = str2num(name(min(zlocs)-1))*mult;
    if isempty(zlocs)
        alldata(i).speed = 37;
    end
    alldata(i).XData = tc.XData;
    alldata(i).YData = tc.YData;
    alldata(i).UData = 2*tc.YNegativeDelta;
end
save('speedvarydata','alldata')