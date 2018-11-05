%% Creating Partials

cd ~/Documents/MATLAB/decel-sim/Partials/
load('endBetweenLast4X.mat')

for i=1:16
    r = rsf(i);
    first = '';
    if r.numleft < 1200
        first = 'W0S1';
    elseif r.numleft < 5000
        first = 'W0SF';
    else
        first = 'W0VSF';
    end
    first = 'W0XSF';
    p = r.phase;
    second = sprintf('%dp%2d',fix(p),round(100*(p-fix(p))));
    second = strrep(second,' ','0');
    save([first second],'r')
end