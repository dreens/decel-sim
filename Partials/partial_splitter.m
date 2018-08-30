%% Creating Partials

cd ~/Documents/MATLAB/decel-sim/Partials/
load('seriousPartials.mat')

for i=1:48
    r = rsf(i);
    first = '';
    switch(fix((i-1)/16))
        case 0
            first = 'Bnorm';
        case 1
            first = 'BSF';
        case 2
            first = 'BVSF';
    end
    p = r.phase;
    second = sprintf('%dp%2d',fix(p),round(100*(p-fix(p))));
    second = strrep(second,' ','0');
    save([first second],'r')
end