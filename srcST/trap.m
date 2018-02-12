function r = trap(r)
%% Step in trapping fields to see who lasts
% Eventually we could implement Majorana loss or other experimental
% sequences in here.
    fprintf('Trapping\n')
    for i=1:round(r.traptime/r.trapt)
        for j=1:round(r.trapt/r.smallt)
            r = smallstep(r,r.smallt);
        end
        
        r = update(r);
        
        if isnan(r.pos(1,3))
            fprintf('No Molecules Remaining\n');
            break
        end
                
        fprintf('Step %d,\t%dmols\n',i,r.numleft);
    end
end
