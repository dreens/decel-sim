function r = inittrap(r)

    % Trapping is beginning so loading has finished:
    r.loadingofftime = r.time;
    
    % save the loading zmax before the fields are updated to loading.
    zmaxload = r.f.zmaxload;

    % Load the mat file, or generate it from a COMSOL .dat file if it
    % doesn't exist yet.
    if ~r.reloadfields && exist(['../Traps/' r.trapname '.mat'],'file');
        r.f = load(['../Traps/' r.trapname '.mat']);
    elseif exist(['../Traps/' r.trapname '.dat'],'file');
        processtrappingfields(r);
        r.f = load(['../Traps/' r.trapname '.mat']);
    else
        error(['File ''slowANDtrap/Traps/' r.trapname '.dat'' not found']);
    end
    
    % This way you can set reloadfields and it will reload once for a
    % collection of runs, but not each time.
    r.reloadfields = false;
    
    % The assumption is made that the trapping fields have z=0 defined in
    % the way that is most sensible for this to be done. It is also assumed
    % that the loading fields are defined up to the same geometric boundary
    % in the positive z direction as are the loading fields. Thus it is
    % possible to infer the coordinate shift between them based on the z
    % coordinates of their respective endpoints.
    r.trapzero = r.loadzero - (r.f.zmaxtrap - zmaxload);
    r.offsetz = r.trapzero;
    
    % Nonetheless, lets check that there is actually a minimum near
    % r.trapzero. If it is too far away, we'll make a shift.
%     zcheck = -10e-3:1e-6:10e-3;
%     [minima, loczero] = min(r.f.vf([zeros(length(zcheck),2) zcheck'],1));
%     if abs(loczero)>1e-3
%         r.trapzero = r.trapzero + loczero;
%         r.offsetz = r.trapzero;
%         disp(['\nWarning, trapping fields found not to have minimum',...
%             ' at z=0, trap center shifted\n'])
%     end
    
    % Knowing the trapping zero in terms of simulation coordinates (z=0 at
    % beginning of the decelerator) we transform the molecules into the
    % loading frame:
    r.pos(:,3) = r.pos(:,3) - (r.trapzero - r.loadzero);
    
    % Undo any voltage scaling during loading.
    r.xtrascale = 1;
    
    % Update the data before beginning trapping
    r = update(r);
end
