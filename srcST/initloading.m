function r = initloading(r)

    % Loading is beginning so deceleration is finished.
    r.decelofftime = r.time;
    
    % Now we need to make sure the molecule positions line up with
    % expectations for the new fields. z must be shifted by the full
    % decelerator length, since the loading fields have +90 of the last
    % stage as their z=0, and the decelerator has -90 of the first stage as
    % its z=0, the full decelerator length defined as numstages *
    % stagelength gives the needed coordinate shift:
    r.pos(:,3) = r.pos(:,3) - r.decellength;
    
    % Store stage length since its useful for checking the length of the
    % loading fields after we overwrite the decel fields with them:
    r.zstagel = r.f.zstagel;
    
    % Load the mat file, or generate it from a COMSOL .dat file if it
    % doesn't exist yet. We overwrite the deceleration fields, since fields
    % use a lot of RAM and aren't saved at the end anyway. 
    if ~r.reloadfields && exist(['../Traploading/' r.loadname '.mat'],'file');
        r.f = load(['../Traploading/' r.loadname '.mat']);
    elseif exist(['../Traploading/' r.loadname '.dat'],'file');
        processloadingfields(r);
        r.f = load(['../Traploading/' r.loadname '.mat']);
    else
        error(['File ''slowANDtrap/Traploading/' r.loadname '.dat'' not found']);
    end
    
    % Since the loading fields include one decelerator stage, and have
    % their zero at the end of this, we can get the total length of the
    % loading fields by adding the maximum z coordinate to the stage
    % length. Note that zstagel is a full pair of stages (360 degrees) so
    % we divide by two.
    r.loadlength = r.f.zmaxload + r.zstagel/2;
    r.offsetz = r.loadzero;
    
    r.xtrascale = r.voltagescaling;
    
end
