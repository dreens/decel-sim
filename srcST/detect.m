function r = detect(r)
%% Detection Sequences
% We may implement a detection step in the experiment where free flight or
% direct expulsion are used in conjunction with LIF.
%
% In the free flight case, we can't turn off the magnetic fields, so we'll
% need to recompute the trapping fields with E turned to zero in the
% Hamiltonian.
%
% In the direct expulsion case, we'll want to combine the E-fields from a
% deceleration stage with the B-fields in the trap.
%
% Eventually, it would be nice to map out what portion of trapped phase
% space a given detection sequence is sensitive to. It may be that allowing
% free flight before detection enables us to sample the cold portion of the
% trap only for example.

    switch(r.detection)
        case 'in situ' % default- no propagation needed
        case 'slowing' % turn off the fields
            r.xtrascale = 0;
        case 'free flight'
            if ~r.reloadfields && ...
                    exist(['../Traps/' r.trapname '_free.mat'],'file')
                r.f = load(['../Traps/' r.trapname '_free.mat']);
            elseif exist(['../Traps/' r.trapname '.dat'],'file')
                processtrappingfields(r);
                r.f = load(['../Traps/' r.trapname '_free.mat']);
            else
                error(['File ''Traps/' r.trapname '.dat'' not found']);
            end
        case 'expel'
            if ~strcmp(r.trapname,'magneticpintrap')
                error(['Trap Expel Detection only implemented in the ',...
                    'magnetic pin trap'])
            end
            if ~r.reloadfields && exist('../Decels/magneticpin.mat','file')
                r.f = load('../Decels/magneticpin.mat');
            elseif exist('../Decels/magneticpin.dat','file')
                r.decel = 'magneticpin';
                processdecelfields(r);
                r.f = load('../Decels/magneticpin.mat');
            else
                error('File ''Decels/magneticpin.dat'' not found');
            end

            % Now we make sure the stage is even so that molecules are on a
            % hill.
            r.numstage = 2;
    end

    fprintf('Detecting\n')
    
    maxnum = 0;
    
    for i=1:round(r.detecttime/r.detectt)
        for j=1:round(r.detectt/r.smallt)
            r = smallstep(r,r.smallt);
        end
        
        r = update(r);
        
        if r.molnumlaser(end) > maxnum
            maxnum = r.molnumlaser(end);
            r.molsatmaxinlaser = r.pos(inlaser(r),:);
        end
        
        if isnan(r.pos(1,3))
            fprintf('No Molecules Remaining\n');
            break
        end
                
        fprintf('Step %d,\t%dmols\n',i,r.numleft);
    end
