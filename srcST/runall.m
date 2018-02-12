%% Loop through the array of runstructs.
% Takes in an array produced in simulate() by a call to unpacker()
function rsf = runall(rs)

    % Get a time string for file saving
    t = datestr(now,'mmmm-dd-yyyy_HH-MM-SS');


    % Decide whether or not to use parfor.
    if rs(1).parallel 
        if rs(1).subparallel
            r = rs(1);
            r.num = r.num/r.cores;
            rs = repmat(r,1,r.cores);
        end
        parfor i=1:length(rs)
            
            % Seed the random number generator
            if rs(i).doseed
                rng(rs(i).seed);
            else
                rng;
            end
            
            % Some status output
            fprintf('run:%3d/%d\n',i,length(rs))
            
            % Initialize all Variables
            rr = init(rs(i));
            
            % Perform a run, populate a struct of run results (rsf)
            rsf(i) = run(rr);
            
            % Don't resave the fields, they're static.
            rsf(i).f = 0;
        end
        if rs(1).subparallel
            % Now put arrays back together in a single struct
            rfinal = rsf(1);
            vcat = {'vel','pos','veli','posi','veliall','posiall','lost'};
            sigma = {'molnum','molnumlaser'};
            for j=2:r.cores
                for i=1:length(vcat)
                    if isfield(rfinal,vcat{i})
                        f = vcat{i};
                        rfinal.(f) = vertcat(rfinal.(f),rsf(j).(f));
                    end
                end
                for i=1:length(sigma)
                    if isfield(rfinal,sigma{i})
                        f = sigma{i};
                        l1 = length(rfinal.(f));
                        l2 = length(rsf(j).(f));
                        new = zeros(1,max(l1,l2));
                        new(1:l1) = rfinal.(f);
                        new(1:l2) = new(1:l2) + rsf(j).(f);
                        rfinal.(f) =new;
                    end
                end
            end
            rsf = rfinal;
        end
    else
        for i=1:length(rs)
            
            % Same Comments as above
            if rs(i).doseed
                rng(rs(i).seed) ;
            else
                rng;
            end
            fprintf('run:%3d/%d\n ',i,length(rs))
            rr = init(rs(i));
            rr = run(rr);
            rr.f = 0;
            rsf(i) = rr;
            
            % Since we're not in parallel, save after each run
            save(['../autosaves/runstructs_' t '_' rs(1).dname '.mat'],'rsf')
        end
    end

    % Save the struct of runs, just in case the fitsim2data or results
    % functions have petty errors. Wouldn't want to lose everything.
    save(['../autosaves/runstructs_' t '_' rs(1).dname '.mat'],'rsf')
    
    % Make a copy of this file, for future reference.
    system(['cp -r ./ ../autosaves/src_' t '_' rs(1).dname '/']);
    
    % Can't plot directly when running parallel- plotting must be done on
    % the head node.
    if ~rs(1).parallel
        
        % Plot results. This function lives in a separate file and is rewritten
        % and renamed frequently.
        if rs(1).plotdata
            eval([rs(1).plotfunc '(rsf)']);
        end

        % Fit results to experimental data. Also a separate function with many
        % variations.
        if rs(1).fitdata
            eval([rs(1).fitfunc '(rsf)'])
        end
    end

end

