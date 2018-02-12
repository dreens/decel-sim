function r = update(r)
%% Data Updating
% This function is called by several other functions in this simulation:
% stage, trapload, and trap.
%
% See the subfunction initvars of init to see where these variables are
% initialized.
%

    if ~r.undo
        
        % Note: time-stepping backwards cannot undo molecules that are lost
        % during the forward step.
        r.lost    = isnan(sum(r.pos,2));
        r.pos     = r.pos(~r.lost,:);
        r.vel     = r.vel(~r.lost,:);
        r.numleft = size(r.pos,1);
        if r.rememberwho
            r.veli = r.veli(~r.lost,:);
            r.posi = r.posi(~r.lost,:);
        end
        
        % If none are left, keep the sim limping along to avoid erroring
        % out.
        if r.numleft == 0
            r.pos = [nan nan nan];
            r.vel = r.pos;
        end
        
        % Update all data variables
        r.zsynch(end+1) = r.pos(1,3) + r.offsetz;
        r.vzsynch(end+1) = r.vel(1,3);
        r.times(end+1) = r.time;
        r.molnum(end+1) = r.numleft;
        r.molnumlaser(end+1) = sum(inlaser(r));
        temp = r.pos; r.pos = [2e-4 0 r.pos(1,3); 0 2e-4 r.pos(1,3)];
        twosix = -r.acc(r)*r.mOH/6.626e-25/1e6/2e-4;
        r.ksynch(end+1,:) = [twosix(1,1) twosix(2,2)];
        r.pos = temp;
        r.pesynch(end+1) = r.f.vf(r.pos(1,:),r.numstage)/6.626e-25;
    else
        % A backward step has occurred, so erase some data.
        r.vzsynch = r.vzsynch(1:end-1);
        r.zsynch = r.zsynch(1:end-1);
        r.times = r.times(1:end-1);
        r.molnum = r.molnum(1:end-1);
        r.molnumlaser = r.molnumlaser(1:end-1);
        r.ksynch = r.ksynch(1:end-1,:);
        r.pesynch = r.pesynch(1:end-1);
        
        r.undo = false;
    end
end