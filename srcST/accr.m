function a = accr(r)
    ay = r.f.dvdx(r.pos(:,[2 1 3]),r.numstage);
    ax = r.f.dvdy(r.pos(:,[2 1 3]),r.numstage);
    az = r.f.dvdz(r.pos(:,[2 1 3]),r.numstage);
    a = [ax ay az]/r.mOH;
end
