function a = accn(r)
    ax = r.f.dvdx(r.pos,r.numstage);
    ay = r.f.dvdy(r.pos,r.numstage);
    az = r.f.dvdz(r.pos,r.numstage);
    a = [ax ay az]/r.mOH;
    
end
