%%
function scalefig(handle,scsh,fn)
lst = get(handle,'Children');
for i=1:length(lst)
    t = lst(i);
    scaleobj(t,scsh,fn)
end
end

