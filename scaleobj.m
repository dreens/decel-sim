%%
function scaleobj(h,scsh,fn,varargin)
figure(fn); hold on
if isa(h,'matlab.graphics.chart.primitive.Line')
    plot(h.XData*scsh(3)+scsh(4),h.YData*scsh(1)+scsh(2),'Color',h.Color,'LineStyle',h.LineStyle,'Marker',h.Marker,'MarkerSize',h.MarkerSize)
elseif isa(h,'matlab.graphics.chart.primitive.ErrorBar')
    errorbar(h.XData*scsh(3)+scsh(4),h.YData*scsh(1)+scsh(2),h.YNegativeDelta*2*scsh(1),'Color',h.Color,'LineStyle',h.LineStyle,'Marker',h.Marker,'MarkerSize',h.MarkerSize)
end
if ~isempty(varargin)
    if varargin{1}
        delete(h)
    end
end
end