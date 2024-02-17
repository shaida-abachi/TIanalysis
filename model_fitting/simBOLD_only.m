function bold = simBOLD_only(traces, W)

deltabar = mean(traces.delta_it,'all');
xbar = mean(traces.x_it,'all');
ybar = xbar;
sbar = mean(diff(traces.x_it,2),'all');
%sbar = mean(traces.x_it(1,:)); % s1 or s1 is

bold = (W.Ws * sbar) + (W.Wx * xbar) + (W.Wy * ybar) + (W.Wdelta * deltabar);

end
