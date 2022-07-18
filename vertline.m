function hl = vertline(ttag, lcolor, lwidth, lstyle)

% hl = vertline(ttag, lcolor, lwidth);
% 
% Adds a vertical line at time <ttag> to the plot. Returns line handle.
% Optional args:
%   lcolor - line color (e.g. 'k')
%   lwidth - line width
%
% Nov 13, 2009 JS

yl = ylim;
hold on;
hl = line([ttag ttag], yl);
if exist('lcolor','var') && ~isempty(lcolor)
	set(hl,'Color', lcolor);
end
if exist('lwidth','var') && ~isempty(lwidth)
	set(hl,'LineWidth', lwidth);
end
if exist('lstyle','var') && ~isempty(lstyle)
	set(hl,'LineStyle', lstyle);
end
hold off;
