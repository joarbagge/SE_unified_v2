function h = sfigure(f)
if nargin==0
    h = figure();
elseif ishandle(f)
    set(0, 'CurrentFigure', f);
    h = gcf();
else
    h = figure(f);
end
