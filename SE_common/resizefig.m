function resizefig(fig, width, height)
  if nargin < 2
    error('resizefig: Both width and height must be given');
  elseif nargin < 3
    handle = gcf;
    % No handle given, only width and height
    height = width;
    width = fig;
  elseif ishandle(fig)
    handle = fig;
  else
    handle = sfigure(fig);
  end
  pos = get(handle, 'Position');
  set(handle, 'Position', [pos(1), pos(2), width, height]);
end
