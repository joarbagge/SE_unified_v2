function SE_draw_box(box, origin)
%SE_draw_box  Draw periodic cell.
%
%   SE_draw_box(box)
%   SE_draw_box(box, origin)
  if nargin < 2 || isempty(origin)
    origin = [0 0 0];
  end
  L = box;
  xy = [
    0      0      0
    0      L(2)   0
    L(1)   L(2)   0
    L(1)   0      0
    0      0      0];

  xy = bsxfun(@plus, xy, origin);

  hstate = ishold();
  plot3(xy(:,1), xy(:,2), xy(:,3), '-k'), hold on
  plot3(xy(:,1), xy(:,2), xy(:,3)+L(3), '-k')
  for i=1:5
    plot3(xy([i i],1), xy([i i],2), xy([i i],3)+[0 L(3)]', '-k')
  end
  axis equal
  if ~hstate
    hold off
  end
end
