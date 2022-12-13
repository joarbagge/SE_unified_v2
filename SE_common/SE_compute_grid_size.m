function [grid, grid_res] = SE_compute_grid_size(box, grid_res, base_factor, periodicity)
% SE_compute_grid_size  Compute grid given the box and desired grid_res
%    [grid, grid_res] = SE_compute_grid_size(box, grid_res, base_factor, periodicity)
%
% Input:
% - box: [L1,L2,L3] side lengths of the box
% - grid_res: desired number of subintervals per unit length
% - base_factor: factor that grid should be a multiple of
% - periodicity: number of periodic directions (3,2,1,0)
%
% Output:
% - grid: [M1,M2,M3] actual number of subintervals
% - grid_res: actual number of subintervals per unit length
%
% This code will adjust (increase) grid_res such that grid is divisible by
% base_factor in all periodic directions and grid./box is the same
% in all periodic directions. The new value of grid_res is returned.

if periodicity == 0
  % Nothing to do in 0P case
  grid = grid_res * box;
elseif periodicity == 1
  % Determine number of "ticks" needed
  n = ceil(grid_res * box(1) / base_factor);
  % Compute number of subintervals
  M1 = n*base_factor;
  % Compute new grid_res
  grid_res = M1/box(1);
  % Collect grid in all directions
  grid = [M1, grid_res*box(2), grid_res*box(3)];
elseif periodicity == 2
  % Simplify the fraction L1/L2 = a1/a2
  [a1, a2] = rat(box(1)/box(2));
  % Determine number of "ticks" needed
  n = ceil(grid_res * box(1) / (a1 * base_factor));
  % Compute number of subintervals
  M1 = n*a1*base_factor;
  M2 = n*a2*base_factor;
  % Compute new grid_res
  grid_res = M1/box(1);
  % Collect grid in all directions
  grid = [M1, M2, grid_res*box(3)];
elseif periodicity == 3
  % First simplify the fraction L1/L2 = c1/c2
  [c1, c2] = rat(box(1)/box(2));
  q = box(1)/c1;
  % Then simplify the fraction q/L3 = b/a3
  [b, a3] = rat(q/box(3));
  a1 = c1*b;
  a2 = c2*b;
  % Determine number of "ticks" needed
  n = ceil(grid_res * box(1) / (a1 * base_factor));
  % Compute number of subintervals
  M1 = n*a1*base_factor;
  M2 = n*a2*base_factor;
  M3 = n*a3*base_factor;
  % Compute new grid_res
  grid_res = M1/box(1);
  % Collect grid in all directions
  grid = [M1, M2, M3];
end

end
