function opt = SE2P_check_options(opt)
%SE2P_check_options  Check option struct for Spectral Ewald.
%
%   opt = SE2P_check_options(opt)
%
%   For documentation of options, see "help SE2P_base_fourier_space"

% Check that required options are present
assert(isfield(opt, 'box'), 'Periodic cell size "box" must be given in opt struct');
assert(isfield(opt, 'grid_res'), 'Grid resolution "grid_res" must be given in opt struct');
assert(isfield(opt, 'xi'), 'Ewald decomposition parameter "xi" must be given in opt struct');
assert(isfield(opt, 'window_P'), 'Window support size "window_P" must be given in opt struct');
assert(isfield(opt, 'kernel'), 'Kernel name "kernel" must be given in opt struct');

% Default values
opt = set_default_option(opt, 'return_pot', true);
opt = set_default_option(opt, 'return_gridval', false);
opt = set_default_option(opt, 'return_fourierval', false);
if strcmp(opt.kernel, 'laplace')
  opt = set_default_option(opt, 'return_potgrad', false);
  opt = set_default_option(opt, 'return_gridvalgrad', false);
  opt = set_default_option(opt, 'return_fouriervalgrad', false);
  opt = set_default_option(opt, 'use_ik_diff', false);
  if ~opt.use_ik_diff
    opt.return_gridvalgrad = false;
    opt.return_fouriervalgrad = false;
  end
end
opt = set_default_option(opt, 'eval_idx', 'default');
opt = set_default_option(opt, 'eval_ext_x', []);
opt = set_default_option(opt, 'base_factor', 4);
opt = set_default_option(opt, 'window', 'kaiser_poly');
opt = set_default_option(opt, 'window_scaling_power', 2);
opt = set_default_option(opt, 'use_fast_gridding', true);
opt = set_default_option(opt, 'EffAbsTol', 1e-16);
opt = set_default_option(opt, 'upsampling_zero', 2);
opt = set_default_option(opt, 'upsampling_global', 1);

% Make sure window_P is an even integer
opt.window_P = 2*ceil(opt.window_P/2);
if strcmp(opt.window, 'kaiser_poly')
  opt.window_P = min(opt.window_P, 20); % maximum P for kaiser_poly
else
  opt.window_P = min(opt.window_P, 32); % maximum P for other windows
end

% Default value for grid_add3
if strcmp(opt.window, 'gaussian')
  opt = set_default_option(opt, 'grid_add3', 0.5*opt.window_P);
else
  opt = set_default_option(opt, 'grid_add3', 1.4*opt.window_P);
end

% Window shape parameter
if strcmp(opt.window, 'gaussian')
  opt = set_default_option(opt, 'window_shape_factor', 0.91*pi/2); % gaussian
else
  opt = set_default_option(opt, 'window_shape_factor', 2.5); % other windows
end
if strcmp(opt.window, 'kaiser_poly') && opt.window_shape_factor ~= 2.5
  warning('SE2P:PolynomialShapeFactor', ...
    'Shape factor must be 2.5 for kaiser_poly, setting it to that now');
  opt.window_shape_factor = 2.5;
end
opt.window_shape = opt.window_shape_factor * opt.window_P;

% Window-specific options
if strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
  opt.kaiser_scaling = 1/besseli(0,opt.window_shape);
end
if strcmp(opt.window, 'kaiser_poly')
  poly_deg = max(min(opt.window_P/2 + 2, 10), 1);
  opt = set_default_option(opt, 'polynomial_degree', poly_deg);
end

% Adjust grid resolution (basic grid)
old_grid_res = opt.grid_res;
[opt.grid, opt.grid_res] = SE_compute_grid_size(opt.box, opt.grid_res, opt.base_factor, 2);
if opt.grid_res / old_grid_res > 2
  warning('SE2P:ChangedGridRes', 'grid_res increased by more than a factor 2 (from %g to %g)', ...
          old_grid_res, opt.grid_res);
end
opt.h = 1/opt.grid_res; % grid step size
opt.window_halfwidth = opt.h * opt.window_P/2;
% Extend grid in free direction (z direction)
if ischar(opt.grid_add3) && strcmp(opt.grid_add3, 'cover_remainder')
  deltaM = cover_remainder_method(opt);
else % default method
  deltaM = opt.window_P + opt.grid_add3;
end
% TODO/FIXME: Maybe deltaM should be adjusted by multiplying with opt.grid_res/old_grid_res?
opt.grid3 = opt.grid(3) + deltaM;
opt.grid3 = opt.base_factor * ceil(opt.grid3 / opt.base_factor); % round up
opt.box3 = opt.h * opt.grid3; % adjust box side length
dbox3_half = (opt.box3 - opt.box(3))/2;
LextOverL = opt.box3 / opt.box(3);
% Store values in standard vectors as well
opt.grid(3) = opt.grid3;
opt.box(3) = opt.box3;
% Check that h is the same in all directions
Eps = eps(opt.h);
Thres = 4*Eps;
diff1 = abs(opt.h - opt.box(1)/opt.grid(1));
assert(diff1 <= Thres, 'Step size mismatch: abs(h-h1) = %g (> %g)', diff1, Thres);
diff2 = abs(opt.h - opt.box(2)/opt.grid(2));
assert(diff2 <= Thres, 'Step size mismatch: abs(h-h2) = %g (> %g)', diff2, Thres);
diff3 = abs(opt.h - opt.box(3)/opt.grid(3));
assert(diff3 <= Thres, 'Step size mismatch: abs(h-h3) = %g (> %g)', diff3, Thres);

wbox = [0 opt.box(1); 0 opt.box(2); -dbox3_half opt.box3+dbox3_half];
opt.free_offset = wbox(3,1);
opt.greens_truncation_R = opt.box3;

% Adaptive upsampling factors
if ~isfield(opt,'A_fun')
  % TODO/FIXME: This default is for Laplace, would perhaps be
  % better to have different default for different kernels.
  % Also, should we estimate the total potential rather than the
  % Fourier-space part?
  co = [12.62, 0.8909, 0.01411, 4.315e-5];
  xiLfun = @(x) exp(-co(1)./x.^2) .* (co(2) + co(3)*x + co(4)*x.^2);
  opt.A_fun = @(Q, xi, L) sqrt(Q) * xiLfun(xi*L) / L;
end
Lperiodic = max(opt.box(1), opt.box(2));
A = opt.A_fun(opt.source_quantity, opt.xi, Lperiodic);
if ~isfield(opt,'upsampling_local') % upsampling factor for local pad
  opt.upsampling_local = (1 - 1/(2*pi) * log(2*opt.EffAbsTol/A)) / LextOverL;
  % TODO/FIXME: This default value may NOT work for a noncubic periodic box!
end
if ~isfield(opt,'local_modes_size') % size of local pad
  opt.local_modes_size = -1/(2*pi) * log(2*opt.EffAbsTol/A) / (LextOverL-1) - 1;
  % TODO/FIXME: This default value may NOT work for a noncubic periodic box!
end
opt.local_modes_size = ceil(max(opt.local_modes_size,0));
% Upsampling factors cannot be below 1
opt.upsampling_zero = max(opt.upsampling_zero, 1);
opt.upsampling_local = max(opt.upsampling_local, 1);
opt.upsampling_global = max(opt.upsampling_global, 1);
% Upsampled grids must be multiple of base_factor
s0 = opt.upsampling_zero;
sl = opt.upsampling_local;
sg = opt.upsampling_global;
grid_s0 = opt.base_factor * ceil(opt.grid3*s0 / opt.base_factor);
grid_sl = opt.base_factor * ceil(opt.grid3*sl / opt.base_factor);
grid_sg = opt.base_factor * ceil(opt.grid3*sg / opt.base_factor);
opt.actual_upsampling_zero = grid_s0 / opt.grid3;
opt.actual_upsampling_local = grid_sl / opt.grid3;
opt.actual_upsampling_global = grid_sg / opt.grid3;
% TODO/FIXME: At the moment opt.upsampling_local and opt.local_modes_size
% apply to both x and y directions, could be more flexible and have
% separate variables.

% Defining the local modes to be upsampled ("local pad")
if opt.actual_upsampling_local ~= opt.actual_upsampling_global
  n = opt.local_modes_size;
  Mx = opt.grid(1);
  My = opt.grid(2);
  % The local pad consists of k = -n:n, and the grid is -M/2:(M/2-1)
  % [if M is even] or -(M-1)/2:(M-1)/2 [if M is odd]. Thus, n can
  % at most be M/2-1 [if M is even] or (M-1)/2 [if M is odd].
  % These are both captured by floor((M-1)/2).
  nx = min(floor((Mx-1)/2), n);
  ny = min(floor((My-1)/2), n);
  % The first index is the zero mode, which should not be included.
  % However, for simplicity we include it here and then overwrite
  % it whenever needed.
  opt.local_modes1 = [1:nx+1 Mx-nx+1:Mx];
  opt.local_modes2 = [1:ny+1 My-ny+1:My];
else
  opt.local_modes1 = 1;
  opt.local_modes2 = 1;
end

end

function opt = set_default_option(opt, field, def_val)
  if ~isfield(opt, field), opt.(field) = def_val; end
end

function deltaM = cover_remainder_method(opt)
% Method to determine deltaM from legacy SE0P code
  deltaM_grid = opt.window_P; % to cover the gridding window
  deltaM_rem = deltaM_grid; % to cover the "remainder Gaussian"
  if strcmp(opt.window, 'gaussian')
    w = opt.window_halfwidth;
    a = opt.window_shape;
    eta = 2*(w*opt.xi)^2 / a;
    if eta < 1
      deltaM_rem = sqrt(4*(1-eta)*a)/opt.xi / opt.h;
    end
  else
    deltaM_rem = sqrt(8*opt.window_shape_factor)/opt.xi / opt.h;
    % TODO/FIXME: LEGACY Stokeslet 2P code had:
    %deltaM_rem = sqrt(4/3*opt.window_shape_factor)/opt.xi / opt.h;
  end
  deltaM = max(deltaM_grid, deltaM_rem);
  if isfield(opt, 'no_extra_support') && opt.no_extra_support == true % undocumented
    deltaM = deltaM_grid;
  end
end
