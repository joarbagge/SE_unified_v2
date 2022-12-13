function opt = SE3P_check_options(opt)
%SE3P_check_options  Check option struct for Spectral Ewald.
%
%   opt = SE3P_check_options(opt)
%
%   For documentation of options, see "help SE3P_base_fourier_space"

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

% Make sure window_P is an even integer
opt.window_P = 2*ceil(opt.window_P/2);
if strcmp(opt.window, 'kaiser_poly')
  opt.window_P = min(opt.window_P, 20); % maximum P for kaiser_poly
else
  opt.window_P = min(opt.window_P, 32); % maximum P for other windows
end

% Window shape parameter
if strcmp(opt.window, 'gaussian')
  opt = set_default_option(opt, 'window_shape_factor', 0.91*pi/2); % gaussian
else
  opt = set_default_option(opt, 'window_shape_factor', 2.5); % other windows
end
if strcmp(opt.window, 'kaiser_poly') && opt.window_shape_factor ~= 2.5
  warning('SE3P:PolynomialShapeFactor', ...
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

% Adjust grid resolution
old_grid_res = opt.grid_res;
[opt.grid, opt.grid_res] = SE_compute_grid_size(opt.box, opt.grid_res, opt.base_factor, 3);
if opt.grid_res / old_grid_res > 2
  warning('SE3P:ChangedGridRes', 'grid_res increased by more than a factor 2 (from %g to %g)', ...
          old_grid_res, opt.grid_res);
end
opt.h = 1/opt.grid_res; % grid step size
% Check that h is the same in all directions
Eps = eps(opt.h);
Thres = 4*Eps;
diff1 = abs(opt.h - opt.box(1)/opt.grid(1));
assert(diff1 <= Thres, 'Step size mismatch: abs(h-h1) = %g (> %g)', diff1, Thres);
diff2 = abs(opt.h - opt.box(2)/opt.grid(2));
assert(diff2 <= Thres, 'Step size mismatch: abs(h-h2) = %g (> %g)', diff2, Thres);
diff3 = abs(opt.h - opt.box(3)/opt.grid(3));
assert(diff3 <= Thres, 'Step size mismatch: abs(h-h3) = %g (> %g)', diff3, Thres);

opt.window_halfwidth = opt.h * opt.window_P/2;

end

function opt = set_default_option(opt, field, def_val)
  if ~isfield(opt, field), opt.(field) = def_val; end
end
