function opt = SE0P_check_options(opt)
%SE0P_check_options  Check option struct for Spectral Ewald.
%
%   opt = SE0P_check_options(opt)
%
%   For documentation of options, see "help SE0P_base_fourier_space"
%
%   The 0P (free-space) code uses four different grids:
%   - inner/original: the domain containing all sources and target;
%                     given by opt.box and opt.grid
%   - extended: extended grid used for gridding, to cancel wrap
%               effects from the periodic gridding code;
%               given by opt.extended_box and opt.extended_grid
%   - padded: FFT grid with upsampling factor 2, used for
%             aperiodic convolution; given by opt.padded_box and
%             opt.padded_grid
%   - upsampled: FFT grid with upsampling factor opt.upsampling,
%                used for truncated Green's function;
%                given by opt.upsampled_box and opt.upsampled_grid

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
opt = set_default_option(opt, 'upsample_all', false);
if strcmp(opt.kernel, 'stokeslet') || strcmp(opt.kernel, 'stresslet') ...
  || (strcmp(opt.kernel, 'rotlet') && isfield(opt, 'vico_variant') && opt.vico_variant == 2)
  % Default for kernels based on the biharmonic
  opt = set_default_option(opt, 'upsampling_padded', 2); % undocumented
else
  opt = set_default_option(opt, 'upsampling_padded', 2); % undocumented
end
% The following two options are used for precomputation of the biharmonic
opt = set_default_option(opt, 'upsampling_convolution', -8); % undocumented
opt = set_default_option(opt, 'transition_level', 1e-2); % undocumented

% Make sure window_P is an even integer
opt.window_P = 2*ceil(opt.window_P/2);
if strcmp(opt.window, 'kaiser_poly')
  opt.window_P = min(opt.window_P, 20); % maximum P for kaiser_poly
else
  opt.window_P = min(opt.window_P, 32); % maximum P for other windows
end

% Default value for grid_add
if strcmp(opt.kernel, 'laplace')
  if strcmp(opt.window, 'gaussian')
    opt = set_default_option(opt, 'grid_add', 0);
  else
    opt = set_default_option(opt, 'grid_add', 0.3*opt.window_P);
  end
elseif strcmp(opt.kernel, 'stokeslet')
  if strcmp(opt.window, 'gaussian')
    opt = set_default_option(opt, 'grid_add', 0.6*max(opt.window_P,8));
  else
    opt = set_default_option(opt, 'grid_add', 1.2*max(opt.window_P,8));
  end
elseif strcmp(opt.kernel, 'rotlet')
  if strcmp(opt.window, 'gaussian')
    opt = set_default_option(opt, 'grid_add', 0);
  else
    opt = set_default_option(opt, 'grid_add', 0.5*opt.window_P);
  end
elseif strcmp(opt.kernel, 'stresslet')
  if strcmp(opt.window, 'gaussian')
    opt = set_default_option(opt, 'grid_add', 0.9*max(opt.window_P,8));
  else
    opt = set_default_option(opt, 'grid_add', 1.4*max(opt.window_P,8));
  end
end

% Window shape parameter
if strcmp(opt.window, 'gaussian')
  opt = set_default_option(opt, 'window_shape_factor', 0.91*pi/2); % gaussian
else
  opt = set_default_option(opt, 'window_shape_factor', 2.5); % other windows
end
if strcmp(opt.window, 'kaiser_poly') && opt.window_shape_factor ~= 2.5
  warning('SE0P:PolynomialShapeFactor', ...
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

% Inner/original grid -- given by opt.box and other input
[opt.grid, opt.grid_res] = SE_compute_grid_size(opt.box, opt.grid_res, opt.base_factor, 0);
opt.h = 1/opt.grid_res; % grid step size
opt.window_halfwidth = opt.h * opt.window_P/2;

% Extended grid -- used for gridding
if ischar(opt.grid_add) && strcmp(opt.grid_add, 'cover_remainder')
  deltaM = cover_remainder_method(opt);
else % default method
  deltaM = opt.window_P + opt.grid_add;
end
opt.extended_grid = opt.grid + deltaM;
opt.extended_grid = opt.base_factor * ceil(opt.extended_grid / opt.base_factor); % round up
opt.extended_box = opt.h * opt.extended_grid; % adjust box side length

% Upsampling factor
dbox_half = (opt.extended_box - opt.box)/2;
opt.free_offset = -dbox_half;
opt.greens_truncation_R = norm(opt.extended_box);
def_upsampling = 1 + opt.greens_truncation_R/min(opt.extended_box);
def_upsampling = ceil(def_upsampling*10)/10; % round up slightly
opt = set_default_option(opt, 'upsampling', def_upsampling);

opt.upsampling = max(opt.upsampling, 1);
%opt.upsampling_padded = max(opt.upsampling_padded, 1);
%opt.upsampling_convolution = max(opt.upsampling_convolution, 1);
if opt.upsample_all % use opt.upsampling factor also for the padded grid
  opt.upsampling_padded = opt.upsampling;
  opt.transition_level = [];
end

% Padded grid -- used for FFTs
s = opt.upsampling_padded;
if s > 0
  grid_s = opt.base_factor * ceil(s * opt.extended_grid / opt.base_factor); % round up
else
  grid_s = opt.base_factor * ceil((2*opt.extended_grid - s) / opt.base_factor); % round up
end
opt.actual_upsampling_padded = grid_s ./ opt.extended_grid;
opt.padded_grid = grid_s;
opt.padded_box = opt.extended_box .* opt.actual_upsampling_padded;

% Upsampled grid -- used for precomputing the Green's function
s = opt.upsampling;
if ~opt.upsample_all
  % Upsampled grid can never be smaller than padded grid
  s = max(s, opt.actual_upsampling_padded);
end
grid_s = opt.base_factor * ceil(s .* opt.extended_grid / opt.base_factor); % round up
opt.actual_upsampling = grid_s ./ opt.extended_grid;
opt.upsampled_grid = grid_s;
opt.upsampled_box = opt.extended_box .* opt.actual_upsampling;

% Check that h is the same in all directions (not done in LEGACY code)
Eps = eps(opt.h);
Thres = 10*Eps; % 4*Eps
diff1 = abs(opt.h - opt.upsampled_box(1)/opt.upsampled_grid(1));
assert(diff1 <= Thres, 'Step size mismatch: abs(h-h1) = %g (> %g)', diff1, Thres);
diff2 = abs(opt.h - opt.upsampled_box(2)/opt.upsampled_grid(2));
assert(diff2 <= Thres, 'Step size mismatch: abs(h-h2) = %g (> %g)', diff2, Thres);
diff3 = abs(opt.h - opt.upsampled_box(3)/opt.upsampled_grid(3));
assert(diff3 <= Thres, 'Step size mismatch: abs(h-h3) = %g (> %g)', diff3, Thres);

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
    % TODO/FIXME: LEGACY Stokeslet 0P code seems to have had
    %deltaM_rem = sqrt(opt.window_shape_factor)/opt.xi / opt.h;
    % (or maybe it was even opt.window_shape, i.e. multiplied by P)
  end
  deltaM = max(deltaM_grid, deltaM_rem);
  if isfield(opt, 'no_extra_support') && opt.no_extra_support == true % undocumented
    deltaM = deltaM_grid;
  end
end
