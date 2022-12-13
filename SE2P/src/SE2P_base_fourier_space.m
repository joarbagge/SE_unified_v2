function out = SE2P_base_fourier_space(x, f, opt)
%SE2P_base_fourier_space  Compute Fourier-space part of a
%2-periodic Ewald sum, using the Spectral Ewald method.
%
%   out = SE2P_base_fourier_space(x, f, opt)
%
%   Input parameters:
%       :param x: source locations (N×3)
%       :param f: source strengths (N×?×...)
%       :param opt: option struct, see below
%
%   Periodicity is assumed in the first two directions (x and y),
%   while the third direction (z) is assumed free.
%
%   The shape of the source strength input depends on the kernel,
%   e.g. N×1 for Laplace, N×3 for the stokeslet and rotlet,
%   and N×3×3 for the stresslet. Here, N is the number of sources.
%
%   Valid options for the opt struct are as follows.
%
%   Options controlling the output:
%       :param opt.return_pot: return potential (default: true)
%       :param opt.return_potgrad: return potential gradient
%           (default: false)
%       :param opt.return_gridval: return values on uniform grid
%           (default: false)
%       :param opt.return_gridvalgrad: return gradient values on
%           uniform grid (default: false)
%       :param opt.return_fourierval: return Fourier coefficients
%           (after scaling step) (default: false)
%       :param opt.return_fouriervalgrad: return Fourier coefficients
%           of gradient (after scaling step) (default: false)
%
%   The gradient options ("return_*grad") are available only for
%   the Laplace kernel. Furthermore, gridvalgrad and fouriervalgrad
%   will be empty unless opt.use_ik_diff is true.
%
%   Options controlling where to evaluate the potential:
%       :param opt.eval_idx: indices of source points to evaluate
%           at, corresponding to rows in x (default: all)
%       :param opt.eval_ext_x: external (non-source) points to
%           evaluate at (N_ext×3) (default: none)
%
%   Options controlling the Spectral Ewald algorithm:
%       :param opt.box: size of periodic cell [L1,L2,L3] (required)
%       :param opt.grid_res: number of grid subintervals per unit
%           length, i.e. 1/(grid step size) (required)
%       :param opt.base_factor: integer that the final grid size,
%           i.e. box * grid_res, should be divisible by (default: 4);
%           the grid size will be rounded up
%       :param opt.xi: Ewald decomposition parameter (required)
%       :param opt.window: window function (default: 'kaiser_poly');
%           valid choices are 'gaussian', 'expsemicirc',
%           'kaiser_exact' and 'kaiser_poly'
%       :param opt.window_P: window function support size,
%           measured in number of grid subintervals (required)
%       :param opt.window_shape_factor: shape parameter (alpha or beta)
%           divided by window_P (default: 0.91*pi/2 for 'gaussian',
%           2.5 for all other windows); cannot be varied for the
%           'kaiser_poly' window
%       :param opt.window_scaling_power: exponent used in the
%           scaling step (default: 2)
%       :param opt.use_fast_gridding: use fast gridding and gathering
%           (default: true)
%       :param opt.use_ik_diff: compute gradients using differentiation
%           in Fourier space and three IFFTs; the default is to
%           use analytical differentiation instead (default: false);
%           available only for the Laplace kernel
%
%   Note that opt.grid_res may be adjusted (increased) to keep
%   the grid consistent. Large adjustments are less likely if the
%   numerator and denominator of L1/L2 are small integers.
%
%   Options controlling the adaptive FFT grid:
%       :param opt.grid_add3: extra number to add to the grid
%           size in the free direction, see below; if set to the
%           special value 'cover_remainder', the grid is instead
%           extended to cover the "remainder Gaussian"
%           (default: 0.5*opt.window_P for 'gaussian' window,
%           1.4*opt.window_P for all other windows)
%       :param opt.EffAbsTol: tolerance for default values of
%           opt.upsampling_local and opt.local_modes_size below
%           (default: 1e-16);
%           NOTE: the error estimates for computing the default
%           values may not be good for a noncubic box!
%       :param opt.upsampling_zero: upsampling factor for the
%           zero mode (default: 2)
%       :param opt.upsampling_local: upsampling factor for the
%           "local pad" of wavenumber vectors (default: depends
%           on opt.EffAbsTol)
%       :param opt.upsampling_global: upsampling factor for the
%           rest of the wavenumber vectors (default: 1)
%       :param opt.local_modes_size: defines the size of the "local
%           pad" of wavenumber vectors (default: depends on
%           opt.EffAbsTol)
%
%   The grid size in the free direction (M3) will be set to
%       M3 = opt.grid_res*opt.box(3) + opt.window_P + opt.grid_add3
%   i.e. both opt.window_P and opt.grid_add3 are added to the
%   value provided by opt.grid_res.
%
%   Options specific to the 'kaiser_poly' window:
%       :param opt.polynomial_degree: polynomial degree (default:
%           min(window_P/2 + 2, 10))
%
%   Options which are set automatically and cannot be set manually:
%       :param opt.kernel: name of the kernel; one of 'laplace',
%           'rotlet', 'stokeslet' and 'stresslet'
%
%   Return values:
%       :returns: **out.pot** (array) -- potential evaluated at source
%           points selected by opt.eval_idx (empty if opt.return_pot
%           is false)
%       :returns: **out.potgrad** (array) -- potential gradient at
%           source points selected by opt.eval_idx (empty if
%           opt.return_potgrad is false); only for the Laplace kernel
%       :returns: **out.extpot** (array) -- potential evaluated at
%           external (non-source) points selected by opt.eval_ext_x
%           (empty if opt.return_pot is false)
%       :returns: **out.extpotgrad** (array) -- potential gradient at
%           external points selected by opt.eval_ext_x (empty if
%           opt.return_potgrad is false); only for the Laplace kernel
%       :returns: **out.gridval** (cell array) -- values on the
%           uniform grid, after the IFFT step but before gathering
%           (empty if opt.return_gridval is false)
%       :returns: **out.gridvalgrad** (cell array) -- gradient values
%           on uniform grid (empty if opt.return_gridvalgrad is false
%           or opt.use_ik_diff is false); only for the Laplace kernel
%       :returns: **out.fourierval** (cell array) -- Fourier
%           coefficients, after scaling step but before IFFT (empty if
%           opt.return_fourierval is false)
%       :returns: **out.fouriervalgrad** (cell array) -- Fourier
%           coefficients of gradient (empty if
%           opt.return_fouriervalgrad is false or opt.use_ik_diff is
%           false); only for the Laplace kernel
%       :returns: **out.walltime** -- struct containing timings
%       :returns: **out.opt** -- final option struct, modified by
%           the program; actual parameters used by the algorithm
%
%   Each cell array output has three components, corresponding to
%   the three spatial components of the vector field: U{1}, U{2},
%   U{3}. (For a scalar field the cell array has only a single
%   component.) Note that the physical meaning of the cell array
%   outputs is affected by opt.window_scaling_power.
%
%   NOTE: This function may give unexpected results if there are
%         sources outside of opt.box in the *free* direction. It is
%         up to the caller to ensure that the sources are inside
%         the box.

opt = SE2P_check_options(opt);
is_laplace = strcmp(opt.kernel, 'laplace');

% Prepare output struct
out.pot = [];
out.extpot = [];
out.gridval = {};
out.fourierval = {};
if is_laplace
  out.potgrad = [];
  out.extpotgrad = [];
  out.gridvalgrad = {};
  out.fouriervalgrad = {};
end
out.walltime = struct('pre_gridding',0, ...
                      'pre_fft',0, ...
                      'gridding',0, ...
                      'fft',0, ...
                      'scaling',0, ...
                      'ifft',0, ...
                      'gathering',0);
out.opt = opt;

if ~(opt.return_pot || opt.return_gridval || opt.return_fourierval ...
    || (is_laplace && (opt.return_potgrad || opt.return_gridvalgrad ...
        || opt.return_fouriervalgrad)))
  return;
end

if is_laplace && opt.return_potgrad && ~opt.use_ik_diff ...
    && strcmp(opt.window, 'expsemicirc')
  error(['Gradient calculations not supported with expsemicirc window' ...
         ' and analytical differentiation']);
end

% Dimensionality of input
size_f = size(f);
N = size_f(1);
dim_in = size_f(2:end);

% 0. Precomputation
ts = tic;
gridding = precompute_gridding(x, opt);
out.walltime.pre_gridding = toc(ts);

if strcmp(opt.window, 'expsemicirc')
  % Precompute window Fourier transform
  ts = tic;
  pre_fft = SE2P_window_fft_precomp(opt);
  out.walltime.pre_fft = toc(ts);
else
  pre_fft = [];
end

% 1. Gridding
ts = tic;
H = cell([dim_in, 1]);
for i=1:prod(dim_in)
  H{i} = gridding.grid_fun(f(:,i));
end
out.walltime.gridding = toc(ts);

% 2. Fourier transform
ts = tic;
Hl = cell([dim_in, 1]);
H0 = cell([dim_in, 1]);
for i=1:prod(dim_in)
  [H{i}, Hl{i}, H0{i}] = SE2P_fftnd(H{i}, opt);
end
out.walltime.fft = toc(ts);

% 3. Scaling
ts = tic;
scaling_fun = select_scaling_fun(opt);
[H, Hl, H0] = scaling_fun(H, Hl, H0, opt, pre_fft);
out.walltime.scaling = toc(ts);
dim_out = numel(H);

if opt.return_fourierval
  if is_laplace && dim_out > 1
    out.fourierval = {H{1}};
  else
    out.fourierval = H;
  end
end
if is_laplace && opt.return_fouriervalgrad && dim_out > 1
  out.fouriervalgrad = {H{2}; H{3}; H{4}};
end
if ~(opt.return_pot || opt.return_gridval ...
    || (is_laplace && (opt.return_potgrad || opt.return_gridvalgrad)))
  out.walltime.total = sum(struct2array(out.walltime));
  return;
end

% 4. Inverse Fourier transform
ts = tic;
for i=1:dim_out
  H{i} = real(SE2P_ifftnd(H{i}, Hl{i}, H0{i}, opt)); % should be real to eps accuracy
end
out.walltime.ifft = toc(ts);

if opt.return_gridval
  if is_laplace && dim_out > 1
    out.gridval = {H{1}};
  else
    out.gridval = H;
  end
end
if is_laplace && opt.return_gridvalgrad && dim_out > 1
  out.gridvalgrad = {H{2}; H{3}; H{4}};
end
if ~(opt.return_pot || (is_laplace && opt.return_potgrad))
  out.walltime.total = sum(struct2array(out.walltime));
  return;
end

% 5. Gathering
ts = tic;
if opt.return_pot
  if is_laplace && dim_out > 1
    u = gridding.gath_fun(H{1});
  else
    u = zeros(N, dim_out);
    for i=1:dim_out
      u(:,i) = gridding.gath_fun(H{i});
    end
  end
  if opt.use_fast_gridding && isnumeric(opt.eval_idx), u = u(opt.eval_idx,:); end
  out.pot = u;
end
if is_laplace && opt.return_potgrad
  if opt.use_ik_diff
    du = zeros(N, dim_out-1);
    for i=1:dim_out-1
      du(:,i) = gridding.gath_fun(H{i+1});
    end
  else
    du = gridding.gath_fun_grad(H{1});
  end
  if opt.use_fast_gridding && isnumeric(opt.eval_idx), du = du(opt.eval_idx,:); end
  out.potgrad = du;
end
out.walltime.gathering = toc(ts);

% Evaluate in external (non-source) points
if numel(opt.eval_ext_x) > 0
  N_ext = size(opt.eval_ext_x, 1);
  % New gridding precomputation for these points needed
  ts = tic;
  ext_gridding = precompute_gridding(opt.eval_ext_x, opt, true);
  out.walltime.pre_gridding = out.walltime.pre_gridding + toc(ts);
  % Gathering in the external points
  ts = tic;
  if opt.return_pot
    if is_laplace && dim_out > 1
      u = ext_gridding.gath_fun(H{1});
    else
      u = zeros(N_ext, dim_out);
      for i=1:dim_out
        u(:,i) = ext_gridding.gath_fun(H{i});
      end
    end
    out.extpot = u;
  end
  if is_laplace && opt.return_potgrad
    if opt.use_ik_diff
      du = zeros(N_ext, dim_out-1);
      for i=1:dim_out-1
        du(:,i) = ext_gridding.gath_fun(H{i+1});
      end
    else
      du = ext_gridding.gath_fun_grad(H{1});
    end
    out.extpotgrad = du;
  end
  out.walltime.gathering = out.walltime.gathering + toc(ts);
end

out.walltime.total = sum(struct2array(out.walltime));

end

% ------------------------------------------------------------------------------
function fun = select_scaling_fun(opt)

if strcmp(opt.kernel, 'laplace')
  fun = @SE2P_Laplace_scaling;
elseif strcmp(opt.kernel, 'rotlet')
  fun = @SE2P_Rotlet_scaling;
elseif strcmp(opt.kernel, 'stokeslet')
  fun = @SE2P_Stokeslet_scaling;
elseif strcmp(opt.kernel, 'stresslet')
  fun = @SE2P_Stresslet_scaling;
else
  error('Unsupported kernel: %s', opt.kernel);
end

end

% ------------------------------------------------------------------------------
function pre_grid = precompute_gridding(x, opt, external_eval)

if nargin < 3
  external_eval = false;
end

% Function selection
if strcmp(opt.window, 'gaussian')
  W_fast_grid = @(x,F,opt,S) SE_fg_grid_split_thrd_mex_2p(x,F,opt,...
                                S.zs,S.zx,S.zy,S.zz,S.idx);
  W_fast_gath = @(F,opt,S) SE_fg_int_split_mex_2p(0,F,opt,...
                                S.zs,S.zx,S.zy,S.zz,S.idx);
  % FIXME: The first argument to SE_fg_int_split_mex_2p is unused.
  % It is not possible to specify the evaluation points, so we
  % have to pick out the correct evaluation points at the end.
  W_fast_gath_grad = @(F,opt,S) SE_fg_int_split_force_mex_2p(x,F,opt,...
                                S.zs,S.zx,S.zy,S.zz,...
                                S.zfx,S.zfy,S.zfz,S.idx);
  W_plain_grid = @SE_fg_grid_mex_2p;
  W_plain_gath = @SE_fg_int_mex_2p;
  W_plain_gath_grad = @SE_fg_int_force_mex_2p;
elseif strcmp(opt.window, 'expsemicirc') || strcmp(opt.window, 'kaiser_exact') ...
    || strcmp(opt.window, 'kaiser_poly')
  W_fast_grid = @(x,F,opt,S) SE_fg_grid_split_kaiser_mex_2p(x,F,opt,...
                                S.zx,S.zy,S.zz,S.idx);
  W_fast_gath = @(F,opt,S) SE_fg_int_split_kaiser_mex_2p(0,F,opt,...
                                S.zx,S.zy,S.zz,S.idx);
  % FIXME: The first argument to SE_fg_int_split_kaiser_mex_2p is unused.
  % It is not possible to specify the evaluation points, so we
  % have to pick out the correct evaluation points at the end.
  W_fast_gath_grad = @(F,opt,S) SE_fg_int_split_kaiser_force_mex_2p(x,F,opt,...
                                S.zx,S.zy,S.zz,...
                                S.zfx,S.zfy,S.zfz,S.idx);
  W_plain_grid = @SE_fg_grid_kaiser_mex_2p;
  W_plain_gath = @SE_fg_int_kaiser_mex_2p;
  W_plain_gath_grad = @SE_fg_int_kaiser_force_mex_2p;
else
  error('Unsupported window function: %s', opt.window);
end

% Perform precomputation and store output
if opt.use_fast_gridding
  S = SE2P_gridding_precomp(x, opt);
  grid_fun = @(F) W_fast_grid(x(S.perm,:), F(S.perm), opt, S);
  iperm = @(u) u(S.iperm,:);
  gath_fun = @(F) iperm(W_fast_gath(F, opt, S));
  gath_fun_grad = @(F) iperm(W_fast_gath_grad(F, opt, S));
else
  warning('SE2P:SlowGridding', 'Using slow gridding routines');
  grid_fun = @(F) W_plain_grid(x, F, opt);
  if strcmp(opt.eval_idx, 'default') || external_eval
    gath_fun = @(F) W_plain_gath(x, F, opt);
    gath_fun_grad = @(F) W_plain_gath_grad(x, F, opt);
  else
    gath_fun = @(F) W_plain_gath(x(opt.eval_idx,:), F, opt);
    gath_fun_grad = @(F) W_plain_gath_grad(x(opt.eval_idx,:), F, opt);
  end
end

pre_grid.grid_fun = grid_fun;
pre_grid.gath_fun = gath_fun;
pre_grid.gath_fun_grad = gath_fun_grad;

end
