function opt = SE_Laplace_params(varargin)
%SE_Laplace_params  Automatic parameter selection for the Laplace solver
%
%   opt = SE_Laplace_params(box, f, ...)
%       Creates a parameter selection struct `opt'.
%
%   opt = SE_Laplace_params(opt, ...)
%       Extends an existing `opt' struct.
%
%   Required input parameters (first syntax):
%       :param box: size of periodic cell [L1,L2,L3]
%       :param f: source charges (N×1)
%
%   In the second syntax, the given `opt' struct is assumed to
%   already contain the fields 'box' and 'source_quantity' (see below).
%
%   Additional parameters may optionally be given as name-value pairs:
%       - 'periodicity': periodicity of simulation (integer between
%                        0 and 3, inclusive), ignored for Laplace
%       - 'AbsTol': absolute RMS error tolerance
%       - 'RelTol': relative RMS error tolerance
%       - 'rc': cutoff radius for the real-space sum
%       - 'xi': Ewald decomposition parameter
%       - 'grid_res': number of grid subintervals per unit length,
%                     for the Fourier-space sum
%       - 'window': window function for Fourier-space sum (one of
%                   'gaussian', 'expsemicirc', 'kaiser_exact' and
%                   'kaiser_poly'); default is 'kaiser_poly';
%                   may also be 'none' if window-related parameters
%                   should not be computed
%       - 'window_P': window function support size, measured in
%                     number of grid subintervals
%       - 'polynomial_degree': degree of polynomial window function
%       - 'fix_policy': policy controlling whether grid_res and
%                       window_P should be further adjusted; one
%                       of 'never', 'always', 'simple' and
%                       'indicator'; default is 'simple'
%       - 'fix_applied': flag to indicate that an adjustment has
%                        been made
%
%   The effective absolute RMS error tolerance is set to
%
%       EffAbsTol = min(AbsTol, A*RelTol),
%
%   where A is an estimate of the rms of the Fourier-space part
%   of the potential. If neither AbsTol nor RelTol are given, the
%   default is an absolute tolerance of 1e-8, i.e. EffAbsTol = 1e-8.
%   If RelTol is specified, xi should also be specified.
%
%   Returns a structure `opt' which contains some of the parameters
%   mentioned above as fields, namely 'box', 'source_quantity',
%   'periodicity, 'AbsTol', 'RelTol', 'EffAbsTol', 'rc', 'xi',
%   'grid_res', 'window', 'window_P', 'polynomial_degree',
%   'fix_policy' and 'fix_applied'.
%   Note that source_quantity = sum(f.^2), i.e. not all charges are
%   stored in `opt' but only the sum of their squares.
%
%   Parameters which cannot be determined are not included in
%   `opt'. In particular:
%       - At least one of 'rc', 'xi' and 'grid_res' must be given
%         in order to compute the remaining ones using the
%         Kolafa-Perram error estimates. Otherwise these
%         parameters are not included.
%       - The parameter 'window_P' is computed from estimates of
%         the approximation error, using the default window shape
%         factors. However, if 'window' is 'none', it is not
%         computed.

% Determine which calling syntax is used
if nargin < 1
  error('Not enough input arguments.');
end
if isstruct(varargin{1}) % second calling syntax
  opt = varargin{1};
  box = opt.box;
  source_quantity = opt.source_quantity;
  idx_start = 2;
elseif nargin < 2
  error('Not enough input arguments.');
else % first calling syntax
  box = varargin{1};
  source_quantity = sum(varargin{2}(:).^2); % sum of squared charges
  opt = struct('box', box, 'source_quantity', source_quantity);
  idx_start = 3;
end

% Go through name-value arguments
for j=idx_start:numel(varargin)
  arg = varargin{j};
  if mod(j-idx_start,2) == 0 % name
    name = arg;
    if ~anystr(name, {'periodicity', 'AbsTol', 'RelTol', 'rc', 'xi', 'grid_res', ...
                      'window', 'window_P', 'polynomial_degree', 'fix_policy', ...
                      'fix_applied'})
      error('Invalid parameter name: %s', name);
    end
  else % value
    opt.(name) = arg;
  end
end

AbsTol_given = isfield(opt, 'AbsTol') && ~isempty(opt.AbsTol);
RelTol_given = isfield(opt, 'RelTol') && ~isempty(opt.RelTol);
rc_given = isfield(opt, 'rc') && ~isempty(opt.rc);
xi_given = isfield(opt, 'xi') && ~isempty(opt.xi);
grid_res_given = isfield(opt, 'grid_res') && ~isempty(opt.grid_res);
if ~isfield(opt, 'window') || isempty(opt.window), opt.window = 'kaiser_poly'; end
window_given = ~strcmp(opt.window, 'none');
if ~window_given, opt = rmfield(opt, 'window'); end
P_given = isfield(opt, 'window_P') && ~isempty(opt.window_P);
deg_given = isfield(opt, 'polynomial_degree') && ~isempty(opt.polynomial_degree);
if ~isfield(opt, 'fix_policy') || isempty(opt.fix_policy), opt.fix_policy = 'simple'; end
fix_applied = isfield(opt, 'fix_applied') && opt.fix_applied;

comp_rc_xi_grid_res = (rc_given || xi_given || grid_res_given) && ~(rc_given && xi_given && grid_res_given);
comp_P = window_given && ~P_given;
comp_deg = window_given && ~deg_given;

V = box(1)*box(2)*box(3); % volume of periodic box
Lmin = min(box);
Vmin = Lmin^3;

% Determine tolerance
L = min(box); % TODO: temporary assumption, check this for non-cubic boxes!
co = [12.18, 1.815e-2, 1.080e-4]; CC = 0.92;
xiLfun = @(t) exp(-co(1)./t.^2) .* (1 + co(2)*t + co(3)*t.^2);
opt.A_fun = @(Q, xi, L) CC * sqrt(Q) .* xiLfun(xi.*L) ./ L;
if isfield(opt, 'EffAbsTol'), opt = rmfield(opt, 'EffAbsTol'); end
if AbsTol_given && RelTol_given
  if ~xi_given
    error('xi must be given if RelTol is given');
  end
  A = opt.A_fun(source_quantity, opt.xi, L);
  EffAbsTol = min(opt.AbsTol, A*opt.RelTol);
  opt.EffAbsTol = EffAbsTol;
elseif RelTol_given
  if xi_given
    A = opt.A_fun(source_quantity, opt.xi, L);
    EffAbsTol = A*opt.RelTol;
    opt.EffAbsTol = EffAbsTol;
  elseif comp_rc_xi_grid_res
    error('xi must be given if RelTol is given');
  else
    % Computing window_P from RelTol does in fact not require xi
    A = 1;
    EffAbsTol = A*opt.RelTol; % don't store this in opt since it's not "correct"
  end
elseif AbsTol_given
  EffAbsTol = opt.AbsTol;
  opt.EffAbsTol = EffAbsTol;
else
  EffAbsTol = 1e-8; % default tolerance
  opt.EffAbsTol = EffAbsTol;
end

if comp_rc_xi_grid_res
  Ffactor = compute_FoverVolTol2(EffAbsTol, source_quantity, V);
  Ffactor_min = compute_FoverVolTol2(EffAbsTol, source_quantity, Vmin);
  % It has been seen that Ffactor_min works better for the
  % Fourier-space relations in the case of a non-cubic box.
  % For the real-space relations, we always use Ffactor.
end

% Compute xi
if comp_rc_xi_grid_res && ~xi_given
  if rc_given
    xi_rc = xi_from_rc(opt.rc, Ffactor);
    opt.xi = xi_rc;
  end
  if grid_res_given
    xi_grid_res = xi_from_grid_res(opt.grid_res, Ffactor_min);
    opt.xi = xi_grid_res;
  end
  if rc_given && grid_res_given
    opt.xi = (xi_rc + xi_grid_res)/2;
  end
  opt.xi = check_positive(opt.xi, 'xi');
end

% Compute rc
if comp_rc_xi_grid_res && ~rc_given
  opt.rc = rc_from_xi(opt.xi, Ffactor);
  opt.rc = check_positive(opt.rc, 'rc');
end

% Compute grid_res
if comp_rc_xi_grid_res && ~grid_res_given
  opt.grid_res = grid_res_from_xi(opt.xi, Ffactor_min);
  opt.grid_res = check_positive(opt.grid_res, 'grid_res');
end

% Compute window_P
has_xi = xi_given || comp_rc_xi_grid_res;
has_A = exist('A', 'var');
if comp_P
  if ~has_A
    if ~has_xi
      error('Cannot compute window_P from AbsTol without xi');
    end
    A = opt.A_fun(source_quantity, opt.xi, L);
  end
  if strcmp(opt.window, 'gaussian')
    % For the Gaussian, it's not clear why there should be a
    % square root here, but it works really well in practice.
    shape_factor = sqrt(0.91)*pi/2; % default value for Gaussian window
    extra_factor = 2; % from the estimates
  else
    shape_factor = 2.5; % default value for Kaiser windows
    extra_factor = 10; % from the estimates
  end
  opt.window_P = -log(EffAbsTol/A/extra_factor) / shape_factor;
  opt.window_P = check_positive(opt.window_P, 'window_P');
end

% Compute polynomial degree
if comp_deg && strcmp(opt.window, 'kaiser_poly')
  opt.polynomial_degree = max(min(ceil(opt.window_P/2)+2, 10), 1);
end

% Adjust grid_res and window_P if necessary
has_grid_res = grid_res_given || comp_rc_xi_grid_res;
has_P = P_given || comp_P;
has_deg = (deg_given || comp_deg) && strcmp(opt.window, 'kaiser_poly');
deg_given = isfield(opt, 'polynomial_degree') && ~isempty(opt.polynomial_degree);
if ~fix_applied && has_grid_res && has_xi && has_P && ~strcmp(opt.fix_policy, 'never')
  if strcmp(opt.fix_policy, 'always')
    I = true(size(opt.window_P));
  elseif strcmp(opt.fix_policy, 'simple')
    I = true(size(opt.window_P));
  elseif strcmp(opt.fix_policy, 'indicator')
    % Check where indicator (window_P > threshold)
    L = min(box); % TODO: temporary assumption, check this for non-cubic boxes!
    Q = source_quantity;
    gr = opt.grid_res;
    X = gr/opt.xi; % grid_res/xi
    tmp = -pi^2*X.^2/4 + 0.5*log(8*Q./(pi^4*X.^2.*gr.*L.^3));
    if strcmp(opt.window, 'gaussian')
      coeff = [-0.91, 1.15];
    else
      coeff = [-0.67, 1.25];
    end
    threshold = coeff(1)*tmp/log(10) + coeff(2);
    I = opt.window_P > threshold;
  end
  % Where the indicator signals, set P=P+4 and grid_res=grid_res*1.05
  % since the values provided by the estimates are too small in
  % that case
  opt.window_P(I) = opt.window_P(I) + 4;
  opt.grid_res(I) = opt.grid_res(I) * 1.05;
  if has_deg
    deg_rule = max(min(ceil(opt.window_P(I)/2)+2, 10), 1);
    opt.polynomial_degree(I) = max(opt.polynomial_degree(I), deg_rule);
  end
  opt.fix_applied = any(I(:));
elseif ~fix_applied
  opt.fix_applied = false;
end

function fac = compute_FoverVolTol2(tol, F, V)
  fac = F./(V*tol.^2);

function rc = rc_from_xi(xi, Ffactor)
  c = Ffactor/xi;
  rc = (sqrt(3)/(2*xi)) * sqrt(lambertw( (4/3)*c.^(2/3) ));

function xi = xi_from_rc(rc, Ffactor)
  c = Ffactor*rc;
  xi = (1/rc) * sqrt(lambertw( sqrt(c) ));

function grid_res = grid_res_from_xi(xi, Ffactor)
  c = Ffactor/(pi*xi);
  grid_res = (sqrt(3)*xi/pi) * sqrt(lambertw( (4/3)*c.^(2/3) ));

function xi = xi_from_grid_res(grid_res, Ffactor)
  c = Ffactor/(pi^2*grid_res);
  xi = (pi*grid_res/sqrt(2)) ./ sqrt(lambertw( 4*c ));

function found = anystr(str, strcell)
  for j=1:numel(strcell)
    if strcmp(str, strcell{j})
      found = true;
      return
    end
  end
  found = false;

function val = check_positive(val, str)
  if imag(val)/real(val) < 1e-16
    val = real(val);
  end
  assert(isreal(val) && val > 0, sprintf('Fail: %s did not become positive', str));
