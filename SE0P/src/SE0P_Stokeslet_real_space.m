function out = SE0P_Stokeslet_real_space(x, f, opt)
%SE0P_Stokeslet_real_space  Compute real-space part of the
%0-periodic (free-space) Ewald sum for the stokeslet potential.
%
%   out = SE0P_Stokeslet_real_space(x, f, opt)
%
%   Computes the Stokesian velocity field from N point forces.
%
%   Input parameters:
%       :param x: source locations (N×3)
%       :param f: source forces (N×3)
%       :param opt: option struct, see below
%
%   Valid options for the opt struct are as follows.
%
%   Options controlling where to evaluate the potential:
%       :param opt.eval_idx: indices of source points to evaluate
%           at, corresponding to rows in x (default: all)
%       :param opt.eval_ext_x: external (non-source) points to
%           evaluate at (N_ext×3) (default: none)
%
%   Options controlling the Ewald decomposition:
%       :param opt.box: size of periodic cell [L1,L2,L3] (required)
%       :param opt.xi: Ewald decomposition parameter (required)
%       :param opt.rc: cutoff radius (required)
%
%   Return values:
%       :returns: **out.pot** (array) -- potential evaluated at source
%           points selected by opt.eval_idx
%       :returns: **out.extpot** (array) -- potential evaluated at
%           external (non-source) points selected by opt.eval_ext_x
%       :returns: **out.walltime** -- struct containing timings
%       :returns: **out.opt** -- final option struct, modified by
%           the program; actual parameters used by the algorithm

assert(isfield(opt, 'box'), 'Periodic cell size "box" must be given in opt struct');
assert(isfield(opt, 'xi'), 'Ewald decomposition parameter "xi" must be given in opt struct');
assert(isfield(opt, 'rc'), 'Cutoff radius "rc" must be given in opt struct');

if ~isfield(opt,'eval_idx'), opt.eval_idx = 'default'; end
if ~isfield(opt,'eval_ext_x'), opt.eval_ext_x = []; end

% Prepare output struct
out.pot = [];
out.extpot = [];
out.walltime = struct('pot',0, 'extpot',0);
out.opt = opt;

% Evaluate
% FIXME: Would be better if the MEX file evaluated only in
% opt.eval_idx instead of having to slice afterwards.
%   Equivalent MATLAB code found at the bottom (in comment)
if strcmp(opt.eval_idx, 'default') || numel(opt.eval_idx) > 0
  ts = tic;
  out.pot = SE0P_Stokeslet_real_cut_cell_mex(x, f, opt.rc, opt.xi, opt.box);
  if ~strcmp(opt.eval_idx, 'default')
    out.pot = out.pot(opt.eval_idx,:);
  end
  out.walltime.pot = toc(ts);
end
% FIXME: Would be better to have fast (cell-based) code for
% external evaluation, using direct summation for now.
if numel(opt.eval_ext_x) > 0
  ts = tic;
  opt.layers = 1;
  out.extpot = SE0P_Stokeslet_direct_real_cut_ext_mex(x, f, opt);
  out.walltime.extpot = toc(ts);
end

out.walltime.total = sum(struct2array(out.walltime));

end

% Equivalent MATLAB code (not external):
%   xi = opt.xi;
%   rc = opt.rc;
%   N = size(x, 1);
%   u = zeros(N, 3);
%   [idx, d] = rangesearch(x, x, rc);
%   for target=opt.eval_idx
%     nnb_idx = idx{target}(2:end);
%     r = d{target}(2:end)';
%     r2 = r.^2;
%     rvec = bsxfun(@minus, x(target,:), x(nnb_idx,:));
%     fsrc = f(nnb_idx,:);
%     xiexp = xi * exp(-xi^2*r2);
%     c1 = 2*( xiexp ./ (sqrt(pi)*r2) + erfc(xi*r) ./ (2*r.^3) );
%     c2 = -4/sqrt(pi)*xiexp;
%     rdotf = sum(rvec.*fsrc, 2);
%     u(target,:) = sum(bsxfun(@times, c1.*r2 + c2, fsrc) + ...
%                       bsxfun(@times, c1.*rdotf, rvec), ...
%                       1);
%   end
