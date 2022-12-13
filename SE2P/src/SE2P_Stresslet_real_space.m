function out = SE2P_Stresslet_real_space(x, q, n, opt)
%SE2P_Stresslet_real_space  Compute real-space part of the
%2-periodic Ewald sum for the stresslet potential.
%
%   out = SE2P_Stresslet_real_space(x, q, n, opt)
%
%   Computes the stresslet velocity field from N point sources.
%   Periodicity is assumed in the first two directions (x and y),
%   while the third direction (z) is assumed free.
%
%   Input parameters:
%       :param x: source locations (N×3)
%       :param q: source strengths (N×3)
%       :param n: source normals (N×3)
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
if strcmp(opt.eval_idx, 'default') || numel(opt.eval_idx) > 0
  ts = tic;
  out.pot = SE2P_Stresslet_real_cut_cell_mex(x, q, n, opt.rc, opt.xi, opt.box);
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
  out.extpot = SE2P_Stresslet_direct_real_cut_ext_mex(x, q, n, opt);
  out.walltime.extpot = toc(ts);
end

out.walltime.total = sum(struct2array(out.walltime));

end
