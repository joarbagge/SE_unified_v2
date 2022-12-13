function out = SE0P_Stresslet_real_space(x, q, n, opt)
%SE0P_Stresslet_real_space  Compute real-space part of the
%0-periodic (free-space) Ewald sum for the stresslet potential.
%
%   out = SE0P_Stresslet_real_space(x, q, n, opt)
%
%   Computes the stresslet velocity field from N point sources.
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
%   Equivalent MATLAB code found at the bottom (in comment)
if strcmp(opt.eval_idx, 'default') || numel(opt.eval_idx) > 0
  ts = tic;
  out.pot = SE0P_Stresslet_real_cut_cell_mex(x, q, n, opt.rc, opt.xi, opt.box);
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
  out.extpot = SE0P_Stresslet_direct_real_cut_ext_mex(x, q, n, opt);
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
%     qsrc = q(nnb_idx,:);
%     nsrc = n(nnb_idx,:);
%     rdotq = sum(rvec .* qsrc, 2);
%     rdotn = sum(rvec .* nsrc, 2);
%     qdotn = sum(qsrc .* nsrc, 2);
%     c = xi^2*r2;
%     C = -2./r2.^2 .* ( 3.0./r.*erfc(xi*r) + 2.0*xi/sqrt(pi)*(3.0+2.0*c).*exp(-c) );
%     D = 4/sqrt(pi)*xi^3.*exp(-c);
%     u(target,:) = sum(bsxfun(@times, C.*rdotn.*rdotq + D.*qdotn, rvec) + ...
%                       bsxfun(@times, D.*rdotq, nsrc) + ...
%                       bsxfun(@times, D.*rdotn, qsrc), ...
%                       1);
%   end
