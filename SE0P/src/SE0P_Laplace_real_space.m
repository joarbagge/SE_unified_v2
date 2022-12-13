function out = SE0P_Laplace_real_space(x, f, opt)
%SE0P_Laplace_real_space  Compute real-space part of the
%0-periodic (free-space) Ewald sum for the electrostatic potential.
%
%   out = SE0P_Laplace_real_space(x, f, opt)
%
%   Computes the electrostatic potential from N point charges.
%
%   Input parameters:
%       :param x: source locations (N×3)
%       :param f: source charges (N×1)
%       :param opt: option struct, see below
%
%   Valid options for the opt struct are as follows.
%
%   Options controlling the output:
%       :param opt.return_pot: return potential (default: true)
%       :param opt.return_potgrad: return potential gradient
%           (default: false)
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
%           points selected by opt.eval_idx (empty if opt.return_pot
%           is false)
%       :returns: **out.potgrad** (array) -- potential gradient at
%           source points selected by opt.eval_idx (empty if
%           opt.return_potgrad is false)
%       :returns: **out.extpot** (array) -- potential evaluated at
%           external (non-source) points selected by opt.eval_ext_x
%           (empty if opt.return_pot is false)
%       :returns: **out.extpotgrad** (array) -- potential gradient at
%           external points selected by opt.eval_ext_x (empty if
%           opt.return_potgrad is false)
%       :returns: **out.walltime** -- struct containing timings
%       :returns: **out.opt** -- final option struct, modified by
%           the program; actual parameters used by the algorithm

assert(isfield(opt, 'box'), 'Periodic cell size "box" must be given in opt struct');
assert(isfield(opt, 'xi'), 'Ewald decomposition parameter "xi" must be given in opt struct');
assert(isfield(opt, 'rc'), 'Cutoff radius "rc" must be given in opt struct');

if ~isfield(opt,'return_pot'), opt.return_pot = true; end
if ~isfield(opt,'return_potgrad'), opt.return_potgrad = false; end
if ~isfield(opt,'eval_idx'), opt.eval_idx = 'default'; end
if ~isfield(opt,'eval_ext_x'), opt.eval_ext_x = []; end

% Prepare output struct
out.pot = [];
out.potgrad = [];
out.extpot = [];
out.extpotgrad = [];
out.walltime = struct('pot',0, 'potgrad',0, 'extpot',0, 'extpotgrad',0);
out.opt = opt;

% Evaluate
% FIXME: Would be better if the MEX file evaluated only in
% opt.eval_idx instead of having to slice afterwards.
%   Equivalent MATLAB code found at the bottom (in comment)
if opt.return_pot && (strcmp(opt.eval_idx, 'default') || numel(opt.eval_idx) > 0)
  ts = tic;
  out.pot = SE0P_Laplace_real_cut_cell_mex(x, f, opt.rc, opt.xi, opt.box);
  if ~strcmp(opt.eval_idx, 'default')
    out.pot = out.pot(opt.eval_idx,:);
  end
  out.walltime.pot = toc(ts);
end
if opt.return_potgrad && (strcmp(opt.eval_idx, 'default') || numel(opt.eval_idx) > 0)
  ts = tic;
  out.potgrad = SE0P_Laplace_real_cut_cell_grad_mex(x, f, opt.rc, opt.xi, opt.box);
  if ~strcmp(opt.eval_idx, 'default')
    out.potgrad = out.potgrad(opt.eval_idx,:);
  end
  out.walltime.potgrad = toc(ts);
end
% FIXME: Would be better to have fast (cell-based) code for
% external evaluation, using direct summation for now.
if opt.return_pot && numel(opt.eval_ext_x) > 0
  ts = tic;
  opt.layers = 1;
  out.extpot = SE0P_Laplace_direct_real_cut_ext_mex(x, f, opt);
  out.walltime.extpot = toc(ts);
end
if opt.return_potgrad && numel(opt.eval_ext_x) > 0
  ts = tic;
  opt.layers = 1;
  out.extpotgrad = SE0P_Laplace_direct_real_cut_ext_grad_mex(x, f, opt);
  out.walltime.extpotgrad = toc(ts);
end

out.walltime.total = sum(struct2array(out.walltime));

end

% Equivalent MATLAB code for potential (not external):
%   N = size(x, 1);
%   u = zeros(N, 1);
%   [idx, d] = rangesearch(x, x, opt.rc);
%   for target=opt.eval_idx
%     nnb_idx = idx{target}(2:end);
%     nnb_r = d{target}(2:end);
%     u(target) = sum(f(nnb_idx)' .* erfc(opt.xi*nnb_r)./nnb_r);
%   end
%
% Equivalent MATLAB code for gradient (not external):
%   N = size(x, 1);
%   du = zeros(N, 3);
%   [idx, d] = rangesearch(x, x, opt.rc);
%   c = 2/sqrt(pi)*opt.xi;
%   for target=opt.eval_idx
%     nnb_idx = idx{target}(2:end);
%     nnb_r = d{target}(2:end);
%     r = bsxfun(@minus, x(target,:), x(idx{target}(2:end),:));
%     u1 = c*exp(-(opt.xi^2*nnb_r.^2));
%     u2 = erfc(opt.xi*nnb_r) ./ nnb_r;
%     u = f(nnb_idx)' .* (u1+u2) ./ nnb_r.^2;
%     du(target,:) = -sum(bsxfun(@times, u', r), 1);
%   end
