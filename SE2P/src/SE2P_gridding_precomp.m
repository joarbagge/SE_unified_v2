function pre_grid = SE2P_gridding_precomp(x, opt)
%SE2P_gridding_precomp  Precomputation of fast gridding data
%
%   pre_grid = SE2P_gridding_precomp(x, opt)
%
%   Input parameters:
%       :param x: source locations (N×3)
%       :param opt: option struct, see "help SE2P_base_fourier_space"
%
%   :returns: **pre_grid** -- structure containing precomputed data

x = recenter_points(x, opt.box, 2);
analytical_diff = strcmp(opt.kernel, 'laplace') ...
                  && opt.return_potgrad ...
                  && ~opt.use_ik_diff;

if strcmp(opt.window, 'gaussian')
  if analytical_diff
    [zx, zy, zz, zfx, zfy, zfz, idx] = SE_fgg_expand_all_force_mex_2p(x, opt);
  else
    [zx, zy, zz, idx] = SE_fgg_expand_all_mex_2p(x, opt);
  end
  pre_grid.zs = SE_fgg_base_gaussian_mex_2p(opt);
else
  if analytical_diff
    [zx, zy, zz, zfx, zfy, zfz, idx] = SE_fkg_expand_all_force_mex_2p(x, opt);
  else
    [zx, zy, zz, idx] = SE_fkg_expand_all_mex_2p(x, opt);
  end
end

[idx, s] = sort(idx);
zx = zx(1:opt.window_P, s);
zy = zy(1:opt.window_P, s);
zz = zz(1:opt.window_P, s);

pre_grid.analytical_diff = analytical_diff;
pre_grid.zx = zx;
pre_grid.zy = zy;
pre_grid.zz = zz;
pre_grid.idx = idx;
pre_grid.perm = s';
pre_grid.iperm(s) = 1:length(s);

if analytical_diff
  zfx = zfx(1:opt.window_P, s);
  zfy = zfy(1:opt.window_P, s);
  zfz = zfz(1:opt.window_P, s);
  pre_grid.zfx = zfx;
  pre_grid.zfy = zfy;
  pre_grid.zfz = zfz;
end

end
