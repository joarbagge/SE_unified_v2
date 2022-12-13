% Spectral Ewald 1P, Stokeslet demo with accuracy check

clear
if isempty(getenv('SE_reset_ROOT')), run('../init.m'); end;
rng(1);

N = 50; % number of source particles
box = [1.0 1.1 1.2]; % side lengths of periodic box
[x, f] = SE_random_system(N, box, 3); % generate random sources

% Ewald parameters
tol = 1e-14; % absolute RMS error tolerance
rc = 0.7; % Set the cutoff radius rc, the rest is automatic
opt = SE_Stokeslet_params(box, f, 'periodicity', 1, 'AbsTol', tol, 'rc', rc);
assert(opt.rc <= min(box), 'rc (%g) cannot be larger than min(box) (%g)', opt.rc, min(box));
opt.stokeslet_k0_constant = 1; % arbitrary choice of constant

fprintf('Computing direct summation for reference...');
u_ref = compute_reference_solution(x, f, box, opt.xi, 1e-16, opt.stokeslet_k0_constant);
fprintf(' done.\n');

fprintf('Computing Spectral Ewald solution...');
ts = tic();
outr = SE1P_Stokeslet_real_space(x, f, opt);
outf = SE1P_Stokeslet_fourier_space(x, f, opt);
us = SE_Stokeslet_self_term(x, f, opt);
u = outr.pot + outf.pot + us;
time = toc(ts);
fprintf(' done in %.6g sec.\n', time);

abs_rms_error = sqrt(3)*rms(u(:) - u_ref(:));
rel_rms_error = abs_rms_error / (sqrt(3)*rms(u_ref(:)));
fprintf('  RMS error: %.6g (abs), %.6g (rel)\n', abs_rms_error, rel_rms_error);

% ------------------------------------------------------------------------------
function u = compute_reference_solution(x, f, box, xi, tol, c)
% Compute reference solution using direct summation
  opt = SE_Stokeslet_params(box, f, 'periodicity', 1, 'AbsTol', tol, 'xi', xi);
  assert(opt.rc <= min(box), 'rc (%g) cannot be larger than min(box) (%g)', opt.rc, min(box));
  % Real space
  opt.layers = 1; % assuming rc <= min(box), this is enough for real-space part
  ur = SE1P_Stokeslet_direct_real_cut_mex(x, f, opt);
  % Fourier space
  opt.layers = ceil(max(opt.grid_res*box)/2); % M/2 layers
  opt.stokeslet_k0_constant = c;
  uf = SE1P_Stokeslet_direct_fourier_mex(x, f, opt);
  uf_k0 = SE1P_Stokeslet_direct_fourier_k0_mex(x, f, opt);
  us = SE1P_Stokeslet_direct_self_mex(x, f, opt);
  % Put together
  u = ur + uf + uf_k0 + us;
end
