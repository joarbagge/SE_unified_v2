% Spectral Ewald 0P (free space), Laplace potential demo with accuracy check

clear
if isempty(getenv('SE_reset_ROOT')), run('../init.m'); end;
rng(1);

N = 100; % number of source particles
box = [1.0 1.1 1.2]; % side lengths of periodic box
[x, f] = SE_random_system(N, box, 'neutral'); % generate random sources

% Ewald parameters
tol = 1e-14; % absolute RMS error tolerance
rc = 0.7; % Set the cutoff radius rc, the rest is automatic
opt = SE_Laplace_params(box, f, 'AbsTol', tol, 'rc', rc);
assert(opt.rc <= min(box), 'rc (%g) cannot be larger than min(box) (%g)', opt.rc, min(box));

fprintf('Computing direct summation for reference...');
u_ref = compute_reference_solution(x, f);
fprintf(' done.\n');

fprintf('Computing Spectral Ewald solution...');
ts = tic();
outr = SE0P_Laplace_real_space(x, f, opt);
pre_kernel = SE0P_Laplace_precompute_kernel_fft(opt);
outf = SE0P_Laplace_fourier_space(x, f, opt, pre_kernel);
us = SE_Laplace_self_term(x, f, opt);
u = outr.pot + outf.pot + us;
time = toc(ts);
fprintf(' done in %.6g sec.\n', time);

abs_rms_error = rms(u - u_ref);
rel_rms_error = abs_rms_error / rms(u_ref);
fprintf('  RMS error: %.6g (abs), %.6g (rel)\n', abs_rms_error, rel_rms_error);

% ------------------------------------------------------------------------------
function u = compute_reference_solution(x, f)
% Compute reference solution using direct summation
  u = SE0P_Laplace_direct_full_mex(x, f, struct());
end
