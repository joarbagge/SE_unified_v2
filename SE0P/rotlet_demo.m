% Spectral Ewald 0P (free space), Rotlet demo with accuracy check

clear
if isempty(getenv('SE_reset_ROOT')), run('../init.m'); end;
rng(1);

N = 100; % number of source particles
box = [1.0 1.1 1.2]; % side lengths of periodic box
[x, t] = SE_random_system(N, box, 3); % generate random sources

% Ewald parameters
tol = 1e-14; % absolute RMS error tolerance
rc = 0.7; % Set the cutoff radius rc, the rest is automatic
opt = SE_Rotlet_params(box, t, 'AbsTol', tol, 'rc', rc);
assert(opt.rc <= min(box), 'rc (%g) cannot be larger than min(box) (%g)', opt.rc, min(box));

fprintf('Computing direct summation for reference...');
u_ref = compute_reference_solution(x, t);
fprintf(' done.\n');

fprintf('Computing Spectral Ewald solution...');
ts = tic();
outr = SE0P_Rotlet_real_space(x, t, opt);
pre_kernel = SE0P_Rotlet_precompute_kernel_fft(opt);
outf = SE0P_Rotlet_fourier_space(x, t, opt, pre_kernel);
u = outr.pot + outf.pot;
time = toc(ts);
fprintf(' done in %.6g sec.\n', time);

abs_rms_error = sqrt(3)*rms(u(:) - u_ref(:));
rel_rms_error = abs_rms_error / (sqrt(3)*rms(u_ref(:)));
fprintf('  RMS error: %.6g (abs), %.6g (rel)\n', abs_rms_error, rel_rms_error);

% ------------------------------------------------------------------------------
function u = compute_reference_solution(x, t)
% Compute reference solution using direct summation
  u = SE0P_Rotlet_direct_full_mex(x, t, struct());
end
