function G = SE0P_Stresslet_scaling(H, opt, pre_kernel, pre_window)
%SE0P_Stresslet_scaling  Fourier-space scaling for 0-periodic (free-space) stresslet
%
%   G = SE0P_Stresslet_scaling(H, opt, pre_kernel, pre_window)
%
%   Input parameters:
%       :param H: input data in Fourier space (cell array)
%       :param opt: structure with Ewald options
%       :param pre_kernel: mandatory precomputation structure:
%       :param pre_kernel.kernel: Fourier transform of the biharmonic Green's function
%       :param pre_window: optional precomputation structure:
%       :param pre_window.window: W^(-pw) where W is the Fourier transform of the window function
%
%   :returns: **G** -- output data in Fourier space (cell array)

if nargin < 4
  pre_window = [];
end
pw = opt.window_scaling_power; % should be 2 or 1 typically
w = opt.window_halfwidth;

% Compute k-vectors
K = cell(3,1);
[k1, k2, k3] = SE_k_vectors(opt.padded_grid, opt.padded_box, 'shifted');
[K{1}, K{2}, K{3}] = ndgrid(k1, k2, k3);
KK = K{1}.^2 + K{2}.^2 + K{3}.^2;

% The scaling is given by: Green × Screening × Window^(-2),
% each factor representing the Fourier transform of a function.
% For the Gaussian window, we combine Screening×Window^(-2) into
% a single factor (Combo).

% For the stresslet, we first use the biharmonic Green's function,
% and then use a relation between it and the stresslet.

Biharmonic = pre_kernel.kernel;
TMP = KK/(4*opt.xi^2);
if strcmp(opt.window, 'gaussian') % special treatment for Gaussian window
  eta = 2*(w*opt.xi)^2 / opt.window_shape;
  % Screening and Window factors combined
  TMP = (1+TMP) .* exp(-TMP*(1-eta*pw/2)); % Combo
  TMP = Biharmonic .* TMP; % multiply Combo by Green to get Scaling
  % MEX file available for Gaussian: SE0P_stresslet_fast_fs_k_scaling
else
  TMP = (1+TMP) .* exp(-TMP); % Screening
  if ~isempty(pre_window)
    Window_m2 = pre_window.window;
  elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
    % Compute Window_m2 on the fly
    b2 = opt.window_shape^2;
    f1 = kaiser_exact_ft(k1.^2, b2, w, opt.kaiser_scaling);
    f2 = kaiser_exact_ft(k2.^2, b2, w, opt.kaiser_scaling);
    f3 = kaiser_exact_ft(k3.^2, b2, w, opt.kaiser_scaling);
    [F1, F2, F3] = ndgrid(f1.^pw, f2.^pw, f3.^pw);
    Window_m2 = 1./(F1.*F2.*F3); % tensor product of spatial directions
    clear F1 F2 F3;
  else
    error('Precomputed Fourier transform necessary for the %s window', opt.window);
  end
  TMP = Biharmonic .* TMP .* Window_m2; % Scaling
  clear Window_m2;
end
clear Biharmonic;
% Use relation between stresslet and biharmonic Green's function
G = cell(3,1);
HdotK = cell(3,1);
for j=1:3
  HdotK{j} = H{j,1}.*K{1} + H{j,2}.*K{2} + H{j,3}.*K{3};
end
KdotHdotK = K{1}.*HdotK{1} + K{2}.*HdotK{2} + K{3}.*HdotK{3};
for j=1:3
  KdotH = K{1}.*H{1,j} + K{2}.*H{2,j} + K{3}.*H{3,j};
  traceH = H{1,1} + H{2,2} + H{3,3};
  vector_part = -1i*(HdotK{j} + KdotH + K{j}.*traceH).*KK + 2i*K{j}.*KdotHdotK;
  G{j} = TMP .* vector_part;
end

% ------------------------------------------------------------------------------
function F = kaiser_exact_ft(k2, b2, w, scaling)

t = sqrt(b2 - k2*w^2);
F = 2*w*sinh(t)./t * scaling;
