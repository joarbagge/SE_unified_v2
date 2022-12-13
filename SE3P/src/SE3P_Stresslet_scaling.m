function G = SE3P_Stresslet_scaling(H, opt, pre_fft)
%SE3P_Stresslet_scaling  Fourier-space scaling for 3-periodic stresslet
%
%   G = SE3P_Stresslet_scaling(H, opt, pre_fft)
%
%   Input parameters:
%       :param H: input data in Fourier space (cell array)
%       :param opt: structure with Ewald options
%       :param pre_fft: optional precomputation structure:
%       :param pre_fft.window: W^(-pw) where W is the Fourier transform of the window function
%
%   :returns: **G** -- output data in Fourier space (cell array)

if nargin < 3
  pre_fft = [];
end
pw = opt.window_scaling_power; % should be 2 or 1 typically
w = opt.window_halfwidth;

% Compute k-vectors
K = cell(3,1);
[k1, k2, k3] = SE_k_vectors(opt.grid, opt.box);
[K{1}, K{2}, K{3}] = ndgrid(k1, k2, k3);
KK = K{1}.^2 + K{2}.^2 + K{3}.^2;

% The scaling is given by: Green × Screening × Window^(-2),
% each factor representing the Fourier transform of a function.
% For the Gaussian window, we combine Screening×Window^(-2) into
% a single factor (Combo).

% For the stresslet, we first use the biharmonic Green's function,
% and then use a relation between it and the stresslet.

% FIXME: Might want to MEX the scaling (then force H complex)

Biharmonic = -8*pi./(KK.*KK);
C = KK/(4*opt.xi^2);
if strcmp(opt.window, 'gaussian') % special treatment for Gaussian window
  eta = 2*(w*opt.xi)^2 / opt.window_shape;
  % Screening and Window factors combined
  Combo = (1+C) .* exp(-C*(1-eta*pw/2));
  Scaling = Biharmonic .* Combo;
else
  Screening = (1+C) .* exp(-C);
  if ~isempty(pre_fft)
    Window_m2 = pre_fft.window;
  elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
    % Compute Window_m2 on the fly
    b2 = opt.window_shape^2;
    f1 = kaiser_exact_ft(k1.^2, b2, w, opt.kaiser_scaling);
    f2 = kaiser_exact_ft(k2.^2, b2, w, opt.kaiser_scaling);
    f3 = kaiser_exact_ft(k3.^2, b2, w, opt.kaiser_scaling);
    [F1, F2, F3] = ndgrid(f1.^pw, f2.^pw, f3.^pw);
    F = F1.*F2.*F3; % tensor product of spatial directions
    Window_m2 = 1./F;
  else
    error('Precomputed Fourier transform necessary for the %s window', opt.window);
  end
  Scaling = Biharmonic .* Screening .* Window_m2;
end
Scaling(1,1,1) = 0; % zero mode skipped here (treated separately)

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
  G{j} = Scaling .* vector_part;
end

% ------------------------------------------------------------------------------
function F = kaiser_exact_ft(k2, b2, w, scaling)

t = sqrt(b2 - k2*w^2);
F = 2*w*sinh(t)./t * scaling;
