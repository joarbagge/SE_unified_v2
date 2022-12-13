function H = SE3P_Laplace_scaling(H, opt, pre_fft)
%SE3P_Laplace_scaling  Fourier-space scaling for 3-periodic electrostatics
%
%   H = SE3P_Laplace_scaling(H, opt, pre_fft)
%
%   Input parameters:
%       :param H: input data in Fourier space (cell array)
%       :param opt: structure with Ewald options
%       :param pre_fft: optional precomputation structure:
%       :param pre_fft.window: W^(-pw) where W is the Fourier transform of the window function
%
%   :returns: **H** -- output data in Fourier space (cell array)

if nargin < 3
  pre_fft = [];
end
pw = opt.window_scaling_power; % should be 2 or 1 typically
w = opt.window_halfwidth;

% Compute k-vectors
[k1, k2, k3] = SE_k_vectors(opt.grid, opt.box);
[K1, K2, K3] = ndgrid(k1, k2, k3);
KK = K1.^2 + K2.^2 + K3.^2;

% The scaling is given by: Green × Screening × Window^(-2),
% each factor representing the Fourier transform of a function.
% For the Gaussian window, we combine Screening×Window^(-2) into
% a single factor (Combo).

% For Laplace, Green = 4*pi/KK.

C = KK/(4*opt.xi^2);
if strcmp(opt.window, 'gaussian') % special treatment for Gaussian window
  eta = 2*(w*opt.xi)^2 / opt.window_shape;
  % Screening and Window factors combined
  % 1. Using Fourier transform of untruncated Gaussian
  Combo = exp(-C*(1-eta*pw/2));
  % 2. Alternative code using Fourier transform of a truncated Gaussian
  %alf = opt.window_shape;
  %f1 = gauss_exact_ft(k1, alf, w);
  %f2 = gauss_exact_ft(k2, alf, w);
  %f3 = gauss_exact_ft(k3, alf, w);
  %[F1, F2, F3] = ndgrid(f1.^pw, f2.^pw, f3.^pw);
  %F = F1.*F2.*F3; % tensor product of spatial directions
  %Window_m2 = 1./F;
  %Combo = Window_m2 .* exp(-C);
  Scaling = 4*pi * Combo ./ KK; % multiply by Green
else
  Screening = exp(-C);
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
  Scaling = 4*pi * Window_m2 .* Screening ./ KK;
end
Scaling(1,1,1) = 0; % zero mode should be zero
H{1} = Scaling .* H{1};

% ik differentiation to compute gradient
if (opt.return_fouriervalgrad || opt.return_gridvalgrad || opt.return_potgrad) ...
    && opt.use_ik_diff
  Hi = 1i * H{1};
  H{2} = Hi .* K1;
  H{3} = Hi .* K2;
  H{4} = Hi .* K3;
end
if ~(opt.return_fourierval || opt.return_gridval || opt.return_pot ...
    || (opt.return_potgrad && ~opt.use_ik_diff))
  H{1} = [];
end

% ------------------------------------------------------------------------------
function F = kaiser_exact_ft(k2, b2, w, scaling)

t = sqrt(b2 - k2*w^2);
F = 2*w*sinh(t)./t * scaling;

% ------------------------------------------------------------------------------
function F = gauss_exact_ft(k, alf, w)

F = exp(-k.^2 * w^2 / (4*alf)); % untruncated Gaussian
F = F .* (0.5*(complex_erf(sqrt(alf)+1i*k*w/(2*sqrt(alf))) + complex_erf(sqrt(alf)-1i*k*w/(2*sqrt(alf)))));
