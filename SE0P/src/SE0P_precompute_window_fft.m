function pre_fft = SE0P_precompute_window_fft(opt)

pw = opt.window_scaling_power; % should be 2 or 1 typically
sh = opt.window_shape;
L = opt.upsampled_box(1);
M = opt.upsampled_grid(1);
h = L/M;
w = opt.window_halfwidth * h/opt.h; % for upsampled grid

if strcmp(opt.window, 'gaussian')
  error('FFT precomputing not supported for Gaussian window');
elseif strcmp(opt.window, 'expsemicirc')
  window = @(x) expsemicirc(x, sh, w);
elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
  window = @(x) kaiser_exact(x, sh, w, opt.kaiser_scaling);
else
  error('Unsupported window function');
end

x = 0:h:(L-h);
f = h*window(x-L/2);
n = opt.padded_grid(1);
N = opt.upsampled_grid(1);
start = round((N-n)/2);
f = f(start + (1:n));
F = 1./fft(f).^pw;
[F1, F2, F3] = ndgrid(F, F, F);
Fprod = F1.*F2.*F3; % tensor product of spatial directions
pre_fft.window = fftshift(Fprod); % TODO: why fftshift?

end

% ------------------------------------------------------------------------------
function z = expsemicirc(x, beta, w)

t = sqrt(1-(x/w).^2);
z = exp(beta*(t-1)) .* (abs(x) <= w);

end

% ------------------------------------------------------------------------------
function z = kaiser_exact(x, beta, w, scaling)

t = sqrt(1-(x/w).^2);
I = besseli(0,beta*t) * scaling;
z = I .* (abs(x) <= w);

end
