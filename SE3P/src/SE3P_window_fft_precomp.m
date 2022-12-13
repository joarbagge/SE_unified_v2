function pre_fft = SE3P_window_fft_precomp(opt)

pw = opt.window_scaling_power; % should be 2 or 1 typically
sh = opt.window_shape;
w = opt.window_halfwidth;
h = opt.h;

if strcmp(opt.window, 'gaussian')
  error('FFT precomputing not supported for Gaussian window');
elseif strcmp(opt.window, 'expsemicirc')
  window = @(x) expsemicirc(x, sh, w);
elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
  window = @(x) kaiser_exact(x, sh, w, opt.kaiser_scaling);
else
  error('Unsupported window function');
end

F = cell(1,3);
for d=1:3
  L = opt.box(d);
  x = 0:h:(L-h);
  f = h*window(fftshift(x-L/2));
  F{d} = fft(f).^pw;
end

[F1, F2, F3] = ndgrid(F{1}, F{2}, F{3});
F = F1.*F2.*F3; % tensor product of spatial directions
pre_fft.window = 1./F;

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
