function ker_fft = SE0P_base_precompute_kernel_fft(opt)
%SE0P_base_precompute_kernel_fft  Precompute FFT of 0P (free-space) kernel
%
%   ker_fft = SE0P_base_precompute_kernel_fft(opt)
%
%   Precomputes the kernel on the upsampled grid, and then
%   truncates the result to the padded grid. This will precompute
%   the harmonic or biharmonic kernel, which the full kernel
%   (e.g. stokeslet) can be related to.
%
%   Input parameters:
%       :param opt: option struct, see "help SE0P_base_fourier_space"
%
%   :returns: **ker_fft.kernel** -- precomputed data

opt = SE0P_check_options(opt);

if isfield(opt, 'pre_precision') && ischar(opt.pre_precision) && strcmp(opt.pre_precision, 'quad')
  harmonic = @(G,B,R) SE0P_precompute_harmonic_kernel_fft_quad_mex(G,B,R);
  biharmonic = @(G,B,R) SE0P_precompute_biharmonic_kernel_fft_quad_mex(G,B,R);
elseif isfield(opt, 'pre_precision') && isnumeric(opt.pre_precision) && opt.pre_precision > 0
  prec = ceil(opt.pre_precision);
  harmonic = @(G,B,R) SE0P_precompute_harmonic_kernel_fft_mpfr_mex(G,B,R,prec);
  biharmonic = @(G,B,R) SE0P_precompute_biharmonic_kernel_fft_mpfr_mex(G,B,R,prec);
else
  harmonic = @(G,B,R) SE0P_precompute_harmonic_kernel_fft_mex(G,B,R);
  biharmonic = @(G,B,R) SE0P_precompute_biharmonic_kernel_fft_mex(G,B,R);
end

% Select appropriate MEX file
if strcmp(opt.kernel, 'laplace')
  mexfun = harmonic; kernel_name = 'harmonic';
elseif strcmp(opt.kernel, 'rotlet')
  if ~isfield(opt, 'vico_variant') || opt.vico_variant == 1
    mexfun = harmonic; kernel_name = 'harmonic';
  elseif opt.vico_variant == 2
    mexfun = biharmonic; kernel_name = 'biharmonic';
  end
elseif strcmp(opt.kernel, 'stokeslet')
  if ~isfield(opt, 'vico_variant') || opt.vico_variant == 1
    mexfun = biharmonic; kernel_name = 'biharmonic';
  elseif opt.vico_variant == 3
    mexfun = harmonic; kernel_name = 'harmonic';
  end
elseif strcmp(opt.kernel, 'stresslet')
  mexfun = biharmonic; kernel_name = 'biharmonic';
else
  error('Unsupported kernel: %s', opt.kernel);
end

% Compute modified Green's function FFT on upsampled grid using MEX file
GR = mexfun(opt.upsampled_grid, opt.upsampled_box, opt.greens_truncation_R);

% Strategy is now IFFT -> truncate -> FFT

% Truncate in real space by removing the center part of GR
% (corresponding to the part furthest away from the origin).
% Indexing GR directly like below is slightly faster than doing fftshift
% to center GR first (which would make the indexing a bit simpler).
GR = ifftn(GR);
M_pad = opt.padded_grid;
N = ceil(M_pad/2);
n = floor(M_pad/2);
GR = GR([1:N(1), end-n(1)+1:end], ...
        [1:N(2), end-n(2)+1:end], ...
        [1:N(3), end-n(3)+1:end]);

% Apply Gaussian to smoothen out the function in the transition region,
% but only for biharmonic kernel
if strcmp(kernel_name, 'biharmonic') && ~isempty(opt.transition_level) ...
    && isfinite(opt.transition_level)
  if opt.upsampling_convolution > 0
    convolution_grid = 2*ceil(opt.upsampling_convolution*opt.extended_grid/2); % round up
  else
    convolution_grid = 2*ceil((opt.padded_grid + opt.upsampling_convolution)/2); % round up
  end
  if opt.transition_level == -1
    myfun = @(n) bump(n);
  else
    myfun = @(n) gaussian(opt.transition_level, n);
  end
  % The convolution grid can never be larger than the padded grid
  convolution_grid = min(convolution_grid, opt.padded_grid);
  M_conv = convolution_grid;
  NC = ceil(M_conv/2);
  nC = floor(M_conv/2);
  npoints = size(GR,1) - (nC(1)+NC(1));
  G = myfun(npoints);
  GR([NC(1)+1:end-nC(1)],:,:) = GR([NC(1)+1:end-nC(1)],:,:) .* ...
        repmat(G(:), [1,size(GR,2),size(GR,3)]);
  npoints = size(GR,2) - (nC(2)+NC(2));
  G = myfun(npoints);
  GR(:,[NC(2)+1:end-nC(2)],:) = GR(:,[NC(2)+1:end-nC(2)],:) .* ...
        repmat(G(:)', [size(GR,1),1,size(GR,3)]);
  npoints = size(GR,3) - (nC(3)+NC(3));
  G = myfun(npoints);
  GR(:,:,[NC(3)+1:end-nC(3)]) = GR(:,:,[NC(3)+1:end-nC(3)]) .* ...
        repmat(reshape(G, [1,1,numel(G)]), [size(GR,1),size(GR,2),1]);
end

% Go to Fourier space
GR = real(fftn(GR)); % result should be real (so save some memory)
ker_fft.kernel = GR;

end

function y = gaussian(trunc_level, npoints)
  x = 0:npoints-1;
  a2 = -4*log(trunc_level)/(npoints-1)^2;
  y = exp(-a2*x.^2) + exp(-a2*(x-npoints+1).^2);
end

function y = bump(npoints)
  sig = (npoints-1)/2;
  x = 0:npoints-1;
  A = exp(1);
  y = zeros(size(x));
  y(1:npoints/2) = A*exp(-1./(1 - (x(1:npoints/2)/sig).^2));
  y(npoints/2+1:end) = A*exp(-1./(1 - ((npoints-1-x(npoints/2+1:end))/sig).^2));
end
