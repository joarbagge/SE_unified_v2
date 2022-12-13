function ker_fft = SE0P_Stokeslet_precompute_kernel_fft(opt)
%SE0P_Stokeslet_precompute_kernel_fft  Precompute FFT of 0P (free-space) kernel
%
%   ker_fft = SE0P_Stokeslet_precompute_kernel_fft(opt)
%
%   Precomputes the kernel on the upsampled grid, and then
%   truncates the result to the padded grid. Stokeslet version.
%
%   Input parameters:
%       :param opt: option struct, see "help SE0P_base_fourier_space"
%
%   :returns: **ker_fft.kernel** -- precomputed data

opt.kernel = 'stokeslet';
ker_fft = SE0P_base_precompute_kernel_fft(opt);

end
