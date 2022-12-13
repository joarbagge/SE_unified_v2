function x = SE2P_ifftnd(x, xl, x0, opt)
%SE2P_ifftnd  3-dimensional inverse discrete Fourier Transform.
%
%   x = SE2P_ifftnd(x, xl, x0, opt)
%      Returns the 3-dimensional inverse DFT of the 3D array x
%      with two periodic directions and one free direction.
%      The input xl contains the "local pad" data and x0
%      contains the zero mode.
%
%   Note: For simplicity, the local pad contains the zero mode as
%   well but it will be updated by the true zero mode data.
%
%   Input parameters:
%       :param x: input vector of size (Mx,My,Mz*sg)
%       :param xl: local pad vector of size (nx,ny,Mz*sl)
%       :param x0: zero mode vector of size (1,1,Mz*s0)
%       :param opt.grid: vector containing the sizes [Mx,My,Mz]
%       :param opt.grid3: size of third dimension of output (i.e. Mz)
%       :param opt.actual_upsampling_zero: upsampling factor for the zero mode
%       :param opt.actual_upsampling_local: upsampling factor on the "local pad"
%       :param opt.actual_upsampling_global: upsampling factor for the whole domain
%       :param opt.local_modes1: list of modes to include in "local pad", x direction
%       :param opt.local_modes2: list of modes to include in "local pad", y direction
%
%   Here, nx=numel(opt.local_modes1), ny=numel(opt.local_modes2),
%   sg=opt.actual_upsampling_global, sl=opt.actual_upsampling_local,
%   s0=opt.actual_upsampling_zero.
%
%   :returns: **x** -- output vector of size (Mx,My,Mz)
%
%   The periodic directions are x and y, while the free direction is z.

if isempty(x), return; end

sg = opt.actual_upsampling_global;
sl = opt.actual_upsampling_local;
s0 = opt.actual_upsampling_zero;

if sg == sl && sg == s0 % all upsampling factors are the same
  x = ifftn(x);
  x = x(1:opt.grid(1), 1:opt.grid(2), 1:opt.grid3); % truncate
else
  % Zero mode
  F0 = ifft(x0, round(opt.grid3*s0));               % 1D ifft in z (zero mode)

  % Local pad
  Fl = ifft(xl, round(opt.grid3*sl), 3);            % 1D ifft in z (local pad)

  % Global domain
  F = ifft(x, round(opt.grid3*sg), 3);              % 1D ifft in z (global)

  % Put together
  L1 = opt.local_modes1;
  L2 = opt.local_modes2;
  Fz = zeros(opt.grid(1), opt.grid(2), opt.grid3);
  % Truncate all data to 1:opt.grid3
  Fz(:,:,:) = F(:,:,1:opt.grid3);
  Fz(L1,L2,:) = Fl(:,:,1:opt.grid3);
  Fz(1,1,:) = F0(1:opt.grid3);

  x = ifft2(Fz, opt.grid(1), opt.grid(2));          % 2D ifft in x and y
end

end
