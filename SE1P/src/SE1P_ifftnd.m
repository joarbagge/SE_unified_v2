function x = SE1P_ifftnd(x, xl, x0, opt)
%SE1P_ifftnd  3-dimensional inverse discrete Fourier Transform.
%
%   x = SE1P_ifftnd(x, xl, x0, opt)
%      Returns the 3-dimensional inverse DFT of the 3D array x
%      with one periodic direction and two free directions.
%      The input xl contains the "local pad" data and x0
%      contains the zero mode.
%
%   Input parameters:
%       :param x: input vector of size (Mx,My*sg(1),Mz*sg(2))
%       :param xl: local pad vector of size (nx,My*sl(1),Mz*sl(2))
%       :param x0: zero mode vector of size (1,My*s0(1),Mz*s0(2))
%       :param opt.grid: vector containing the sizes [Mx,My,Mz]
%       :param opt.grid2: size of second dimension of output (i.e. My)
%       :param opt.grid3: size of third dimension of output (i.e. Mz)
%       :param opt.actual_upsampling_zero: upsampling factor for the zero mode
%       :param opt.actual_upsampling_local: upsampling factor on the "local pad"
%       :param opt.actual_upsampling_global: upsampling factor for the whole domain
%       :param opt.local_modes1: list of modes to include in "local pad"
%
%   Here, nx=numel(opt.local_modes1), sg=opt.actual_upsampling_global,
%   sl=opt.actual_upsampling_local, s0=opt.actual_upsampling_zero.
%
%   :returns: **x** -- output vector of size (Mx,My,Mz)
%
%   The periodic direction is x, while the free directions are y and z.

if isempty(x), return; end

sg = opt.actual_upsampling_global;
sl = opt.actual_upsampling_local;
s0 = opt.actual_upsampling_zero;

if all(sg == sl) && all(sg == s0) % all upsampling factors are the same
  x = ifftn(x);
  x = x(1:opt.grid(1), 1:opt.grid2, 1:opt.grid3);           % truncate
else
  % Zero mode
  F0 = ifft2(x0, round(opt.grid2*s0(1)), round(opt.grid3*s0(2))); % 2D ifft in y and z (zero mode)

  % Local pad
  Fl = ifft(ifft(xl, round(opt.grid2*sl(1)), 2), ...
            round(opt.grid3*sl(2)), 3);                     % 2D ifft in y and z (local pad)

  % Global domain
  F = ifft(ifft(x, round(opt.grid2*sg(1)), 2), ...
           round(opt.grid3*sg(2)), 3);                      % 2D ifft in y and z (global)

  % Put together
  L1 = opt.local_modes1;
  Glob = setdiff(1:opt.grid(1), [1 L1]);                    % "global pad"
  Fyz = zeros(opt.grid(1), opt.grid2, opt.grid3);
  Fyz(L1,:,:) = Fl(:, 1:opt.grid2, 1:opt.grid3);            % truncate to (nx,My,Mz)
  Fyz(1,:,:) = F0(1:opt.grid2, 1:opt.grid3);                % truncate to (1,My,Mz)
  Fyz(Glob,:,:) = F(Glob, 1:opt.grid2, 1:opt.grid3);        % truncate to (*,My,Mz)

  x = ifft(Fyz, opt.grid(1));                               % 1D ifft in x
end

end
