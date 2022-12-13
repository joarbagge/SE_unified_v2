function [x,f] = SE_random_system(N, box, varargin)
%SE_random_system  Generate particle systems of different dimensions.
%
%   [x,f] = SE_random_system(N, [L1 L2 L3]) or
%   [x,f] = SE_random_system(N, [L1 L2 L3], 1) generates a
%   one-dimensional random system on a box of size [L1 L2 L3].
%
%   [x,f] = SE_random_system(N, [L1 L2 L3], DIM) generates a
%   DIM-dimensional random system.
%
%   [x,f] = SE_random_system(..., 'neutral') generates a charge
%   neutral random system, i.e. a system such that sum(f) is zero.

% Default values
dim = 1;
charge_neutral = false;

for k=1:numel(varargin)
  if isa(varargin{k}, 'double')
    dim = varargin{k};
  elseif strcmp(varargin{k}, 'neutral')
    charge_neutral = true;
  end
end

x = bsxfun(@times, rand(N, 3), box);
f = 1 - 2*rand(N, dim);

if charge_neutral
  f = f - repmat(mean(f,1), N, 1); % neutrality
end

end
