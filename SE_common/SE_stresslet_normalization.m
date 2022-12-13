function [q, n] = SE_stresslet_normalization(q, n)
%SE_stresslet_normalization  Normalize stresslet sources.
%
%   [q, n] = SE_stresslet_normalization(q, n)
%   modifies q such that sum(q(:).*n(:)) is zero.
%
%   This guarantees that the resulting flow field is divergence-free.

% Strategy: modify only q(1,k), where k maximizes abs(n(1,k)).

N = size(q,1);
if N == 0, return; end

[~,k] = max(abs(n(1,:)));

% Compute sum excluding component (1,k)
qn = q .* n;
qn(1,k) = 0;
S = sum(qn(:));

% Modify q(1,k)
q(1,k) = -S/n(1,k);

end
