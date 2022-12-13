function xvec = recenter_points(xvec, box, periodicity)
% xvec = recenter_points(xvec, box, periodicity)
%
% Recenter points in primary box in periodic directions

    if numel(xvec) == 0
        return;
    end

    if nargin < 3
        periodicity = numel(box);
    end

    for d=1:min(periodicity, numel(box))
        xvec(:,d) = xvec(:,d) + box(d);
        xvec(:,d) = mod(xvec(:,d), box(d));
    end

end
