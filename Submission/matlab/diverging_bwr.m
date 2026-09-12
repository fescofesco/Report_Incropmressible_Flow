function cmap = diverging_bwr(n)
% DIVERGING_BWR  Blue-white-red diverging colormap (no toolbox dependency).
%
%   cmap = diverging_bwr(n) returns an n-by-3 colormap running from blue
%   (low) through white (mid) to red (high), suitable for signed fields
%   (e.g. radial velocity) plotted with a symmetric color axis centred on
%   zero via caxis([-lim, lim]).

    if nargin < 1
        n = 256;
    end

    blue  = [0.13, 0.40, 0.85];
    white = [1.00, 1.00, 1.00];
    red   = [0.80, 0.10, 0.10];

    half = floor(n / 2);
    lower = [linspace(blue(1), white(1), half)', ...
              linspace(blue(2), white(2), half)', ...
              linspace(blue(3), white(3), half)'];
    upper = [linspace(white(1), red(1), n - half)', ...
              linspace(white(2), red(2), n - half)', ...
              linspace(white(3), red(3), n - half)'];
    cmap = [lower; upper];
end
