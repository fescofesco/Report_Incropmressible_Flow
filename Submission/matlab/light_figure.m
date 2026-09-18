function f = light_figure()
% LIGHT_FIGURE  Create an invisible figure forced to the light (printable) theme.
%
%   f = light_figure()
%
%   MATLAB R2025a introduced figure themes, and when MATLAB runs headless
%   (-batch, as used to regenerate the report figures) the exported PNG comes
%   out with a BLACK background and near-white axis text. Embedded in the
%   report that is both unreadable in print and a waste of toner.
%
%   Setting 'Color','w' alone is not enough: it whitens the background but
%   leaves the dark theme's light-grey text, so the title and tick labels
%   become almost invisible. theme(f,'light') sets background *and*
%   foreground colours together, which is what we want.
%
%   The theme() call is guarded because releases before R2025a have no such
%   function -- there the light appearance is already the default.

    f = figure('Visible', 'off', 'Color', 'w');

    % The figures are exported large and then scaled down to ~0.8\textwidth
    % in the report, so MATLAB's default 10 pt axis font ends up barely
    % legible in print. Bump it for every axes created in this figure.
    set(f, 'DefaultAxesFontSize', 13);
    try
        theme(f, 'light');
    catch
        % Pre-R2025a: no themes, figures are light by default. Nothing to do.
    end
end
