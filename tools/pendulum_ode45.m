% =========================================================================
%  pendulum_ode45.m  --  Figure 1 (top left) of the paper
%
%  Undamped pendulum      x' = y,   y' = -sin(x)
%
%  Four initial values (x,y) = (0, y0) with y0 in {-2.5, 1.5, 2.0, 2.5}
%  are integrated by MATLAB's ode45 with its default tolerances.
%
%  Three of the four orbits are qualitatively correct.  The one started at
%  y0 = 2.0 (red) makes one loop and then drifts onto a different orbit:
%  the conserved energy  E = 1 - cos(x) + y^2/2  equals 2 for that orbit,
%  which is below the separatrix value 4, so the exact solution is periodic
%  and x must stay bounded.  ode45 loses that.
%
%  Run tools/../repro/figure1.sh to produce the validated counterpart, then
%  tools/plot_pendulum.py (or load the traj_y*.txt files here).
% =========================================================================

f  = @(t, u) [u(2); -sin(u(1))];
y0 = [-2.5, 1.5, 2.0, 2.5];
T  = 20;
col = {'b', 'g', 'r', 'm'};

figure('Color', 'w'); hold on; grid on
for k = 1:numel(y0)
    [~, U] = ode45(f, [0 T], [0; y0(k)]);
    plot(U(:,1), U(:,2), col{k}, 'LineWidth', 1.0, ...
         'DisplayName', sprintf('y(0) = %g', y0(k)));
end
xlabel('x'); ylabel('y');
title('ode45, default tolerances');
legend('Location', 'best');

% Overlay the validated enclosures, if repro/figure1.sh has been run.
d = fullfile(fileparts(mfilename('fullpath')), '..', 'repro', 'out', 'fig1', 'pendulum');
if isfolder(d)
    figure('Color', 'w'); hold on; grid on
    for k = 1:numel(y0)
        fn = fullfile(d, sprintf('traj_y%g.txt', y0(k)));
        if ~isfile(fn), continue; end
        M = readmatrix(fn, 'NumHeaderLines', 1, 'FileType', 'text');
        xm = (M(:,2) + M(:,3)) / 2;
        ym = (M(:,4) + M(:,5)) / 2;
        plot(xm, ym, col{k}, 'LineWidth', 1.0, ...
             'DisplayName', sprintf('y(0) = %g', y0(k)));
    end
    xlabel('x'); ylabel('y');
    title('validated enclosures (Encoretraj)');
    legend('Location', 'best');
end
