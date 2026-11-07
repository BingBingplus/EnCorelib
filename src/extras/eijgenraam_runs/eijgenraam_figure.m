% eijgenraam_figure.m
%
% Figure for Eijgenraam's wrapping-effect example (Ph.D. thesis, Example
% 6.4.2, pp. 127-128; reproduced in Nedialkov's Ph.D. thesis), tested
% with our IVP enclosure solver (ivp2026, EndCover mode).
%
%   x' = 0,  y' = x^2,   initial set:  x in [-lambda,lambda], y = 0.
%   Exact solution set at time T:  y = T*x^2,  x in [-lambda,lambda].
%
% Panel (a): the single tightest axis-aligned convex box enclosing the
%            end set -- illustrates Eijgenraam's point that ANY single
%            convex enclosure of a curved set has unavoidable
%            (quadratic-order) excess area.
% Panels (b)-(d): our EndCover eps-cover of the end set for
%            eps = 0.3, 0.1, 0.02 -- many small boxes, each intersecting
%            the true curve, hugging the parabola tightly and refining
%            as eps -> 0.  This is NOT subject to Eijgenraam's obstruction
%            because we never enclose the whole curved set by one box.
%
% Data: eijgenraam_runs/eps<VAL>/E1.txt produced by
%   ivp2026 (build/examples/version_2026), mode=0 (EndCover),
%   src/examples/eg9Eijgenraam.mk  (n=2, ff = 0, x^2, cen=(0,0), wid=(1,0))
%
% The first N lines of each E1.txt are the actual eps-cover boxes
% (N = the "E1 boxes=" count printed by the solver, see cover_count.txt);
% the remaining lines are auxiliary E0/corner-tracking boxes appended by
% the solver for visualisation and are NOT part of the cover.

clear; clc; close all;

lambda = 1.0;
T      = 1.0;

epsList   = [0.3, 0.1, 0.02];
epsLabels = {'eps0.3','eps0.1','eps0.02'};

% Analytic solution set: y = T*x^2, x in [-lambda,lambda]
xx = linspace(-lambda, lambda, 2000);
yy = T * xx.^2;

fig = figure('Color','w','Position',[100 100 1500 420]);

% ── Panel (a): naive single convex (box) enclosure ─────────────────────
subplot(1,4,1); hold on; axis equal; box on;
ybox = [0, T*lambda^2];
hBox = patch([-lambda lambda lambda -lambda], [ybox(1) ybox(1) ybox(2) ybox(2)], ...
    [1 0.85 0.85], 'EdgeColor',[0.7 0 0], 'LineWidth',1.8, 'DisplayName','excess region (no trajectory points)');
hCurve = plot(xx, yy, 'k-', 'LineWidth', 2.2, 'DisplayName','true end set  y = Tx^2');
xlim([-1.15*lambda, 1.15*lambda]); ylim([-0.15*T*lambda^2, 1.25*T*lambda^2]);
xlabel('x'); ylabel('y');
title({'(a) Single convex enclosure','(Eijgenraam''s obstruction)'});
legend([hBox, hCurve], 'Location','southoutside','FontSize',7);

% ── Panels (b)-(d): our eps-cover (EndCover) ────────────────────────────
for k = 1:numel(epsList)
    epsv = epsList(k);
    dname = fullfile('.', epsLabels{k});
    boxes = read_boxes_file(fullfile(dname,'E1.txt'));
    N = read_count(fullfile(dname,'cover_count.txt'));
    if isnan(N) || N > numel(boxes)
        N = numel(boxes);
    end
    cover = boxes(1:N);

    subplot(1,4,1+k); hold on; axis equal; box on;
    for i = 1:numel(cover)
        b = cover{i};
        x1=b(1); x2=b(2); y1=b(3); y2=b(4);
        rectangle('Position',[x1,y1,x2-x1,y2-y1], ...
            'FaceColor',[0.85 0.92 1], 'EdgeColor',[0 0.3 0.75], 'LineWidth',0.6);
    end
    plot(xx, yy, 'k-', 'LineWidth', 1.6);
    xlim([-1.15*lambda, 1.15*lambda]); ylim([-0.15*T*lambda^2, 1.25*T*lambda^2]);
    xlabel('x'); ylabel('y');
    title(sprintf('(%c) Our eps-cover, eps=%.3g\n(%d boxes)', char('a'+k), epsv, N));
end

sgtitle(sprintf(['Eijgenraam''s example:  x''=0, y''=x^2,  x_0\\in[-%.0f,%.0f], T=%.1f\n' ...
    'IVP enclosure solver (ivp2026, EndCover mode)'], lambda, lambda, T));

set(fig,'PaperPositionMode','auto');
print(fig, 'eijgenraam_figure', '-dpng', '-r300');
savefig(fig, 'eijgenraam_figure.fig');
fprintf('Saved eijgenraam_figure.png and .fig\n');

% ── helpers ──────────────────────────────────────────────────────────
function boxes = read_boxes_file(filename)
    fid = fopen(filename, 'r');
    if fid == -1
        boxes = {};
        return;
    end
    boxes = {};
    tline = fgetl(fid);
    while ischar(tline)
        tokens = regexp(tline, ...
            'Box\s+\d+:\s+\[([+\-\d\.eE]+),\s*([+\-\d\.eE]+)\]\s*x\s*\[([+\-\d\.eE]+),\s*([+\-\d\.eE]+)\]', ...
            'tokens');
        if ~isempty(tokens)
            nums = str2double(tokens{1});
            if numel(nums) == 4 && all(isfinite(nums))
                boxes{end+1} = nums; %#ok<AGROW>
            end
        end
        tline = fgetl(fid);
    end
    fclose(fid);
end

function n = read_count(filename)
    fid = fopen(filename,'r');
    if fid == -1
        n = NaN;
        return;
    end
    n = str2double(strtrim(fgetl(fid)));
    fclose(fid);
end
