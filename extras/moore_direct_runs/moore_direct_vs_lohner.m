% moore_direct_vs_lohner.m
%
% Honest head-to-head: the "pure direct method" (stepB=3, CAPD's
% C0HORect2Set -- plain rigorous Taylor + interval remainder, re-enclosed
% in an axis-aligned box every step, NO Lohner/doubleton wrapping control)
% vs our actual solver pipeline (stepB=0, CR-Lohner doubleton via
% CnRect2Set, run through EndCover as in moore_runs/).
%
% IMPORTANT CONTEXT (see conversation): ivp2026's EndCover/Boundary loop
% unconditionally intersects every step with a fresh single-shot
% stepBcrlohner(...) call ("Tighten end-enclosure with a global one-step
% bound", boundary.h L57-59 / ivp2026.cpp), REGARDLESS of the stepB the
% user picks. That is why setting stepB=3 in eg10Moore.mk gave the exact
% same near-machine-precision result as stepB=0: the direct method's
% real (looser) result was being silently overridden by that hardcoded
% Lohner cross-check every step.
%
% To see the genuine, uncorrected direct-method behaviour, a standalone
% program (moore_direct_test.cpp, compiled via `make moore_direct_test`)
% calls puredirectmethod() directly in a fixed-step loop (h = pi/4),
% completely bypassing that intersection. This is real CAPD output (not
% a hand-simulated model): with a 45 degree rotation per step, the
% axis-aligned re-enclosure grows by a factor of exactly sqrt(2) every
% step -- confirmed in the data: wmax = 0.1, 0.1414, 0.2, 0.283, 0.4,
% 0.566, 0.8, 1.131, 1.6 (one period, 8 steps), then 25.6 after 2
% periods, 409.6 after 3 periods (0.1*sqrt(2)^n exactly).
%
% Data:
%   moore_direct_runs/direct_wmax.txt   -- t wmax xlo xhi ylo yhi  (direct method, real CAPD data)
%   moore_runs/<label>/E1.txt           -- our Lohner-based EndCover eps-cover (real CAPD data)

clear; clc; close all;

%% Load direct-method data (t, wmax, box bounds) -- genuine wrapping effect
D = readmatrix(fullfile('.', 'direct_wmax.txt'));
tD = D(:,1); wD = D(:,2);
period = 2*pi;

%% Load our Lohner-EndCover data at the 9 snapshots (T=0,pi/4,...,2pi)
lohnerDir  = fullfile('..','moore_runs');
snapLabels = {'T0','Tpi4','Tpi2','T3pi4','Tpi','T5pi4','T3pi2','T7pi4','T2pi'};
snapAngles = [0, pi/4, pi/2, 3*pi/4, pi, 5*pi/4, 3*pi/2, 7*pi/4, 2*pi];
nSnap = numel(snapLabels);

cen0 = [1.0, 0.0]; w0 = [0.05, 0.05];
corners0 = [cen0(1)-w0(1), cen0(2)-w0(2);
            cen0(1)+w0(1), cen0(2)-w0(2);
            cen0(1)+w0(1), cen0(2)+w0(2);
            cen0(1)-w0(1), cen0(2)+w0(2)];

wLohner = zeros(nSnap,1);
lohnerBoxes = cell(nSnap,1);
for k = 1:nSnap
    boxes = read_boxes_file(fullfile(lohnerDir, snapLabels{k}, 'E1.txt'));
    N = min(numel(boxes), 85);
    boxes = boxes(1:N);
    lohnerBoxes{k} = boxes;
    lo = [inf inf]; hi = [-inf -inf];
    for i = 1:N
        b = boxes{i};
        lo = min(lo, [b(1) b(3)]);
        hi = max(hi, [b(2) b(4)]);
    end
    wLohner(k) = max(hi - lo);
end

% ════════════════════════════════════════════════════════════════════
% Figure 1: wmax vs elapsed periods (log scale) -- the headline result
% ════════════════════════════════════════════════════════════════════
fig1 = figure('Color','w','Position',[100 100 900 650]);
hold on; box on;
semilogy(tD/period, wD, 'r.-', 'LineWidth', 1.6, 'MarkerSize', 10);
semilogy(snapAngles/period, wLohner, 'bo-', 'LineWidth', 1.8, ...
    'MarkerFaceColor','b', 'MarkerSize', 6);
yline(norm(2*w0)/sqrt(2)*sqrt(2), 'k--', 'LineWidth', 1);
xlabel('elapsed time (periods = t / 2\pi)');
ylabel('enclosure diameter  wmax  (log scale)');
legend({'stepB=3: pure direct method (real wrapping effect)', ...
        'stepB=0: our Lohner-based EndCover pipeline', ...
        'true diameter (constant, = 0.1 or 0.1\surd2)'}, ...
        'Location','northwest', 'FontSize', 9);
title({'Moore''s rotation example: x''=y, y''=-x', ...
       'direct method (stepB=3, no doubleton) vs our Lohner pipeline (stepB=0)', ...
       'both computed directly with CAPD -- real data, not simulated'});
set(gca,'YScale','log');

set(fig1,'PaperPositionMode','auto');
print(fig1, 'moore_direct_vs_lohner_wmax', '-dpng', '-r300');
savefig(fig1, 'moore_direct_vs_lohner_wmax.fig');
fprintf('Saved moore_direct_vs_lohner_wmax.png and .fig\n');

% ════════════════════════════════════════════════════════════════════
% Figure 2: boxes at the 9 snapshots within ONE period, same axes
% (direct method's real growth vs Lohner's real tightness, drawn to the
%  SAME scale so the contrast is visible directly)
% ════════════════════════════════════════════════════════════════════
fig2 = figure('Color','w','Position',[100 100 1500 750]);

% -- panel (a): direct method boxes (one period, 9 snapshots) --
subplot(1,2,1); hold on; axis equal; box on;
cmap = lines(nSnap);
idxOnePeriod = 1:9; % rows in direct_wmax.txt for t=0..2*pi (h=pi/4 run)
for k = idxOnePeriod
    xlo=D(k,3); xhi=D(k,4); ylo=D(k,5); yhi=D(k,6);
    rectangle('Position',[xlo,ylo,xhi-xlo,yhi-ylo], ...
        'FaceColor',[cmap(k,:) 0.25], 'EdgeColor', cmap(k,:), 'LineWidth', 1.2);
end
axis equal; xl = xlim; yl = ylim;
xlabel('x'); ylabel('y');
title({'(a) stepB=3 pure direct method:', 'real wrapping effect over ONE period', ...
    sprintf('(box grows %.1fx per period)', wD(9)/wD(1))});

% -- panel (b): Lohner-EndCover boxes, SAME axis limits as (a) --
subplot(1,2,2); hold on; axis equal; box on;
th = linspace(0,2*pi,400);
plot(cos(th), sin(th), 'Color',[0.7 0.7 0.7], 'LineStyle',':','LineWidth',1);
for k = 1:nSnap
    boxes = lohnerBoxes{k};
    for i = 1:numel(boxes)
        b = boxes{i};
        rectangle('Position',[b(1),b(3),b(2)-b(1),b(4)-b(3)], ...
            'FaceColor',[cmap(k,:) 0.45], 'EdgeColor', cmap(k,:), 'LineWidth', 0.5);
    end
end
xlim(xl); ylim(yl); % SAME scale as panel (a) -- Lohner boxes barely visible at this scale
xlabel('x'); ylabel('y');
title({'(b) stepB=0 our Lohner-based EndCover:', 'same axis limits as (a) --', ...
       'true size is imperceptible at this scale'});

sgtitle('Moore''s rotation example -- same physical scale, drawn from real CAPD output');

set(fig2,'PaperPositionMode','auto');
print(fig2, 'moore_direct_vs_lohner_boxes', '-dpng', '-r300');
savefig(fig2, 'moore_direct_vs_lohner_boxes.fig');
fprintf('Saved moore_direct_vs_lohner_boxes.png and .fig\n');

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
