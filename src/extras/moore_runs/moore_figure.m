% moore_figure.m
%
% Figure for R. Moore's classical rotation example (the standard
% wrapping-effect benchmark), tested with our IVP enclosure solver
% (ivp2026, EndCover mode).
%
%   x' = y,  y' = -x,   initial box: x in [0.95,1.05], y in [-0.05,0.05].
%
% The flow is an exact rigid rotation (Jacobian [[0,1],[-1,0]] is
% skew-symmetric, log-norm = 0): the TRUE end set at time T is just the
% initial square rotated by angle -T about the origin -- same size and
% shape, forever.
%
% Snapshots taken every pi/4 from T=0 to T=2*pi (one full period, 9 points;
% T=2*pi coincides with T=0, confirming exact periodicity / no drift).
%
% Output 1 (moore_overall.png): ONE plot with all 9 snapshots overlaid --
%   the small square orbiting the unit circle, our eps-cover boxes drawn
%   at each position (looks solid at this scale -- that's the point).
% Output 2 (moore_zoom_grid.png): a 3x3 grid, one zoomed panel per
%   snapshot, each showing the individual eps-boxes tiling the true
%   rotated-square outline (black) at that instant.
%
% Data: moore_runs/<label>/E1.txt produced by ivp2026 (EndCover mode),
%   src/examples/eg10Moore.mk  (n=2, ff = y,-x, cen=(1,0), wid=(0.05,0.05))

clear; clc; close all;

cen0 = [1.0, 0.0];
w0   = [0.05, 0.05];   % half-widths
corners0 = [cen0(1)-w0(1), cen0(2)-w0(2);
            cen0(1)+w0(1), cen0(2)-w0(2);
            cen0(1)+w0(1), cen0(2)+w0(2);
            cen0(1)-w0(1), cen0(2)+w0(2)];

snapLabels = {'T0','Tpi4','Tpi2','T3pi4','Tpi','T5pi4','T3pi2','T7pi4','T2pi'};
snapNames  = {'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'};
snapAngles = [0, pi/4, pi/2, 3*pi/4, pi, 5*pi/4, 3*pi/2, 7*pi/4, 2*pi];
nSnap = numel(snapLabels);
cmap = lines(nSnap);

% Load all snapshot box sets once.
snapBoxes = cell(nSnap,1);
for k = 1:nSnap
    dname = fullfile('.', snapLabels{k});
    boxes = read_boxes_file(fullfile(dname,'E1.txt'));
    N = min(numel(boxes), 85); % keep only genuine cover boxes (E1 boxes= count)
    snapBoxes{k} = boxes(1:N);
end

% ════════════════════════════════════════════════════════════════════
% Output 1: overall orbit figure (all 9 snapshots together)
% ════════════════════════════════════════════════════════════════════
fig1 = figure('Color','w','Position',[100 100 800 800]);
hold on; axis equal; box on;

th = linspace(0,2*pi,400);
plot(cos(th), sin(th), 'Color',[0.7 0.7 0.7], 'LineStyle',':','LineWidth',1);

for k = 1:nSnap
    boxes = snapBoxes{k};
    for i = 1:numel(boxes)
        b = boxes{i};
        x1=b(1); x2=b(2); y1=b(3); y2=b(4);
        rectangle('Position',[x1,y1,x2-x1,y2-y1], ...
            'FaceColor', [cmap(k,:) 0.45], 'EdgeColor', cmap(k,:), 'LineWidth', 0.5);
    end
    ang = -snapAngles(k); % flow rotates by -t
    R = [cos(ang) -sin(ang); sin(ang) cos(ang)];
    rc = (R*corners0')';
    plot([rc(:,1); rc(1,1)], [rc(:,2); rc(1,2)], '-', 'Color', cmap(k,:)*0.6, 'LineWidth', 1.4);
    ctr = mean(rc,1);
    % T=0 and T=2*pi coincide exactly (confirms periodicity) -- offset
    % their labels vertically so both remain legible.
    if k == 1
        yoff = 0.16;
    elseif k == nSnap
        yoff = -0.16;
    else
        yoff = 0.10;
    end
    text(ctr(1), ctr(2)+yoff, ['T=' snapNames{k}], 'FontSize', 8, ...
        'FontWeight','bold', 'HorizontalAlignment','center', ...
        'BackgroundColor','w', 'Margin',0.5);
end
xlabel('x'); ylabel('y');
xlim([-1.3 1.3]); ylim([-1.3 1.3]);
title({'Moore''s rotation example: x''=y, y''=-x', ...
    '(x_0,y_0)\in[0.95,1.05]\times[-0.05,0.05], eps=0.02', ...
    'our eps-cover at T = 0,\pi/4,...,2\pi -- orbits with no growth'});

set(fig1,'PaperPositionMode','auto');
print(fig1, 'moore_overall', '-dpng', '-r300');
savefig(fig1, 'moore_overall.fig');
fprintf('Saved moore_overall.png and .fig\n');

% ════════════════════════════════════════════════════════════════════
% Output 2: 3x3 grid, zoomed internal structure per snapshot
% ════════════════════════════════════════════════════════════════════
fig2 = figure('Color','w','Position',[100 100 1400 1400]);
for k = 1:nSnap
    subplot(3,3,k); hold on; axis equal; box on;
    boxes = snapBoxes{k};
    for i = 1:numel(boxes)
        b = boxes{i};
        x1=b(1); x2=b(2); y1=b(3); y2=b(4);
        rectangle('Position',[x1,y1,x2-x1,y2-y1], ...
            'FaceColor',[0.85 0.92 1], 'EdgeColor',[0 0.3 0.75], 'LineWidth', 0.6);
    end
    ang = -snapAngles(k);
    R = [cos(ang) -sin(ang); sin(ang) cos(ang)];
    rc = (R*corners0')';
    plot([rc(:,1); rc(1,1)], [rc(:,2); rc(1,2)], 'k-', 'LineWidth', 1.6);

    ctr = mean(rc,1);
    half = norm(2*w0)/2*1.15;
    xlim([ctr(1)-half, ctr(1)+half]); ylim([ctr(2)-half, ctr(2)+half]);
    xlabel('x'); ylabel('y');
    title(sprintf('T = %s   (%d boxes)', snapNames{k}, numel(boxes)));
end
sgtitle({'Moore''s rotation example -- internal eps-cover structure at each snapshot', ...
    'each panel: our eps-boxes (blue) hugging the true rotated square (black outline)'});

set(fig2,'PaperPositionMode','auto');
print(fig2, 'moore_zoom_grid', '-dpng', '-r300');
savefig(fig2, 'moore_zoom_grid.fig');
fprintf('Saved moore_zoom_grid.png and .fig\n');

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
