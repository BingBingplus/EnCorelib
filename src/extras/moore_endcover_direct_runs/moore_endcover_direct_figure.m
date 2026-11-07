% moore_endcover_direct_figure.m
%
% Our EndCover method (the real subdivide-until-wmax<=eps algorithm,
% endcover::EndCover + Extendnew/Refinenew), but with the per-box
% propagation done by the PURE DIRECT METHOD (stepB=3, CAPD's
% C0HORect2Set) -- and with the hardcoded CR-Lohner "global tightening"
% intersection that ivp2026.cpp normally applies every step REMOVED, so
% stepB=3 genuinely drives the result (see moore_endcover_direct.cpp).
%
%   x' = y,  y' = -x,  initial box [0.95,1.05]x[-0.05,0.05], eps = 0.02.
%
% Because the direct method has no wrapping control, EndCover's
% subdivide-until-eps loop must work harder as elapsed time grows: the
% NUMBER OF BOXES needed to reach the same eps = 0.02 is different at
% every snapshot (85 at T=0,pi/4; 341 at pi/2..pi; 1365 at 5pi/4,3pi/2;
% 5461 at 7pi/4; ... at 2pi) -- this is exactly what one would expect
% from a method with no wrapping mitigation: it is still SOUND (each
% box is still a valid enclosure, epsilon-small), it just needs
% increasingly many boxes to stay epsilon-small as time (and hence
% accumulated wrapping) grows.
%
% Data: E0_<label>.txt / E1_<label>.txt written by moore_endcover_direct.exe
%   (built via `make moore_endcover_direct`), one pair per snapshot.

clear; clc; close all;

cen0 = [1.0, 0.0]; w0 = [0.05, 0.05];
corners0 = [cen0(1)-w0(1), cen0(2)-w0(2);
            cen0(1)+w0(1), cen0(2)-w0(2);
            cen0(1)+w0(1), cen0(2)+w0(2);
            cen0(1)-w0(1), cen0(2)+w0(2)];

snapLabels = {'T0','Tpi4','Tpi2','T3pi4','Tpi','T5pi4','T3pi2','T7pi4','T2pi'};
snapNames  = {'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'};
snapAngles = [0, pi/4, pi/2, 3*pi/4, pi, 5*pi/4, 3*pi/2, 7*pi/4, 2*pi];
nSnap = numel(snapLabels);
cmap = lines(nSnap);

snapBoxes = cell(nSnap,1);
snapCount = zeros(nSnap,1);
for k = 1:nSnap
    boxes = read_boxes_file(fullfile('.', ['E1_' snapLabels{k} '.txt']));
    snapBoxes{k} = boxes;
    snapCount(k) = numel(boxes);
end

% ════════════════════════════════════════════════════════════════════
% Figure 1: overall orbit -- all 9 snapshots, box count annotated
% ════════════════════════════════════════════════════════════════════
fig1 = figure('Color','w','Position',[100 100 900 900]);
hold on; axis equal; box on;

th = linspace(0,2*pi,400);
plot(cos(th), sin(th), 'Color',[0.7 0.7 0.7], 'LineStyle',':','LineWidth',1);

for k = 1:nSnap
    boxes = snapBoxes{k};
    for i = 1:numel(boxes)
        b = boxes{i};
        rectangle('Position',[b(1),b(3),b(2)-b(1),b(4)-b(3)], ...
            'FaceColor', [cmap(k,:) 0.45], 'EdgeColor', cmap(k,:), 'LineWidth', 0.4);
    end
    ang = -snapAngles(k);
    R = [cos(ang) -sin(ang); sin(ang) cos(ang)];
    rc = (R*corners0')';
    ctr = mean(rc,1);
    if k == 1, yoff = 0.16; elseif k == nSnap, yoff = -0.16; else, yoff = 0.11; end
    text(ctr(1), ctr(2)+yoff, sprintf('T=%s\n(%d boxes)', snapNames{k}, snapCount(k)), ...
        'FontSize', 8, 'FontWeight','bold', 'HorizontalAlignment','center', ...
        'BackgroundColor','w', 'Margin', 0.5);
end
xlabel('x'); ylabel('y');
xlim([-1.3 1.3]); ylim([-1.3 1.3]);
title({'Our EndCover method, stepped by the pure direct method (stepB=3)', ...
       'x''=y, y''=-x,  [0.95,1.05]\times[-0.05,0.05],  eps=0.02', ...
       'box count needed to reach eps grows with elapsed time'});

set(fig1,'PaperPositionMode','auto');
print(fig1, 'moore_endcover_direct_overall', '-dpng', '-r300');
savefig(fig1, 'moore_endcover_direct_overall.fig');
fprintf('Saved moore_endcover_direct_overall.png and .fig\n');

% ════════════════════════════════════════════════════════════════════
% Figure 2: 3x3 grid, zoomed internal structure per snapshot
% ════════════════════════════════════════════════════════════════════
fig2 = figure('Color','w','Position',[100 100 1400 1400]);
for k = 1:nSnap
    subplot(3,3,k); hold on; axis equal; box on;
    boxes = snapBoxes{k};
    for i = 1:numel(boxes)
        b = boxes{i};
        rectangle('Position',[b(1),b(3),b(2)-b(1),b(4)-b(3)], ...
            'FaceColor',[0.85 0.92 1], 'EdgeColor',[0 0.3 0.75], 'LineWidth', 0.4);
    end
    ang = -snapAngles(k);
    R = [cos(ang) -sin(ang); sin(ang) cos(ang)];
    rc = (R*corners0')';
    plot([rc(:,1); rc(1,1)], [rc(:,2); rc(1,2)], 'k-', 'LineWidth', 1.6);

    ctr = mean(rc,1);
    lo = [inf inf]; hi = [-inf -inf];
    for i = 1:numel(boxes)
        b = boxes{i};
        lo = min(lo, [b(1) b(3)]); hi = max(hi, [b(2) b(4)]);
    end
    half = max(hi-lo)/2*1.15;
    xlim([ctr(1)-half, ctr(1)+half]); ylim([ctr(2)-half, ctr(2)+half]);
    xlabel('x'); ylabel('y');
    title(sprintf('T = %s   (%d boxes)', snapNames{k}, snapCount(k)));
end
sgtitle({'Our EndCover method with pure direct-method stepping (stepB=3, no doubleton)', ...
    'blue: eps-cover boxes (eps=0.02)   black: true rotated square outline'});

set(fig2,'PaperPositionMode','auto');
print(fig2, 'moore_endcover_direct_zoom_grid', '-dpng', '-r300');
savefig(fig2, 'moore_endcover_direct_zoom_grid.fig');
fprintf('Saved moore_endcover_direct_zoom_grid.png and .fig\n');

% ════════════════════════════════════════════════════════════════════
% Figure 3: box count needed vs elapsed time -- the headline number
% ════════════════════════════════════════════════════════════════════
fig3 = figure('Color','w','Position',[100 100 800 600]);
semilogy(snapAngles/(2*pi), snapCount, 'o-', 'Color',[0 0.3 0.75], ...
    'LineWidth', 1.8, 'MarkerFaceColor',[0 0.3 0.75], 'MarkerSize', 7);
for k = 1:nSnap
    text(snapAngles(k)/(2*pi), snapCount(k)*1.25, num2str(snapCount(k)), ...
        'HorizontalAlignment','center', 'FontSize', 8);
end
xlabel('elapsed time (periods = T / 2\pi)');
ylabel('number of eps-cover boxes needed (log scale)');
set(gca,'YScale','log');
box on; grid on;
title({'Our EndCover method, direct-method stepping (stepB=3):', ...
       'box count needed to reach eps=0.02 grows with elapsed time'});

set(fig3,'PaperPositionMode','auto');
print(fig3, 'moore_endcover_direct_boxcount', '-dpng', '-r300');
savefig(fig3, 'moore_endcover_direct_boxcount.fig');
fprintf('Saved moore_endcover_direct_boxcount.png and .fig\n');

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
