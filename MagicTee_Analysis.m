function MagicTee_Analysis()
    clc; clear; close all;

    project.title   = 'Magic Tee (Hybrid Tee) Formation and Analysis';
    project.subject = 'Microwave and Waveguide';
    project.team    = {'Chrish Zeyeden - URK23EC1036', ...
                       'Pranov Raghavendhra - UTK23EC1001', ...
                       'Shivabalan - URK23EC1055', ...
                       'Joel Jones - URK23EC1017'};

    cfg.a         = 22.86e-3;
    cfg.b         = 10.16e-3;
    cfg.f0        = 10e9;
    cfg.fSweep    = linspace(8.2e9, 12.4e9, 400);
    cfg.insLossDB = 0.15;
    cfg.isoDB     = 28;
    cfg.rlDB      = 22;

    printHeader(project);

    S_ideal = idealMagicTee();
    S_real  = nonIdealMagicTee(cfg);

    reportMatrix('IDEAL S-MATRIX', S_ideal);
    reportMatrix('NON-IDEAL S-MATRIX', S_real);
    verifyProperties(S_ideal, 'Ideal');
    verifyProperties(S_real,  'Non-ideal');

    cases = defineExcitations();
    P_ideal = powerDistribution(S_ideal, cases);
    P_real  = powerDistribution(S_real,  cases);
    printPowerTable(cases, P_ideal, P_real);

    plotGeometry3D(project, cfg);
    plotSMatrix(S_ideal, S_real);
    plotPowerDistribution(cases, P_ideal, P_real);
    plotSumDifference(S_ideal, S_real);
    plotWaveguideDispersion(cfg);
    plotPortCharacteristics(S_real);
end

function printHeader(p)
    fprintf('\n====================================================\n');
    fprintf('  %s\n', p.title);
    fprintf('  Subject: %s\n', p.subject);
    fprintf('====================================================\n');
    fprintf('  Done by:\n');
    for k = 1:numel(p.team)
        fprintf('    %d. %s\n', k, p.team{k});
    end
    fprintf('====================================================\n\n');
end

function S = idealMagicTee()
    S = (1/sqrt(2)) * [0  0  1  1;
                       0  0  1 -1;
                       1  1  0  0;
                       1 -1  0  0];
end

function S = nonIdealMagicTee(cfg)
    S = idealMagicTee();
    iso = 10^(-cfg.isoDB/20);
    rl  = 10^(-cfg.rlDB/20);
    S(1,2) = iso;  S(2,1) = iso;
    S(3,4) = iso;  S(4,3) = iso;
    S(1,1) = rl;   S(2,2) = rl;
    S(3,3) = rl;   S(4,4) = rl;
    S = S * 10^(-cfg.insLossDB/20);
    sMax = max(svd(S));
    if sMax > 1
        S = S / sMax;
    end
end

function reportMatrix(name, S)
    fprintf('%s\n', name);
    fprintf('%-6s %-18s %-18s %-18s %-18s\n', '', 'Port 1', 'Port 2', 'Port 3 (H)', 'Port 4 (E)');
    for r = 1:4
        fprintf('Port %d ', r);
        for c = 1:4
            fprintf('%-18s ', sprintf('%+.4f%+.4fj', real(S(r,c)), imag(S(r,c))));
        end
        fprintf('\n');
    end
    fprintf('\n');
end

function verifyProperties(S, label)
    symErr  = norm(S - S.', 'fro');
    uniErr  = norm(S' * S - eye(4), 'fro');
    maxSV   = max(svd(S));
    isoDB   = -20*log10(abs(S(1,2)) + eps);
    isoDB34 = -20*log10(abs(S(3,4)) + eps);
    fprintf('PROPERTY CHECK (%s)\n', label);
    fprintf('  Symmetry error ||S - S.''||      : %.3e\n', symErr);
    fprintf('  Unitarity error ||S''S - I||     : %.3e\n', uniErr);
    fprintf('  Largest singular value (<=1)    : %.4f\n', maxSV);
    fprintf('  Isolation Port1-Port2           : %.2f dB\n', isoDB);
    fprintf('  Isolation Port3-Port4           : %.2f dB\n', isoDB34);
    fprintf('  Passive                         : %s\n\n', string(maxSV <= 1 + 1e-9));
end

function cases = defineExcitations()
    cases(1).name = 'H-arm only (a3)';
    cases(1).a    = [0; 0; 1; 0];
    cases(2).name = 'E-arm only (a4)';
    cases(2).a    = [0; 0; 0; 1];
    cases(3).name = 'Port 1 only (a1)';
    cases(3).a    = [1; 0; 0; 0];
    cases(4).name = 'Ports 1+2 in phase';
    cases(4).a    = [1; 1; 0; 0] / sqrt(2);
    cases(5).name = 'Ports 1+2 anti-phase';
    cases(5).a    = [1; -1; 0; 0] / sqrt(2);
end

function P = powerDistribution(S, cases)
    P = zeros(numel(cases), 4);
    for k = 1:numel(cases)
        b = S * cases(k).a;
        P(k,:) = abs(b).'.^2;
    end
end

function printPowerTable(cases, Pi, Pr)
    fprintf('OUTPUT POWER DISTRIBUTION (normalised to unit incident power)\n');
    fprintf('%-24s | %-9s %-9s %-9s %-9s | %-9s\n', 'Excitation', 'P1', 'P2', 'P3', 'P4', 'Total');
    for k = 1:numel(cases)
        fprintf('%-24s | %-9.4f %-9.4f %-9.4f %-9.4f | %-9.4f  [Ideal]\n', cases(k).name, Pi(k,:), sum(Pi(k,:)));
        fprintf('%-24s | %-9.4f %-9.4f %-9.4f %-9.4f | %-9.4f  [Non-ideal]\n', '', Pr(k,:), sum(Pr(k,:)));
    end
    fprintf('\n');
end

function plotGeometry3D(p, cfg)
    a  = cfg.a * 1e3;
    b  = cfg.b * 1e3;
    Lm = 70;
    Lh = 60;
    Le = 60;

    cMain = [0.62 0.78 0.94];
    cH    = [0.98 0.80 0.60];
    cE    = [0.70 0.90 0.70];

    fig = figure('Name', 'Magic Tee 3D Geometry', 'Color', 'w', 'Position', [40 40 1100 720]);
    ax  = axes('Parent', fig, 'Position', [0.05 0.12 0.90 0.80]);
    hold(ax, 'on');

    drawBox(ax, [-Lm Lm],     [-a/2 a/2],        [-b/2 b/2],     cMain, 0.35);
    drawBox(ax, [-a/2 a/2],   [a/2 a/2+Lh],      [-b/2 b/2],     cH,    0.40);
    drawBox(ax, [-a/2 a/2],   [-b/2 b/2],        [b/2 b/2+Le],   cE,    0.40);

    patch(ax, [-a/2 a/2 a/2 -a/2], [a/2 a/2 a/2 a/2], [-b/2 -b/2 b/2 b/2], [1.0 0.55 0.15], 'FaceAlpha', 0.85, 'EdgeColor', 'k');
    patch(ax, [-a/2 a/2 a/2 -a/2], [-b/2 -b/2 b/2 b/2], [b/2 b/2 b/2 b/2], [0.20 0.65 0.30], 'FaceAlpha', 0.85, 'EdgeColor', 'k');

    drawPortFace(ax, 'x', -Lm, a, b, [0.15 0.35 0.75]);
    drawPortFace(ax, 'x',  Lm, a, b, [0.15 0.35 0.75]);
    drawPortFace(ax, 'y',  a/2+Lh, a, b, [0.85 0.40 0.05]);
    drawPortFace(ax, 'z',  b/2+Le, a, b, [0.10 0.50 0.15]);

    text(ax, -Lm-8, 0, 0, 'Port 1', 'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'right', 'Color', [0.15 0.35 0.75]);
    text(ax,  Lm+8, 0, 0, 'Port 2', 'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'left',  'Color', [0.15 0.35 0.75]);
    text(ax, 0, a/2+Lh+8, 0, 'Port 3 (H-arm / Sum)', 'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'center', 'Color', [0.85 0.40 0.05]);
    text(ax, 0, 0, b/2+Le+8, 'Port 4 (E-arm / Difference)', 'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'center', 'Color', [0.10 0.50 0.15]);

    text(ax, 0, a/2+Lh/2, b/2+3, 'H-plane tee (shunt) stub', 'FontAngle', 'italic', 'FontSize', 10, 'HorizontalAlignment', 'center');
    text(ax, a/2+4, 0, b/2+Le/2, 'E-plane tee (series) stub', 'FontAngle', 'italic', 'FontSize', 10, 'HorizontalAlignment', 'left');
    text(ax, -Lm/2, -a/2-4, -b/2-4, 'Collinear arms (Ports 1 and 2)', 'FontAngle', 'italic', 'FontSize', 10, 'HorizontalAlignment', 'center');
    text(ax, 0, 0, -b/2-10, 'Junction', 'FontWeight', 'bold', 'FontSize', 10, 'HorizontalAlignment', 'center');

    text(ax, -Lm+10, -a/2-2, b/2+3, sprintf('a = %.2f mm', a), 'FontSize', 10, 'Color', [0.3 0.3 0.3]);
    text(ax, -Lm+10, -a/2-2, -b/2-4, sprintf('b = %.2f mm', b), 'FontSize', 10, 'Color', [0.3 0.3 0.3]);

    axis(ax, 'equal');
    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, 'x (mm)');
    ylabel(ax, 'y (mm)');
    zlabel(ax, 'z (mm)');
    xlim(ax, [-Lm-45 Lm+45]);
    ylim(ax, [-a/2-20 a/2+Lh+45]);
    zlim(ax, [-b/2-25 b/2+Le+45]);
    view(ax, 38, 24);
    camlight(ax, 'headlight');
    lighting(ax, 'gouraud');
    title(ax, sprintf('%s - 3D Geometry (WR-90)', p.title), 'FontSize', 14, 'FontWeight', 'bold');

    teamText = [{'Done by:'}, p.team];
    annotation(fig, 'textbox', [0.02 0.005 0.55 0.115], 'String', teamText, ...
        'EdgeColor', 'none', 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'none', 'VerticalAlignment', 'bottom');
    annotation(fig, 'textbox', [0.62 0.005 0.36 0.06], 'String', ['Subject: ' p.subject], ...
        'EdgeColor', 'none', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'Interpreter', 'none', 'VerticalAlignment', 'bottom');
    rotate3d(ax, 'on');
end

function drawBox(ax, x, y, z, col, alpha)
    V = [x(1) y(1) z(1);
         x(2) y(1) z(1);
         x(2) y(2) z(1);
         x(1) y(2) z(1);
         x(1) y(1) z(2);
         x(2) y(1) z(2);
         x(2) y(2) z(2);
         x(1) y(2) z(2)];
    F = [1 2 3 4;
         5 6 7 8;
         1 2 6 5;
         2 3 7 6;
         3 4 8 7;
         4 1 5 8];
    patch(ax, 'Vertices', V, 'Faces', F, 'FaceColor', col, 'FaceAlpha', alpha, ...
        'EdgeColor', 'k', 'LineWidth', 1.3);
end

function drawPortFace(ax, axisName, pos, a, b, col)
    switch axisName
        case 'x'
            X = [pos pos pos pos];
            Y = [-a/2 a/2 a/2 -a/2];
            Z = [-b/2 -b/2 b/2 b/2];
        case 'y'
            X = [-a/2 a/2 a/2 -a/2];
            Y = [pos pos pos pos];
            Z = [-b/2 -b/2 b/2 b/2];
        case 'z'
            X = [-a/2 a/2 a/2 -a/2];
            Y = [-b/2 -b/2 b/2 b/2];
            Z = [pos pos pos pos];
    end
    patch(ax, X, Y, Z, col, 'FaceAlpha', 0.9, 'EdgeColor', 'k', 'LineWidth', 1.5);
end

function plotSMatrix(Si, Sr)
    figure('Name', 'S-Matrix Magnitudes', 'Color', 'w', 'Position', [100 80 1000 430]);
    labels = {'1', '2', '3', '4'};
    subplot(1,2,1);
    imagesc(abs(Si)); colormap(parula); caxis([0 1]); colorbar;
    set(gca, 'XTick', 1:4, 'XTickLabel', labels, 'YTick', 1:4, 'YTickLabel', labels);
    title('Ideal |S_{ij}|'); xlabel('Incident port j'); ylabel('Output port i');
    annotateCells(abs(Si));
    subplot(1,2,2);
    imagesc(abs(Sr)); colormap(parula); caxis([0 1]); colorbar;
    set(gca, 'XTick', 1:4, 'XTickLabel', labels, 'YTick', 1:4, 'YTickLabel', labels);
    title('Non-ideal |S_{ij}|'); xlabel('Incident port j'); ylabel('Output port i');
    annotateCells(abs(Sr));
end

function annotateCells(M)
    for r = 1:size(M,1)
        for c = 1:size(M,2)
            text(c, r, sprintf('%.3f', M(r,c)), 'HorizontalAlignment', 'center', 'Color', 'w', 'FontWeight', 'bold');
        end
    end
end

function plotPowerDistribution(cases, Pi, Pr)
    figure('Name', 'Power Distribution', 'Color', 'w', 'Position', [120 100 1100 460]);
    names = {cases.name};
    subplot(1,2,1);
    bar(Pi, 'grouped'); grid on;
    set(gca, 'XTickLabel', names, 'XTickLabelRotation', 20);
    ylabel('Normalised output power'); title('Ideal Magic Tee');
    legend({'P_1', 'P_2', 'P_3', 'P_4'}, 'Location', 'northeastoutside'); ylim([0 1.1]);
    subplot(1,2,2);
    bar(Pr, 'grouped'); grid on;
    set(gca, 'XTickLabel', names, 'XTickLabelRotation', 20);
    ylabel('Normalised output power'); title('Non-ideal Magic Tee');
    legend({'P_1', 'P_2', 'P_3', 'P_4'}, 'Location', 'northeastoutside'); ylim([0 1.1]);
end

function plotSumDifference(Si, Sr)
    phi = linspace(0, 2*pi, 361);
    Pi3 = zeros(size(phi)); Pi4 = Pi3; Pr3 = Pi3; Pr4 = Pi3;
    for k = 1:numel(phi)
        a  = [1; exp(1j*phi(k)); 0; 0] / sqrt(2);
        bi = Si * a;  br = Sr * a;
        Pi3(k) = abs(bi(3))^2;  Pi4(k) = abs(bi(4))^2;
        Pr3(k) = abs(br(3))^2;  Pr4(k) = abs(br(4))^2;
    end
    figure('Name', 'Sum and Difference Operation', 'Color', 'w', 'Position', [140 120 900 480]);
    plot(rad2deg(phi), Pi3, 'b-', 'LineWidth', 2); hold on;
    plot(rad2deg(phi), Pi4, 'r-', 'LineWidth', 2);
    plot(rad2deg(phi), Pr3, 'b--', 'LineWidth', 1.6);
    plot(rad2deg(phi), Pr4, 'r--', 'LineWidth', 1.6);
    grid on; xlim([0 360]); xticks(0:60:360);
    xlabel('Phase difference between inputs at Ports 1 and 2 (degrees)');
    ylabel('Normalised output power');
    title('Sum (Port 3) and Difference (Port 4) Outputs vs Input Phase');
    legend({'P_3 ideal', 'P_4 ideal', 'P_3 non-ideal', 'P_4 non-ideal'}, 'Location', 'best');
end

function plotWaveguideDispersion(cfg)
    c  = 299792458;
    fc = c / (2*cfg.a);
    f  = cfg.fSweep;
    lg = c ./ (f .* sqrt(1 - (fc ./ f).^2));
    vp = c ./ sqrt(1 - (fc ./ f).^2);
    vg = c * sqrt(1 - (fc ./ f).^2);
    figure('Name', 'Waveguide Dispersion (WR-90)', 'Color', 'w', 'Position', [160 140 1000 440]);
    subplot(1,2,1);
    plot(f/1e9, lg*1e3, 'k-', 'LineWidth', 2); grid on;
    xlabel('Frequency (GHz)'); ylabel('\lambda_g (mm)');
    title(sprintf('Guide Wavelength - TE_{10}, f_c = %.3f GHz', fc/1e9));
    hold on; xline(cfg.f0/1e9, 'r--', 'f_0');
    subplot(1,2,2);
    plot(f/1e9, vp/c, 'b-', 'LineWidth', 2); hold on;
    plot(f/1e9, vg/c, 'r-', 'LineWidth', 2); grid on;
    xlabel('Frequency (GHz)'); ylabel('Velocity / c');
    title('Phase and Group Velocity'); legend({'v_p', 'v_g'}, 'Location', 'east');
    fprintf('WAVEGUIDE DATA (WR-90, TE10)\n');
    fprintf('  Cutoff frequency      : %.4f GHz\n', fc/1e9);
    fprintf('  Guide wavelength @f0  : %.3f mm\n', interp1(f, lg, cfg.f0)*1e3);
    fprintf('  Free-space wavelength : %.3f mm\n\n', c/cfg.f0*1e3);
end

function plotPortCharacteristics(S)
    figure('Name', 'S-Parameter Magnitude in dB', 'Color', 'w', 'Position', [180 160 900 450]);
    pairs = {'S_{31}', S(3,1); 'S_{41}', S(4,1); 'S_{32}', S(3,2); 'S_{42}', S(4,2); ...
             'S_{21} (isolation)', S(2,1); 'S_{43} (isolation)', S(4,3); 'S_{11} (return loss)', S(1,1); 'S_{33} (return loss)', S(3,3)};
    vals = cellfun(@(x) 20*log10(abs(x) + eps), pairs(:,2));
    barh(vals, 'FaceColor', [0.2 0.45 0.75]); grid on;
    set(gca, 'YTick', 1:size(pairs,1), 'YTickLabel', pairs(:,1), 'YDir', 'reverse');
    xlabel('Magnitude (dB)'); title('Non-ideal Magic Tee - Key S-Parameters');
    for k = 1:numel(vals)
        text(vals(k) - 1, k, sprintf('%.1f dB', vals(k)), 'HorizontalAlignment', 'right', 'Color', 'w', 'FontWeight', 'bold');
    end
end
