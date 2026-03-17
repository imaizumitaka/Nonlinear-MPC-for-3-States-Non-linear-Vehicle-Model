%=== figファイルからデータを抽出するサンプル ===

% 1. figファイルを開く（非表示モード）
fig = openfig('default_result2.fig', 'invisible');

% 2. Figure内のAxesオブジェクトを取得
ax = findall(fig, 'Type', 'axes');

% 3. 各AxesからLineオブジェクトを取得してデータ抽出
for i = 1:numel(ax)
    lines = findall(ax(i), 'Type', 'line');
    for j = 1:numel(lines)
        xData = get(lines(j), 'XData');
        yData3 = get(lines(j), 'YData');
        
        % データを表示（必要に応じて保存）
        fprintf('Axes %d, Line %d: %d points\n', i, j, numel(xData));
        disp(table(xData(:), yData3(:), 'VariableNames', {'X', 'Y'}));
    end
end

% 4. 後処理（不要なら閉じる）
close(fig);

save('slope_default.mat','xData','yData3')
