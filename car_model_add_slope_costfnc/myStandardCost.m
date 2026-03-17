function J = myStandardCost(X, U, e, data, varargin)
% paramsはdata.Parameters{1}から取得。なければローカル定義。
if isfield(data, 'Parameters') && ~isempty(data.Parameters)
    params = data.Parameters{1};
else
    % ローカル定義
    params.Weights.Y   = [1.2 1.2];
    params.Weights.MV  = [0.095 0.095];
    params.Weights.DMV = [0.1 0.1];
    params.Weights.ECR = 1e6;
    params.Weights.Slope = 1e-7; % 傾き誤差の重み（必要に応じて調整）
    params.Scale.Y     = [1 1];
    params.Scale.MV    = [1 1];
    params.Ts          = 0.05;
end
% 再現対象： J = Jy + Ju + Jdu + Jeps   （標準MPCコスト）
% 参照：MathWorks ドキュメント（Linear MPC の標準コスト定義／
% 非線形MPCのカスタムコスト仕様）
% - https://www.mathworks.com/help/mpc/ug/optimization-problem.html
% - https://www.mathworks.com/help/mpc/ug/specify-cost-function-for-nonlinear-mpc.html

% 必要パラメータ（例）
% params.Weights.Y   : [p x ny] または [1 x ny] の行ベクトルを p 回使い回す
% params.Weights.MV  : [p x nMV] または [1 x nMV]
% params.Weights.DMV : [p x nMV] または [1 x nMV]
% params.Weights.ECR : スカラー
% params.Scale.Y     : [1 x ny] 出力スケール（LSyの対角）
% params.Scale.MV    : [1 x nMV] MVスケール（LSuの対角）
%
% 備考：
% - data.References は [p x ny] を想定（nlmpcmove の ref と同じ並び）
% - data.MVTarget   は [p x nMV] を想定（nlmpcmoveopt.MVTarget）
% - 予測出力 Y はモデルの出力関数から生成（非線形MPCでは自前で求める）
% - X(1,:) は現在状態（コスト微分の対象外）、U(end,:) は U(end-1,:) の複製

% 予測ホライズン
p  = data.PredictionHorizon;
ny = data.NumOfOutputs;
nMV = numel(data.MVIndex);

% ========= 1) 出力軌道 Y をモデル出力関数から算出 =========
% 非線形MPCでは Y を X, U から都度生成します（ドキュメント記載）
% data から出力関数ハンドルやパラメータを取得するように適宜変更してください。
% ここでは例として 'myOutputFcn' を使う想定です。
Y = zeros(p+1, ny);
for i = 1:p+1
    % U(i,:) には全入力（MV, MD, UD）が含まれる点に注意
    Y(i,:) = myOutputFcn(X(i,:)', U(i,:)', params)';  %#ok<AGROW>
end

% ========= 2) 参照/ターゲット類の取り出し =========
% 出力参照 r は、ホライズン長 p ぶん（ステップ1..p）を使う
r = data.References;           % [p x ny]
% MV ターゲット ut は、ステップ0..p-1を使う
if isfield(data, 'MVTarget') && ~isempty(data.MVTarget)
    ut = data.MVTarget;        % [p x nMV]
else
    ut = zeros(p, nMV);        % 指定なしならゼロに寄せる
end

% ========= 3) 重みとスケールを p ステップぶんに展開 =========
W_y   = expandWeights(params.Weights.Y  , p, ny);
W_u   = expandWeights(params.Weights.MV , p, nMV);
W_du  = expandWeights(params.Weights.DMV, p, nMV);
ECR   = params.Weights.ECR;

Sy = params.Scale.Y(:)';   % [1 x ny]
Su = params.Scale.MV(:)';  % [1 x nMV]

% ========= 4) 各項の計算 =========

% (a) 出力追従 J_y ： i=1..p を使用（Yは2..p+1が未来出力）
Y_pred = Y(2:p+1, :);        % [p x ny]

Ey = (r - Y_pred) .* Sy;     % スケールして誤差を取る
Jy = sum(sum((W_y .* Ey).^2));

% (a2) 傾き誤差 J_slope
% 車両軌跡の傾き
dx = Y_pred(2:end,1) - Y_pred(1:end-1,1);
dy = Y_pred(2:end,2) - Y_pred(1:end-1,2);
dx(abs(dx)<1e-6) = 1e-6; % 0除算防止
theta_car = dy ./ dx; % [p-1 x 1] 傾きそのもの

% 目標軌跡の傾き
dx_ref = r(2:end,1) - r(1:end-1,1);
dy_ref = r(2:end,2) - r(1:end-1,2);
dx_ref(abs(dx_ref)<1e-6) = 1e-6; % 0除算防止
theta_ref = dy_ref ./ dx_ref; % [p-1 x 1] 傾きそのもの

% 傾き誤差
Etheta = theta_car - theta_ref;
%J_slope = params.Weights.Slope * sum(Etheta.^2); % 二乗誤差
J_slope = params.Weights.Slope * sum(Etheta); % 二乗から一次に変更（必要に応じて調整）

% (b) MV ターゲット追従 J_u ： i=0..p-1 を使用（Uは1..pが最適化変数）
U_mv = U(1:p, data.MVIndex); % [p x nMV]
Eu = (U_mv - ut) .* Su;
Ju = sum(sum((W_u .* Eu).^2));

% (c) MV 変化抑制 J_du ： U(1..p) - U(0..p-1) を MV 列で取り出す
% U(0|k) は data.LastMV（プラント前回出力）を使う
lastMV = data.LastMV(:)';              % [1 x nMV]
U_prev = [lastMV; U(1:p-1, data.MVIndex)];  % [p x nMV]
dU = (U_mv - U_prev) .* Su;
Jdu = sum(sum((W_du .* dU).^2));

% (d) スラック J_eps
Jeps = ECR * (e^2);

% ========= 5) 総和 =========
J = Jy + Ju + Jdu + Jeps + J_slope;

end

% --- ユーティリティ：重みを [p x n] に展開 ---
function W = expandWeights(w, p, n)
    % w が [1 x n] or [p x n] を許可
    if isvector(w) && numel(w)==n
        W = repmat(reshape(w,1,[]), p, 1);
    elseif isequal(size(w), [p n])
        W = w;
    else
        error('Weights must be size [1 x %d] or [%d x %d].', n, p, n);
    end
end