clear all;
clc;

nx = 3; % x, y, theta
ny = 2; % x, y
nu = 2; % v, delta
nlobj = nlmpc(nx, ny, nu);

Ts = 0.05; % 50ms sampling
nlobj.Ts = Ts;
nlobj.PredictionHorizon = 20;
nlobj.ControlHorizon = 5;

nlobj.Model.StateFcn = @(xk, uk, params) carDT(xk, uk, params); % give discrete car model
nlobj.Model.IsContinuousTime = false;
nlobj.Model.NumberOfParameters = 1;
nlobj.Model.OutputFcn = @(x,u,Ts) [x(1); x(2)]; % x, y


% --- myStandardCost用パラメータ設定 ---
params.Weights.Y   = [1 1];   % 出力(x, y)の重み
params.Weights.MV  = [0.1 0.1]; % 入力(v, delta)の重み
params.Weights.DMV = [0.1 0.1]; % 入力変化の重み
params.Weights.ECR = 1e6;     % スラック変数の重み
params.Scale.Y     = [1 1];   % 出力スケール
params.Scale.MV    = [1 1];   % 入力スケール
params.Ts          = Ts;      % サンプリング周期もparamsに含める

% initial values for x, u
x0 = [0; 0; 0];
u0 = [0; 0];
validateFcns(nlobj, x0, u0, [], {params});

nlobj.Optimization.CustomCostFcn = 'myStandardCost';

% Kalman filter to estimate car state in executing
EKF = extendedKalmanFilter(@carStateFcn,@carMeasurementFcn);
x = x0;
y = x0(1:2);
EKF.State = x; % initial value for EKF
mv = [0; 0]; % = u (v; delta)

nloptions = nlmpcmoveopt; % object of options for MPC
nloptions.Parameters = {params};

load('xypath.mat'); % reference car (x,y) path

Duration = 10; % simulation time
tLength = Duration/Ts;
xHistory = zeros(length(x(:,1)), tLength); % x state result
xHistory(:, 1) = x0;

for t = 1:tLength
    t
    % Correct previous prediction
    xk = correct(EKF,y);
    % Compute optimal control moves
    % --- 参照軌跡をPredictionHorizon分生成 ---
    if t+nlobj.PredictionHorizon-1 <= size(xypath,1)
        yrefH = xypath(t:t+nlobj.PredictionHorizon-1,:);
    else
        % 末尾は最後の点で埋める
        yrefH = [xypath(t:end,:); repmat(xypath(end,:), nlobj.PredictionHorizon-(size(xypath,1)-t+1), 1)];
    end
    [mv,nloptions] = nlmpcmove(nlobj,xk,mv,yrefH,[],nloptions);
    % Predict prediction model states for the next iteration
    predict(EKF,[mv; Ts], params);
    % Implement first optimal control move
    x = carDT(x, mv, params);
    % Generate sensor data (with noise)
    y = x([1 2]) + randn(2,1)*0.01;
    % Save plant states
    xHistory(:, t + 1) = x;
end

save ("xHistory", "xHistory");

%%

load("xHistory.mat");
figure;
plot(xypath(:,1), xypath(:,2));
hold on;
plot(xHistory(1,:), xHistory(2,:));
legend("reference", "actual");
title("Vehicle Control by Non-linear MPC");
xlabel("x");
ylabel("y");
hold off;

%%
% 軌跡追従の偏差を計算・プロット
% xHistory: シミュレーションで得られた車両状態 [x; y; theta]
% xypath: 参照軌跡 [x, y]

deviation = zeros(1, tLength); % 偏差格納用
for t = 1:tLength
    car_pos = xHistory(1:2, t)'; % [x, y]
    ref_pos = xypath(t, :);     % [x_ref, y_ref]
    deviation(t) = norm(car_pos - ref_pos); % ユークリッド距離
end

figure;
plot((0:tLength-1)*Ts, deviation, 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Deviation [m]');
title('軌跡追従の偏差');
grid on;
