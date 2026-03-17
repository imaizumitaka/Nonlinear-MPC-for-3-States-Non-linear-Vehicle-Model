function y = myOutputFcn(x, u, params)
% x: 状態ベクトル [x; y; theta]
% u: 入力ベクトル [v; delta]
% params: パラメータ（未使用）

% 出力は [x, y] のみ
% paramsはmyStandardCostのインターフェースのために受け取るが、ここでは使わない

y = [x(1); x(2)];
end
