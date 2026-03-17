function xk1 = carStateFcn(xk, u, params)
% paramsはdata.Parameters{1}から取得してください（EKF用）
if nargin < 3
	params.Ts = 0.05; % デフォルト値
end
uk = u(1:2);
xk1 = carDT(xk, uk, params);