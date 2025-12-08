% Sample data
%X = randn(100, 3);          % 100 observations, 3 predictors
%beta_true = [2; -1; 0.5];   % true coefficients (excluding intercept)
%y = 1 + X * beta_true + randn(100,1);  % intercept = 1, added Gaussian noise

DI = [0.71383013	1.076801997	2.179470155	2.475532632	1.390707877	1.336111172	1.882954204	2.457886055	0.564935512	1.542106938	0.835725102	0.87349666];
GG = [4.229739602	3.091520556	2.359538225	1.132795818	2.204975712	2.739220425	1.868538631	1.69935869	3.452308202	1.995905641	3.012676539	2.596460037];
RE = [0.4632, 0.3783, 0.3029, 0.2619, 0.3579, 0.3140, 0.2635, 0.2323, 0.1787, 0.1160, 0.1889, 0.1512];
FF = [6.67, 6.67, 6.67, 6.67, 6.67, 6.67, 6.67, 6.67, 3.33, 3.33, 3.33, 3.33];

%coef_names = {'$D_i^R(0)$', 'Re', 'F'};

% Compute regression uncertainty
[beta_mean, beta_std] = regression_uncertainty([RE; FF]', GG');

disp('Estimated coefficients (mean):');
disp(beta_mean);

disp('Standard deviation (uncertainty):');
disp(beta_std);
