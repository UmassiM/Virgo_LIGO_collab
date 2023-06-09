function [fitresult, gof] = createFit(h0, crvec)
%CREATEFIT(H0,CRVEC)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : h0
%      Y Output: crvec
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 20-Jul-2021 18:30:46


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( h0, crvec );

% Set up fittype and options.
ft = 'linearinterp';
excludedPoints = excludedata( xData, yData, 'Indices', [19 20 21 22 23 24 25] );
opts = fitoptions( 'Method', 'LinearInterpolant' );
opts.Normalize = 'on';
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData, excludedPoints );
legend( h, 'crvec vs. h0', 'Excluded crvec vs. h0', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'h0', 'Interpreter', 'none' );
ylabel( 'crvec', 'Interpreter', 'none' );
grid on


