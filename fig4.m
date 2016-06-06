%% load contour w from file, and generate contour f as a cubic spline
% w and f do not have to be splines, but they need to be continuse and differntiable

datafile = 'fig4_data';  % fig5_data


load(datafile);

fR2C = @(x) complex(x(:,1), x(:,2));
fC2R = @(x) [real(x) imag(x)];

fGenPP = @(c) cscvn(fC2R(c([1:end 1]))');
fEvalCurve = @(pp) fR2C( ppval(pp, linspace(pp.breaks(1), pp.breaks(end), numel(w)))' );

f = fEvalCurve(fGenPP(cf1));

invGammaDemo(w, f, nv0)
