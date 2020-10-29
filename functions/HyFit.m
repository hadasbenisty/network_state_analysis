function [fitResult,Coeff,CoeffCI]=HyFit(x,y)
ft=fittype( '(a*x^b)/(x^b+c^b)+d', 'independent', 'x', 'dependent', 'y');%,'coefficients',{'a','b','c'});
coeffnames(ft);
opts=fitoptions('Method', 'NonlinearLeastSquares');
opts.Lower=[0 1 0.5 0];
opts.Upper=[max(y) 3 30 0.2];
opts.StartPoint=[0.7 2 10 0];
%c is c50, a is rmax, b is the power, and d is baseline
fitResult=fit(x,y,ft,opts);
%fitResult=fit(x,y,ft,'Weight',weights);
Coeff=coeffvalues(fitResult);
CoeffCI=confint(fitResult);
end