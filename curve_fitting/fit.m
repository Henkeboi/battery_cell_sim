run('soc.m');
run('ocv.m');
x = soc_LUT;
y = flip(ocv_LUT);

constant = lsqcurvefit(@model, [0; 0], x, y);
m = constant(1);
c = constant(2);

xfit = 0 : .0001 : 1;
yfit = model(constant, xfit);
