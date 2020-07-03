opts = delimitedTextImportOptions("NumVariables", 8);
opts.DataLines = [9, Inf];
opts.Delimiter = ";";
opts.VariableNames = ["time", "voltage", "current", "power", "temp", "ErrRateCAN1", "ErrRateCAN2", "VarName8"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "string"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, "VarName8", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "VarName8", "EmptyFieldRule", "auto");
data = readtable("data/em.log", opts);

delta_time = data.time;

% Transform to delta time
for i = 1 : size(delta_time, 1) - 1
    delta_time(i, 1) = delta_time(i + 1, 1) - delta_time(i, 1);
end
delta_time(end, 1) = 0;
load_cycle = [delta_time, data.current];
