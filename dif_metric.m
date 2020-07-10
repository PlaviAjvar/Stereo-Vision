function dif = dif_metric(d, GT)

delta = double(d) - double(GT);
entries = ~isnan(delta);
dif = sum(abs(delta(entries))) / sum(entries(:));

end