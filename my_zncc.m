function m = my_zncc(w1, w2)
	w1 = w1 - mean(w1(:));
	w2 = w2 - mean(w2(:));
	denom = sqrt( sum(sum(w1.^2))*sum(sum(w2.^2)) );
	m = sum(sum((w1.*w2))) / denom;
end