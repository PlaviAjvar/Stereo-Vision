function m = my_zncc(w1, w2)
	w1 = w1 - mean(w1(:));
	w2 = w2 - mean(w2(:));
	denom = sqrt( sum(sum(w1.^2))*sum(sum(w2.^2)) );

	if denom < 1e-10
        if sum >= 0
            m = Inf;
        else
            m = -Inf;
        end
	else
		m = sum(sum((w1.*w2))) / denom;
    end
end