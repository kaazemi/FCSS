function baseline = calc_base(sys)
% calculates baseline of the traces
[n,T] = size(sys.y);
baseline = min(sys.y,[],2);
sigma = sqrt(diag(sys.R));
p = length(sigma);
if isequal(n,p)
    y = sys.y-repmat(baseline,1,T);
    for i = 1:n
        a = y(i,:);
        a = a(a<3*sigma(i));
    baseline(i) = mean(a);
    end
else
    disp('Warning: number of measurements (n) is less than the ambient dimension (p)')
    disp('Baseline set to minimum of y !')
end

end

