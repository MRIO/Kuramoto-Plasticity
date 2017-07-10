function values = bimodalGaussian(N,M,mu, sigma)
threshold = 0.0001;
trim = 0.0001;

pdf = @(x) (normpdf(x,mu(1),sigma) + normpdf(x,mu(2),sigma))./2;

range = [min(mu)-4.*sigma max(mu)+4.*sigma 100000];
t = linspace(range(1),range(2),range(3));

% Check pdf
validate(integral(pdf,range(1),range(2)), threshold, 1);

pdfSample = pdf(t);
cdfSample = cumsum(pdfSample).*((range(2)-range(1))./range(3));

% Check cdf
validate(cdfSample(end),threshold,2);

idx = (cdfSample<trim|cdfSample>(1-trim)); % Trim edges of cdf to ensure monotonicity, required for interp1
cdfSample(idx) = [];
t(idx) = [];

data = rand(N,M);
values = interp1(cdfSample,t,data,'linear');
end

function validate(value, threshold, mode)
if(abs(1-value) > threshold)
    if mode== 1
        warning('Integral over the domain of the generated pdf deviates significantly from expected 1.0000')
    elseif mode == 2
        warning('Final value of generated CDF deviates significantly fom expected 1.0000')
    end
end
end