function values = trimodalGaussian(N, M, mu, sigma)
% values = trimodalGaussian(N, M, mu, sigma)
% Returns an N-by-M matrix containing pseudorandom values drawn from the
% trimodal Gaussian distribution, with parameters as set by the user.
% 
% Input:
%       mu    : Array of length 3, each value representing the mean for one of
%               the three modes of the distribution.
%       sigma : Scalar representing the standard deviation of the modes.
%               Each mode has the same standard deviation.

errorThreshold = 0.0001;

pdf       = @(x) (normpdf(x, mu(1), sigma) + normpdf(x, mu(2), sigma) + normpdf(x, mu(3), sigma)) ./3;
range     = [min(mu)-4.*sigma max(mu)+4.*sigma 100000];
t         = linspace(range(1), range(2), range(3));

% Check pdf
validate(integral(pdf, range(1), range(2)), errorThreshold, 1);

pdfSample = pdf(t);
cdfSample = cumsum(pdfSample) .* ((range(2) - range(1)) ./ range(3));

% Check cdf
validate(cdfSample(end), errorThreshold, 2);

% Ensure monotonically-increasing requirement of interp1
[cdfUniqueSample, idx] = unique(cdfSample); 

data      = rand(N, M);
values    = interp1(cdfUniqueSample, t(idx), data, 'linear');
end

function validate(value, errorThreshold, mode)
    if(abs(1 - value) > errorThreshold)
            if mode == 1
                warning('Integral over the domain of the generated pdf deviates significantly from expected 1.0000')
        elseif mode == 2
                warning('Final value of generated CDF deviates significantly fom expected 1.0000')
            end
    end
end