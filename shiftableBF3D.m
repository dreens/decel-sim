function [ f , param ]  =  shiftableBF3D(f0, sigmas, sigmar, tol, varargin)
% make image real so complex arithmetic is supported:
f0 = double(f0);

% User can specify dynamic range directly as the last argument if they have
% an alternative estimate for this.
if ~isempty(varargin)
    T = varargin{1};
else
    T = approxDynamicRange(f0);
end

assert(sigmar~=0,'Range Kernel of zero not solvable')
assert(sigmas~=0,'Spatial Kernel of zero will return original image')
assert(tol<1&&tol>=0,'Tolerance between 0 and 1, smaller more accurate')

% set range interval and the order of raised cosine
N  =  ceil( 0.405 * (T / sigmar)^2 );
gamma    =  1 / (sqrt(N) * sigmar);
twoN     =  2^N;

% compute truncation
if tol == 0
    M = 0;
else
    if sigmar > 40
        M = 0;
    elseif sigmar > 10
        sumCoeffs = 0;
        for k = 0 : round(N/2)
            sumCoeffs = sumCoeffs + nchoosek(N,k)/twoN;
            if sumCoeffs > tol/2
                M = k;
                break;
            end
        end
    else
        M = ceil( 0.5*( N - sqrt(4*N*log(2/tol)) ) );
    end
end
% filter
h = waitbar(0, 'Computing filter ...');
warning('off'); %#ok<WNOFF>
ss   =  size(f0);
fnum  =  zeros(ss);
fdenom  =  zeros(ss);
f  =  zeros(ss);
ii = sqrt(-1);
if N < 50
    for k = M : N - M
        waitbar( (k - M + 1) / (N - 2*M + 1), h);
        omegak = (2*k - N)*gamma;
        bk = nchoosek(N,k) / twoN;
        H  = exp(-ii*omegak*f0);
        G  = conj(H);
        F  = G.*f0;
        barF = imgaussfilt3(real(F), sigmas)...
             + imgaussfilt3(imag(F), sigmas)*1i;
        barG = imgaussfilt3(real(G), sigmas)...
             + imgaussfilt3(imag(G), sigmas)*1i;
        fnum =  fnum + bk * H .* barF;
        fdenom  = fdenom + bk * H .* barG;
    end
    close(h);
else
    for k = M : N - M
        waitbar((k - M + 1) / (N - 2*M + 1), h);
        omegak = (2*k - N)*gamma;
        % use Sterling's approximation
        bk = exp(logfactorial(N) - logfactorial(k) ...
            - logfactorial(N-k) - N*log(2));
        H  = exp(-ii*omegak*f0);
        G  = conj(H);
        F  = G.*f0;
        barF = imgaussfilt3(real(F), sigmas)...
             + imgaussfilt3(imag(F), sigmas)*1i;
        barG = imgaussfilt3(real(G), sigmas)...
             + imgaussfilt3(imag(G), sigmas)*1i;
        fnum =  fnum + bk * H .* barF;
        fdenom  = fdenom + bk * H .* barG;
    end
    close(h);
end
% check: avoid division by zero
idx1 = find( fdenom < 1e-3);
idx2 = find( fdenom > 1e-3);
f( idx1 ) = f0( idx1 );
f( idx2 ) = real(fnum( idx2 )./fdenom (idx2 ));
% save parameters
param.T  = T;
param.N  = N;
param.M  = M;
end

% Stirling's Approximation to the log of the factorial, with the small n
% results done normally. (Otherwise n=0 diverges)
function val = logfactorial(n)
   val = n*(log(n)-1) + 0.5*log(2*pi*n);
   if n<10
       val = log(factorial(n));
   end
end

% Instead of looping over every possible f(x-y) within the gaussian range
% kernel, I find the largest adjacent pixel difference and multiply it by
% the range kernel width. If this ends up being larger than the max
% possible dynamic range, I return the latter.
function dynR = approxDynamicRange(f0)
    mx3 = @(x,d) abs(max(max(max(diff(x,d)))));
    dynR = max([mx3(f0,1) mx3(f0,2) mx3(f0,3)]) * sigmaR;
    
    dynR = min( dynR, max(max(max(f0))) - min(min(min(f0)))  );
end
