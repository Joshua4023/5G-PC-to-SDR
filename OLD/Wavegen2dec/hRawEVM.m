function evm = hRawEVM(varargin)
%   hRawEVM Error vector magnitude (EVM) calculation
%   EVM = hRawEVM(...) returns a structure EVM containing error vector magnitude information.
%   EVM is a structure with the fields:
%   RMS  - Root Mean Square (RMS) EVM, the square root of the mean square
%          of the EVM across all input values
%   Peak - Peak EVM, the largest EVM value calculated across all 
%          input values
%   EV   - The normalized error vector
%
%   EVM = hRawEVM(X,R) returns a structure EVM for the input array X given
%   the reference signal array R. The EVM is defined using the error
%   (difference) between the input values X and the reference signal R.
%
%   EVM = hRawEVM(EV) returns a structure EVM for the input array EV which
%   is taken to be the normalized error vector
%   EV=(X-R)/sqrt(mean(abs(R.^2))). This allows for peak and RMS EVM
%   calculation for pre-existing normalized error vectors. This can be used
%   for example to calculate the EVM across an array of previous EVM
%   results, by extracting and concatenating the EV fields from the array
%   to form the EV input.

% Copyright 2020-2022 The MathWorks, Inc.

    % Create a comm.EVM object to measure the EVM by averaging EVM across
    % symbols and layers
    evmComm = comm.EVM(AveragingDimensions=[1 2]);

    % Set number of outputs to 1
    evmComm.MaximumEVMOutputPort = 1;

    if (nargin == 2)
        x = varargin{1};
        r = varargin{2};
        errorVector = x-r;
        p = sqrt(mean(abs(r(:).^2)));
        if (p == 0)
            p = 1;
        end
        evnorm = errorVector/p;
    else
        ev = varargin{1};
        evnorm = ev;
    end
    
    evm.EV = evnorm;
    evmsignal = abs(evnorm(:));

    if nargin == 2
        [rms, peak] = evmComm(r,x);

        % Convert EVM percentages to decimal values
        evm.RMS = rms/100;
        evm.Peak = peak/100;
    else
        evm.Peak = max(evmsignal);
        evm.RMS = sqrt(mean(evmsignal.^2));
    end
end