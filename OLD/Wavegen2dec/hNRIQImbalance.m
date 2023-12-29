function [out,ampImbEst,phImbEst] = hNRIQImbalance(varargin)
% hNRIQImbalance I/Q imbalance estimation and compensation
%   [OUT,AMPIMBEST,PHIMBEST] = hNRIQImbalance(IN) returns a time domain IQ
%   baseband waveform OUT compensated with estimates of amplitude imbalance
%   AMPIMBEST (dB) and phase imbalance PHIMBEST (degress). IN is a time
%   domain baseband IQ samples input.
%   [OUT,AMPIMBEST,PHIMBEST] = hNRIQImbalance(...,STEPSIZE) allows an
%   optional input adaption size STEPSIZE, specified as a scalar, that the
%   algorithm uses when estimating the I/Q imbalance.
%   comm.IQImbalanceCompensator Communications Toolbox(TM) System object is
%   used for this purpose.
%
%   See also comm.IQImbalanceCompensator.

% Copyright 2022-2023 The MathWorks, Inc.

    narginchk(1,2);

    in = varargin{1};

    % Default adaptation step size
    stepSize = 1e-5;

    % Override default stepsize if needed
    if nargin == 2
        stepSize = varargin{2};
    end

    % Initialize parameters
    ampImbEst = 0;
    phImbEst = 0;
    out = in;

    % Set the maximum number of iterations needed for the IQ imbalance
    % estimation algorithm to converge. This value is a tradeoff between
    % the estimation accuracy versus processing time and has been
    % determined heuristically.
    maxIterations = 4;

    % Accumulate the IQ imbalance estimates for each iteration
    for idx = 1:maxIterations

        % Obtain amplitude and phase imbalance estimates
        [a,p] = imbalanceEstimation(out,stepSize);
        ampImbEst = ampImbEst+a;
        phImbEst = phImbEst+p;

        % Compensate for IQ imbalance
        out = imbalanceCorrection(out,a,p);

        % Stop here if the accumulated estimates are lower than a
        % threshold. 1e-5 has been heuristically determined.
        if abs(ampImbEst) < 1e-5 && abs(phImbEst) < 1e-5
            break;
        end
    end
end

function [ampImbEst,phImbEst] = imbalanceEstimation(in,stepSize)
%   [AMPIMBEST,PHIMBEST] = imbalanceEstimation(IN,STEPSIZE) estimates the
%   amplitude imbalance AMPIMBEST and phase imbalance PHIMBEST from a
%   signal IN using comm.IQImbalanceCompensator with an adaptive stepsize
%   STEPSIZE

    iqComp = comm.IQImbalanceCompensator('StepSizeSource','Input Port','CoefficientOutputPort',true);

    % Normalize the power of the signal
    in = in./std(in);

    % Get the number of receive antennas R in the waveform
    R = size(in,2);

    % Initialize parameters
    coef = zeros(size(in,1),R);
    ampImbEst = 0;
    phImbEst = 0;

    % Accumulate estimates across receive antennas
    for r = 1:R

        % Reset iqComp for each receive antenna
        if r > 1
            reset(iqComp);
        end
        % Calculate the compensator coefficients
        [~,coefAnt] = iqComp(in(:,r),stepSize);
        coef(:,r) = coefAnt;

        % Compute the amplitude imbalance and phase imbalance for a given compensator coefficient
        [a,p] = iqcoef2imbal(coefAnt(end));

        % Convert to linear before accumulating
        ampImbEst = ampImbEst+db2mag(a);
        phImbEst = phImbEst+p;
    end

    % Average the estimates across number of receive antennas
    ampImbEst = ampImbEst/R;

    % Convert to dB
    ampImbEst = mag2db(ampImbEst);

    % Average across receive antennas
    phImbEst = phImbEst/R;
end

function out = imbalanceCorrection(in,a,p)
%   OUT = imbalanceCorrection(IN,A,P) returns a waveform OUT corrected
%   for a specified amplitude imbalance A and phase imbalance P in the
%   time-domain waveform IN using comm.IQImbalanceCompensator.
%   A is in dB, specified as a real-valued scalar.
%   P is in degrees, specified as a real-valued scalar.

    iqcomp = comm.IQImbalanceCompensator(CoefficientSource="Input port");

    % Obtain compensator coefficient from A and P
    c = iqimbal2coef(a,p);

    % Correct for each receive antenna
    R = size(in,2);
    out = zeros(size(in,1),R);
    for r = 1:R
        out(:,r) = iqcomp(in(:,r),c);
    end
end
