function H = hChannelEstimateEVM3GPP(rxGrid,refGrid,varargin)
%hChannelEstimateEVM3GPP Channel estimation customized for 3GPP EVM mode
%   H = hChannelEstimateEVM3GPP(...) Estimates the channel returning the
%   channel coefficients H. H is a K-by-N-by-R-by-P array where K is the
%   number of subcarriers, N is the number of symbols, and R is the number
%   of receive antennas and P is the number of reference signal ports
%
%   H = hChannelEstimateEVM3GPP(RXGRID,REFGRID) Performs the channel
%   estimation customized for 3GPP EVM mode on the received resource grid
%   RXGRID. RXGRID must be an array of size K-by-N-by-R. REFGRID is a
%   predefined reference array with nonzero elements representing the
%   reference symbols in their appropriate locations. It is of size
%   K-by-N-by-P.  REFGRID can span multiple slots. K is the number of
%   subcarriers, N is the number of OFDM symbols, R is the number of
%   receive antennas and P is the number of reference signal ports.
%
%   H = hChannelEstimateEVM3GPP(...,FREQAVGMETHOD) FREQAVGMETHOD is an
%   optional parameter to select smoothening of the channel coefficients in
%   the frequency direction, using a moving average filter. Smoothening is
%   performed as described in TS 38.104 Annex B.6 (FR1) or C.6 (FR2)
%   When this parameter is not specified, the time averaging across the
%   duration of REFGRID is enabled.
%
%   H = hChannelEstimateEVM3GPP(...,NAME,VALUE) specifies
%   additional options as NAME,VALUE pairs:
%
%   'CDMLengths'      - A 2-element row vector [FD TD] specifying the
%                       length of frequency domain code domain multiplexing
%                       (FD-CDM) and time domain code domain multiplexing
%                       (TD-CDM) despreading. A value of 1 for an element
%                       indicates no CDM and a value greater than 1
%                       indicates the length of the CDM. For example, [2 1]
%                       indicates FD-CDM2 and no TD-CDM. The default is [1
%                       1] (no orthogonal despreading)

% Copyright 2020-2023 The MathWorks, Inc.

    firstOptArg = 1;
    if nargin == 3 || nargin == 5
        freqAvgMethod = varargin{1};
    else
        freqAvgMethod = [];
    end
    if ~isempty(freqAvgMethod)
        firstOptArg = firstOptArg + 1;
    end
    fcnName = 'hChannelEstimateEVM3GPP';
    optNames = {'CDMLengths'};
    opts = coder.const(nr5g.internal.parseOptions( ...
            fcnName,optNames,varargin{firstOptArg:end}));

    fdCDM = double(opts.CDMLengths(1));
    tdCDM = double(opts.CDMLengths(2));

    % Channel estimation is performed as defined in TS 38.104 Annex
    % B.6(FR1)/C.6(FR2). The pilots are averaged in time and frequency and
    % the resulting vector is linearly interpolated.
    [K,N,P] = size(refGrid);

    % Extract reference indices and symbols from reference grid
    refInd = find(refGrid(:)~=0);
    refSym = refGrid(refInd);

    % Get the channel estimate output dimensions. The final channel
    % estimate will be of size K-by-N-by-R-by-P
    R = size(rxGrid,3);

    % Create the channel estimate grid
    H = complex(zeros([K N R P]));

    % ---------------------------------------------------------------------
    % LS estimation
    % ---------------------------------------------------------------------

    % For each transmit port
    for p = 1:P
        % Get frequency (subcarrier k) and time (OFDM symbol n) subscripts
        % of reference signal for the current port. 'thisport' is a logical
        % indexing vector for the current port, which extracts the
        % corresponding reference symbols
        [ksubs,nsubs,thisport] = getPortSubscripts(K,N,P,refInd,p);
        refSymThisPort = refSym(thisport);

        % For each OFDM symbol
        un = unique(nsubs).';
        for uni = 1:numel(un)

            % Get frequency and OFDM symbol subscripts
            n = un(uni);
            k = ksubs(nsubs==n);
            % For each receive antenna
            for r = 1:R
                % Perform least squares (LS) estimate of channel in the
                % locations of the reference symbols. 'H_LS' is a column
                % vector containing the LS estimates for all subcarriers
                % for the current port, OFDM symbol, and receive antenna
                H_LS = rxGrid(k,n,r) ./ refSymThisPort(nsubs==n);
                % Perform FD-CDM despreading if required
                if (fdCDM>1)
                    % 'm' is zero if the LS estimates divide evenly into
                    % FD-CDM groups, otherwise 'm' is the number of LS
                    % estimates in the final partial group
                    m = mod(numel(H_LS),fdCDM);
                    for a = 0:double(m~=0)
                        if (~a)
                            % whole CDM lengths (may be empty)
                            k_LS = 1:(numel(H_LS)-m);
                            nkCDM = fdCDM;
                        else
                            % part CDM length (may be empty)
                            k_LS = numel(H_LS) + (-m+1:0);
                            nkCDM = m;
                        end
                        % Extract the LS estimates and reshape so that each
                        % column contains an FD-CDM group
                        x = reshape(H_LS(k_LS),nkCDM,[]);
                        % Average across first dimension (i.e. across the
                        % FD-CDM group) and repeat the averaged value
                        x = repmat(mean(x,1),[nkCDM 1]);
                        % Reshape back into a single column
                        H_LS(k_LS) = reshape(x,[],1);
                    end
                end

                % Assign the estimates into the appropriate region of
                % the overall channel estimate array
                H(k,n,r,p) = H_LS;
            end
        end
    end

    % Perform TD-CDM despreading if required
    if (tdCDM>1)
        % 'm' is zero if the estimates divide evenly into TD-CDM
        % groups, otherwise 'm' is the number of estimates in the final
        % partial group
        m = mod(numel(un),tdCDM);
        for a = 0:double(m~=0)
            if (~a)
                % whole CDM lengths (may be empty)
                unCDM = un(1:end-m);
                nlCDM = tdCDM;
            else
                % part CDM length (may be empty)
                unCDM = un(end-m+1:end);
                nlCDM = m;
            end
            for r = 1:R
                % Extract the estimates and reshape so that each row
                % and plane contains a TD-CDM group
                x = reshape(H(:,unCDM,r,p),K,nlCDM,[]);
                % 'mx' is the mean of 'x' across the 2nd dimension
                % (i.e. across the TD-CDM group), excluding zeros
                mx = sum(x,2) ./ max(sum(x~=0,2),1);
                % Repeat the averaged value and reshape back into a
                % matrix
                x = repmat(mx,[1 nlCDM]);
                H(:,unCDM,r,p) = reshape(x,[],numel(unCDM));
            end
        end
    end

    % Get phase (theta) and magnitude(radius) of complex estimates
    [theta,radius] = cart2pol(real(H),imag(H));

    for p = 1:P
        for r = 1:R
            % Perform averaging in time and frequency for phase and
            % magnitude separately and recombine these seperate components
            
            % Declare vector to temporarily store averaged phase and then
            % magnitude
            phasemagVec = zeros(size(H,1),2);
            for phaseOrMag = 1:2
                % Define temporary H vector and set estimates to phase then
                % magnitude. Initialize an averaging vector used as a
                % placeholder for both time and frequency averaging
                % directions
                avgVec = complex(zeros(size(H,1),1));
                if phaseOrMag==1
                    H_EST_tmp = theta(:,:,r,p);
                else
                    H_EST_tmp = radius(:,:,r,p);
                end

                if phaseOrMag == 1
                    % First obtain all valid columns
                    col = find(any(abs(H_EST_tmp),1));
                    for colIdx = 1:length(col)
                        currentcol = col(colIdx);
                        colInd = find(H_EST_tmp(:,currentcol));
                        colPhase = unwrap(H_EST_tmp(colInd,currentcol));
                        H_EST_tmp(colInd,currentcol) = colPhase;
                    end
                end
                % First obtain all valid rows
                row = find(sum(H_EST_tmp,2));
                for rowIdx = 1:length(row)
                    currentRow = row(rowIdx);
                    avgVec(currentRow)= (sum(H_EST_tmp(currentRow,:)))./length(find(H_EST_tmp(currentRow,:)));
                end
                vec = find(avgVec);

                % Remove subcarriers with no pilot symbols and then
                % perform frequency averaging
                avgVec = avgVec(vec);

                if strcmp(freqAvgMethod,'movingAvgFilter')
                    % 'isEven' indicates whether the length of 'avgVec' is
                    % odd (0) or even (1). It is used to extract the
                    % correct elements after filtering and to correctly
                    % calculate the normalization factor
                    isEven = (mod(length(avgVec),2)==0);

                    % Perform averaging in frequency direction using a
                    % moving window of size N. N = 19 as defined in TS
                    % 38.104 Annex B.6(FR1) or C.6(FR2). At the edges,
                    % where the number of samples are less than N, the
                    % window size is reduced to span 1, 3, 5, 7 ... samples
                    N = 19;

                    % Use filter to perform part of the averaging (not
                    % normalized)
                    zeroPad = zeros(N-1,1);
                    data = [avgVec; zeroPad];
                    weights = ones(N,1);
                    freqAvg = filter(weights,1,data);

                    % 'isEvenOrLong' configures the appropriate odd or even
                    % symmetry in the elements of 'freqAvg' extracted after
                    % filtering. 'isEvenOrLong' is additionally used in
                    % the calculation of the normalization factor
                    isEvenOrLong = (isEven || length(avgVec)>=N);

                    % Remove unwanted elements
                    if (isEvenOrLong)
                        X = N;
                    else
                        X = length(data);
                    end
                    removeIdx = [(2:2:X) (length(data) - ((X-N+2):2:N) + isEvenOrLong)];
                    freqAvg(removeIdx) = [];

                    % Normalization factor. As stated above and according
                    % to the standard, only odd reference subcarriers 1, 3,
                    % 5, 7 ... are kept when averaging less than the full
                    % 19 reference subcarriers. 'M' here determines the
                    % (1-based) index of the last reference subcarrier kept
                    % (i.e. reference subcarriers 1:2:M are kept); for
                    % length(avgVec)>=N, M = N-2 = 17, as this is the last
                    % odd reference subcarrier index <19. For
                    % length(avgVec)<N, M = length(avgVec)-isEven (M is the
                    % last odd index whether length(avgVec) is odd or
                    % even).
                    M = min(length(avgVec)-isEven,N-2);
                    normFactor = [1:2:M N*ones(1,length(freqAvg)-(N-1)) (M-(1-isEvenOrLong)*2):-2:1]';
                    freqAvg = freqAvg./normFactor;

                    % Place frequency averaged pilots into temporary storage
                    % vector
                    if ~isempty(freqAvg)
                        phasemagVec(vec,phaseOrMag) = freqAvg;
                    end
                else
                    if ~isempty(avgVec)
                        phasemagVec(vec,phaseOrMag) = avgVec;
                    end
                end
            end

            % Convert averaged symbols back to Cartesian coordinates
            [X,Y] = pol2cart(phasemagVec(:,1),phasemagVec(:,2));
            H_tmp = complex(X,Y);

            % Interpolate coefficients to account for missing pilots
            if length(avgVec) > 1
                interpEqCoeff = interp1(find(H_tmp~=0),H_tmp(H_tmp~=0),(1:length(H_tmp)).','linear','extrap');
            else
                interpEqCoeff = H_tmp;
                if sum(abs(H_tmp))
                    interpEqCoeff(:,1) = H_tmp(H_tmp~=0);
                end
            end
            H(:,:,r,p) = repmat(interpEqCoeff,1,size(H,2));
        end
    end
end

function [ksubs,nsubs,thisport] = getPortSubscripts(K,N,P,ind,port)
    % Extract a frequency-time grid (along with the associated port) from
    % the given linear indices
    [ksubs,nsubs,psubs] = ind2sub([K N P],ind(:));
    thisport = (psubs==port);
    ksubs = ksubs(thisport);
    nsubs = nsubs(thisport);
end