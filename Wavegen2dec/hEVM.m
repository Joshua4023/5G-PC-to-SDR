function evmInfo = hEVM(varargin)
%hEVM EVM calculation
%   EVMINFO = hEVM(CARRIER,EQGRID,REFGRID)
%   Calculates the error vector magnitude (EVM) of a slot-based grid
%   containing equalized symbols EQGRID, using the grid of reference
%   modulated symbols REFGRID and returns the EVM statistics in structure
%   EVMINFO. 
%   3GPP specified EVM grid selection is optionally performed when the
%   fourth dimension of input EQGRID is 2
%   EQGRID and REFGRID are 4-dimensional matrices of size S-by-N-by-L-by-E.
%   S is the number of subcarriers, given by CARRIER.NSizeGrid * 12. N is
%   the number of symbols spanning the duration of the waveform.  L is the
%   number of transmitted layers. E is the number of EVM window positions
%
%   EVMINFO is a struct containing EVM statistics with fields.
%      SubcarrierRMS   - Root mean square (RMS) EVM per subcarrier
%                        (Column vector of N subcarriers)
%      SubcarrierPeak  - Peak EVM per subcarrier
%                        (Column vector of N subcarriers)
%      SymbolRMS       - RMS EVM per symbol
%                        (Column vector of S symbols)
%      SymbolPeak      - Peak EVM per symbol
%                        (Column vector of S symbols)
%      SlotRMS         - RMS EVM per slot
%                        (Column vector of L slots)
%      SlotPeak        - Peak EVM per slot
%                        (Column vector of L slots)
%      EVM             - Struct array of EVM statistics. Each array element
%                        contains EVM per PDSCH slot and per window
%                        location, peak and RMS EVM for the overall
%                        waveform). When EVM is not measured according to
%                        3GPP standard, only one window location is used.
%                        (E x Number of slots)
%      EVMGRID         - Error vector matrix for the selected window edge
%                        EVM is averaged across layers where EVMGRID is a
%                        (2-D array of dimensions N subcarriers x symbols)
%      OverallEVM      - Structure containing EVM statistics for the
%                        overall waveform. It contains these fields:
%       EV             - Error vector for the overall waveform
%                        (Array of J x L). J is of number of all equalized
%                        samples in the waveform.
%       Peak           - Peak EVM for the overall waveform
%                        (Scalar)
%       RMS            - RMS EVM for the overall waveform
%                        (Scalar)
%
%   CARRIER            - Carrier configuration object, <a
%                        href="matlab:help('nrCarrierConfig')"
%                        >nrCarrierConfig</a>
%
%   EVMINFO = hEVM(NRB,SCS,EQGRID,REFGRID,NAME,VALUE)
%   Specifies an additional option as a NAME,VALUE pair to set cyclic
%   prefix
%   'CyclicPrefix'     - Cyclic prefix ('normal' (default), 'extended')
%   NRB                - Number of RBs in carrier resource grid
%   SCS                - Subcarrier spacing in kHz
%
%   EVMINFO = hEVM(CARRIER,EQ,REF,IND)
%   Calculates the error vector magnitude (EVM) of equalized IQs EQ 
%   using reference IQs REF whose locations are given by IND
%   EQ                 - Column vector of size N. N is the length of EQ
%   REF                - Column vector of size N. N is the length of REF
%   IND                - Column vector of 1-based linear indices for a slot
%
%   EVMINFO = hEVM(NRB,SCS,EQ,REF,IND,NAME,VALUE)
%   Specifies an additional option as a NAME,VALUE pair to set cyclic
%   prefix
%   'CyclicPrefix'     - cyclic prefix ('normal' (default), 'extended')
%
%   NRB                - Number of RBs in carrier resource grid
%   SCS                - Subcarrier spacing in kHz
%
%   % Example:
%   % This example illustrates the various possible signatures for 'hEVM'
%
%   rng('default');
%
%   % Create a column vector of random reference IQs
%   refIQs = 1-3i*rand(624,1);
% 
%   % Copy refIQs to eqIQs . Add noise to eqIQs
%   eqIQs = awgn(refIQs,30);
% 
%   % Create a column vector of indices
%   ind = 1:624;
% 
%   % Prepare a carrier object
%   carrier = nrCarrierConfig;
%   nrb = carrier.NSizeGrid;
%   scs = carrier.SubcarrierSpacing;
% 
%   % Calculate EVM using a carrier based syntax along with a column vector for indices
%   out = hEVM(carrier,eqIQs,refIQs,ind);
%
%   % Calculate EVM using slot-based indices and without a carrier based syntax 
%   out = hEVM(nrb,scs,eqIQs,refIQs,ind);
% 
%   % Prepare two grids. Populate the first grid with reference IQs with 'refIQs'. Fill the second grid with 'eqIQs'.
%   carrier = nrCarrierConfig('SubcarrierSpacing',60,'CyclicPrefix','extended');
%   refGrid = nrResourceGrid(carrier);
%   eqGrid = refGrid;
%   refGrid(ind) = refIQs;
%   eqGrid(ind) = eqIQs;
% 
%   % Calculate EVM using a carrier based syntax with reference and equalized IQ grids
%   out = hEVM(carrier,refGrid,eqGrid,'CyclicPrefix','extended');
% 
%   % Calculate EVM using reference IQ grid, equalized IQ grid and additional information needed to internally form a reference grid
%   out = hEVM(nrb,scs,refGrid,eqGrid,'CyclicPrefix','extended');

%   Copyright 2020-2022 The MathWorks, Inc.

    narginchk(3,7);
    fcnName = 'hEVM';

    % Determine if carrier syntax is used
    isCarrierSyntax = isa(varargin{1},'nrCarrierConfig');
    if (isCarrierSyntax)
        carrier = varargin{1};
        rxSlotGrid = varargin{2};
        refSlotGrid = varargin{3};
        if nargin == 4
            ind = varargin{4};
        end
        ofdmInfo = nrOFDMInfo(carrier);
    else
        % Parse for optional name-value pair (cyclicPrefix)
        nrb = varargin{1};
        scs = varargin{2};
        rxSlotGrid = varargin{3};
        refSlotGrid = varargin{4};
        optNames = {'CyclicPrefix'};
        firstOptArg = 5;
        cp = 'normal';

        % Check for presence of the optional parameter 'ind'
        if nargin >= firstOptArg && isnumeric(varargin{firstOptArg})
            ind = varargin{firstOptArg};
            firstOptArg = firstOptArg + 1;
        end
        if nargin == (firstOptArg+1) && (ischar(varargin{firstOptArg}) || isstring(varargin{firstOptArg}))
             opts = coder.const(nr5g.internal.parseOptions( ...
            fcnName,optNames,varargin{firstOptArg:end}));
            cp = opts.CyclicPrefix;
        end
        ofdmInfo = nrOFDMInfo(nrb,scs,'CyclicPrefix',cp);
    end

    % If input is a column vector, construct a slot-grid
    if iscolumn(rxSlotGrid)
        % Construct grids, rxSlot & refSlot
        if isCarrierSyntax
            nrb = carrier.NSizeGrid;
        end
        numSCs = nrb*12;
        rxSlot = zeros(numSCs,ofdmInfo.SymbolsPerSlot,1,1);
        refSlot = zeros(numSCs,ofdmInfo.SymbolsPerSlot,1,1);
        rxSlot(ind) = rxSlotGrid;
        refSlot(ind) = refSlotGrid;
        rxSlotGrid = rxSlot;
        refSlotGrid = refSlot;

        % Initialize grid parameters
        nEdge = 1;
        nSlots = 1;
        nLayers = 1;
        evmGridEdge = zeros(1,numSCs,ofdmInfo.SymbolsPerSlot);
    else
        numSCs = size(rxSlotGrid,1);
        numSym = size(rxSlotGrid,2);
        nSlots = floor(size(rxSlotGrid,2)/ofdmInfo.SymbolsPerSlot);
        if nSlots == 0
            nSlots = 1;
        end
        % Sum across subcarriers, layers, and evm window dimensions to
        % locate the range of allocated symbols
        nEdge = size(rxSlotGrid,4);
        nLayers = size(rxSlotGrid,3);
        evmGridEdge = zeros(2,numSCs,numSym);
    end

    % Locate allocated symbols in grid
    symRange = find(sum(rxSlotGrid,[1 3 4]))-1;
    slotRange = unique(floor(symRange/ofdmInfo.SymbolsPerSlot));
    if iscolumn(slotRange)
        slotRange = slotRange.';
    end

    % evm3GPP is used to choose the grid with higher EVM, when 2 grids are
    % present with EVM statistics
    evm3GPP = nEdge==2;
    evm = repmat(hRawEVM([]), 1+evm3GPP , nSlots);
    nFrames = floor(nSlots/(10*ofdmInfo.SlotsPerSubframe));
    frameEVM = repmat(hRawEVM([]), 1, max(nFrames,1));
    L = ofdmInfo.SymbolsPerSlot;

    % Loop over each valid slot with these steps
    % - Extract allocated REs, reference IQs per layer, per edge
    % - Calculate raw error vector
    % - Construct up to 2 EVM grids corresponding to the allocated slots
    for slotIdx = slotRange
        for e = 1:nEdge
            symIdx = (slotIdx)*L+(1:L);
            rxSymb = rxSlotGrid(:,symIdx,:,e);
            refSymb = refSlotGrid(:,symIdx,:,e);
            rxSymbols = reshape(rxSymb,numSCs*ofdmInfo.SymbolsPerSlot,nLayers);
            refSymbols = reshape(refSymb,numSCs*ofdmInfo.SymbolsPerSlot,nLayers);
            ind  = find(rxSymbols);
            ind = reshape(ind,length(ind)/nLayers,nLayers);
            rxSymbols = rxSymbols(ind);
            refSymbols = refSymbols(ind);
            evm(e,slotIdx+1) = hRawEVM(rxSymbols,refSymbols);
        end
        % Build low and high edge EVM grids across all slots
        evmSlotEdge = zeros(2,numSCs,ofdmInfo.SymbolsPerSlot);
        % In case of non-3GPP case, only a single EVM grid is used
        numSymbPerSlot = ofdmInfo.SymbolsPerSlot;
        evmSlotEdge(1,ind(:,1)) = mean(abs(evm(1,slotIdx+1).EV)*100,2);
        evmGridEdge(1,:,numSymbPerSlot*slotIdx+1:numSymbPerSlot*(slotIdx+1)) =  evmSlotEdge(1,:,:);
        if evm3GPP
            evmSlotEdge(2,ind(:,1)) = mean(abs(evm(2,slotIdx+1).EV)*100,2);
            evmGridEdge(2,:,numSymbPerSlot*slotIdx+1:numSymbPerSlot*(slotIdx+1)) =  evmSlotEdge(2,:,:);
        end
    end
	% The low or high edge timing is chosen for plotting
	% automatically based on whichever has the largest RMS across
	% all slots (Largest RMS chosen as mentioned in TS 38.101-1/2 , Annex F.6)
	evmGrid = squeeze(evmGridEdge(1,:,:));
    if evm3GPP
        evmMaxLow = max([evm(1,:).RMS]);
        evmMaxHigh = max([evm(2,:).RMS]);
        if evmMaxHigh > evmMaxLow
            evmGrid = squeeze(evmGridEdge(2,:,:));
        end
    end
    % RMS and Peak EVM versus subcarrier
	evmSubcarrierRMS = sqrt(sum(evmGrid.^2,2)./sum(evmGrid~=0,2));
	evmSubcarrierPeak = max(evmGrid,[],2)./any(evmGrid,2);

    % RMS and Peak EVM versus OFDM symbol
	evmSymbolRMS = sqrt(sum(evmGrid.^2,1)./sum(evmGrid~=0,1)).';
	evmSymbolPeak = (max(evmGrid,[],1)./any(evmGrid,1)).';

    % RMS and Peak EVM versus slot
	evmSlotRMS = [];
	evmSlotPeak = [];
    for sIdx = 1:nSlots
        evmSlotRMSTmp = max(evm(1:nEdge,sIdx).RMS)*100;
        evmSlotPeakTmp = max(evm(1:nEdge,sIdx).Peak)*100;
        if isempty(evmSlotPeakTmp)
            evmSlotPeakTmp = NaN;
        end
        evmSlotRMS= [evmSlotRMS;evmSlotRMSTmp]; %#ok<*AGROW>
        evmSlotPeak = [evmSlotPeak;evmSlotPeakTmp];
    end

    evmInfo.SubcarrierRMS = evmSubcarrierRMS;
    evmInfo.SubcarrierPeak = evmSubcarrierPeak;
    evmInfo.SymbolRMS = evmSymbolRMS;
    evmInfo.SymbolPeak = evmSymbolPeak;
    evmInfo.SlotRMS = evmSlotRMS;
    evmInfo.SlotPeak = evmSlotPeak;
    evmInfo.EVM = evm;
    evmInfo.EVMGrid = evmGrid;

    % Below loop needs to run at least once
    if (nFrames == 0)
        nFrames = 1;
    end

    % 1-based indexing for accessing evm
    % Limit frame-averaging to complete frames only
    slotRange = slotRange+1;
    slotRange = slotRange(slotRange <= (nFrames*10*ofdmInfo.SlotsPerSubframe));
    % Loop through each frame, selecting the frames with higher RMS (when
    % measuring 3GPP EVM)
    for frameIdx = 0:nFrames-1
        frameLowEVM = hRawEVM(cat(1,evm(1,slotRange).EV));
        frameEVM(frameIdx+1) = frameLowEVM;
        if evm3GPP
            frameHighEVM = hRawEVM(cat(1,evm(2,slotRange).EV));
            if frameHighEVM.RMS > frameLowEVM.RMS
                frameEVM(frameIdx+1) = frameHighEVM;
            end
        end
    end

    overallEVM = hRawEVM(cat(1,frameEVM(:).EV));
    evmInfo.OverallEVM = overallEVM;
end