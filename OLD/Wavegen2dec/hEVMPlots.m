function plots = hEVMPlots(varargin)
% hEVMPlots Function for EVM related plots
%   PLOTS = hEVMPlots(EVMINFO,EQSYM,REFSYM,BWPIDX)
%   returns an array of objects PLOTS, used for plotting the EVM statistics.
%   EVMINFO is a struct containing EVM statistics with fields.
%   SubcarrierRMS      - Root Mean square EVM (RMS) per subcarrier
%                        (Column vector of N subcarriers)
%   SubcarrierPeak     - Per subcarrier peak EVM
%                        (Column vector of N subcarriers)
%   SymbolRMS          - Root Mean square EVM per symbol
%                        (Column vector of S symbols)
%   SymbolPeak         - Per symbol peak EVM
%                        (Column vector of S symbols)
%   SlotRMS            - Root Mean square EVM per slot
%                        (Column vector of x slots)
%   SlotPeak           - Per slot peak EVM
%                        (Column vector of x slots)
%   OverallEVM         - structure containing EVM statistics for the
%                        overall waveform. It contains these fields:
%       EV             - Error vector for the overall waveform
%                        (Array of 1-by-n layers)
%       RMS            - RMS EVM for the overall waveform
%                        (Scalar)
%       Peak           - Peak EVM for the overall waveform
%                        (Scalar)
%   EVMGRID            - Error vector for the selected window edge
%                        EVM is averaged across layers where EVMGRID is a
%                        (2-D array of dimensions N subcarriers-by-Symbols)
%   EVM                - Cell array of EVM statistics. The first
%                        cell contains EVM per PDSCH slot and per window
%                        location, while the second cell contains peak and
%                        RMS EVM for the overall waveform). When EVM is not
%                        measured according to 3GPP standard, only one
%                        window location is used
%                        (2-by-Number of Slots)
%
%   EQSYM              - Cell array of IQ constellations for low & high EVM
%                        window locations
%                        (2-by-Number of Slots)
%   REFSYM             - Cell array of reference IQ constellations for low
%                        and high EVM window locations
%                        (2-by-Number of Slots)
%   LABEL              - Character string representing the physical channel
%                        with these values: 'PDSCH', 'PDCCH' or 'PUSCH'
%
%   PLOTS = hEVMPlots(...,EMISSIONS)
%   EMISSIONS is an optional structure containing in-band emission
%   statistics as defined in TS 38.101-1 (FR1) / TS 38.101-2 (FR2), Annex
%   F.3. It contains these fields.
%   DeltaRB            - Array of unallocated PUSCH RBs
%                        Each RB is expressed as an offset relative to the
%                        nearest edge of the PUSCH allocation
%                        Row vector of delta RBs
%   Absolute           - Absolute in-band emission in each delta RB vs slot.
%                        Array of deltaRBs vs number of slots
%   Relative           - Relative in-band emission in each delta RB vs slot.
%                        Array of deltaRBs vs number of slots

% Copyright 2020-2023 The MathWorks, Inc.

    evmInfo = varargin{1};
    eqSym= varargin{2};
    refSym = varargin{3};
    label = varargin{4};
    if nargin == 5
        emissions= varargin{5};
    end

    % Obtain plot positions
    [repositionPlots,plotPositions] = hNRPlotPositions();

    % Copy EVM statistics to be used later for plot purposes
    evmGrid = evmInfo.EVMGrid;
    evmSymbolRMS = evmInfo.SymbolRMS;
    evmSymbolPeak = evmInfo.SymbolPeak;
    evmSlotRMS = evmInfo.SlotRMS;
    evmSlotPeak = evmInfo.SlotPeak;
    evmSubcarrierRMS = evmInfo.SubcarrierRMS;
    evmSubcarrierPeak = evmInfo.SubcarrierPeak;
    if ~ischar(label)
        bwpIdx = label;
        label = '';
    else
        bwpIdx = evmInfo.BandwidthPartID;
    end

    % Plot EVM versus OFDM symbol
    fh = figure;
    if repositionPlots
        fh.Position = plotPositions(2,:);
    end
    subplot(3,1,1)
    plot(0:length(evmSymbolRMS)-1,evmSymbolRMS,'.-',0:length(evmSymbolPeak)-1,evmSymbolPeak,'.:');
    title([label,' EVM vs OFDM symbol, BWP index : ',num2str(bwpIdx)]);
    grid on; xlabel('Symbol number'); ylabel('EVM (%)'); legend('rms EVM','peak EVM','Location','bestoutside')
    yLim = max(evmSymbolPeak)*1.1;
    yTickStep = linspace(0,yLim,3);
    yticks(yTickStep)
    axis([0 length(evmSymbolRMS)-1 0 yLim])

    %Plot EVM versus slot
    subplot(3,1,2)
    slotIdx = 0:length(evmSlotRMS)-1;
    plot(slotIdx,evmSlotRMS,'.-',slotIdx,evmSlotPeak,'.:');
    title([label,' EVM vs Slot, BWP index : ',num2str(bwpIdx)]);
    grid on; xlabel('Slot number'); ylabel('EVM (%)'); legend('rms EVM','peak EVM','Location','bestoutside')
    yLim = max(evmSlotPeak)*1.1;
    yTickStep = linspace(0,yLim,3);
    yticks(yTickStep)
    axis([0 slotIdx(end)+1 0 yLim])

    % Plot EVM versus subcarrier
    subplot(3,1,3)
    scIdx = 0:length(evmSubcarrierRMS)-1;
    plot(scIdx,evmSubcarrierRMS,'.-',scIdx,evmSubcarrierPeak,'.:')
    title([label,' EVM vs Subcarrier, BWP index : ',num2str(bwpIdx)]);
    grid on; xlabel('Subcarrier number'); ylabel('EVM (%)'); legend('rms EVM','peak EVM','Location','bestoutside')
    yLim = max(evmSubcarrierPeak)*1.1;
    yTickStep = linspace(0,yLim,3);
    yticks(yTickStep)
    axis([0 length(evmSubcarrierPeak)-1 0 yLim])

    % Plot EVM resource grid
    % Display a 3-D surface plot using evmGrid
    % EVM resource grid
    evmGridFigure = figure('Visible','off');
    evmGridFigure.Name = 'EVM (%)';
    evmGridFigure.NumberTitle = 'off';
    figure(evmGridFigure);
    surf(evmGrid);
    shading flat;
    if repositionPlots
        evmGridFigure.Position = plotPositions(1,:);
    end
    evmGridFigure.CurrentAxes.Title.String = {[label,' EVM Resource Grid, BWP index : ',num2str(bwpIdx)]}; 
    evmGridFigure.CurrentAxes.XGrid = 'on';
    evmGridFigure.CurrentAxes.YGrid = 'on';
    evmGridFigure.CurrentAxes.ZGrid = 'on';
    evmGridFigure.CurrentAxes.XLabel.String = 'OFDM symbols';
    evmGridFigure.CurrentAxes.XLabel.FontSize = 8;
    evmGridFigure.CurrentAxes.YLabel.String = 'Subcarriers';
    evmGridFigure.CurrentAxes.YLabel.FontSize = 8;
    evmGridFigure.CurrentAxes.ZLabel.String = 'EVM (%)';
    evmGridFigure.CurrentAxes.ZLabel.FontSize = 8;
    evmGridFigure.CurrentAxes.View = [-30 60];
    evmGridFigure.CurrentAxes.XLim = [0 size(evmGrid,2)+1];
    evmGridFigure.CurrentAxes.YLim = [0 size(evmGrid,1)+1];
    evmGridFigure.CurrentAxes.ZLim = [0 max(evmGrid(:))*1.1];

    % Concatenate all plots related structures
    plots = {fh; evmGridFigure};

    % For in-band emissions, prepare the location of the plot 
    if nargin == 5 && (~isempty(emissions.DeltaRB) && ~isempty(emissions.Absolute))

        % Setup figure handle
        emissionsFigure = figure;

        % Plot in-band emissions
        plot(emissions.Absolute.','*-');

        % Set plot styles
        emissionsFigure.NumberTitle = 'off';
        emissionsFigure.Name = 'Absolute in-band emissions';
        if repositionPlots
            emissionsFigure.Position = plotPositions(5,:);
        end
        emissionsFigure.CurrentAxes.Title.String = {['Absolute In-band Emissions For Each Unallocated RB, BWP : ',num2str(bwpIdx)]};
        emissionsFigure.CurrentAxes.YLim(1) = -0.5e-07;
        emissionsFigure.CurrentAxes.XGrid = 'on';
        emissionsFigure.CurrentAxes.YGrid = 'on';
        emissionsFigure.CurrentAxes.XLabel.String = 'Slot';
        emissionsFigure.CurrentAxes.YLabel.String = 'Absolute in-band emissions';

        % Add legend
        nRB = length(emissions.DeltaRB);
        s = cell(1,nRB);
        for k = 1:nRB
            s{k} = sprintf('\\Delta_R_B = %d',emissions.DeltaRB(k));
        end
        legend(s,'Location','Best','TextColor',[175 175 175]/255,'EdgeColor',[175 175 175]/255);

        plots = {plots;emissionsFigure};
    end

    % Plot the constellation
    constPlotFh = figure;
    if repositionPlots
        constPlotFh.Position = plotPositions(3,:);
    end
    eqPlot = plot(eqSym,'.'); hold on
    refPlot = plot(refSym,'+');
    title([label,' Equalized Symbols Constellation, BWP index : ',num2str(bwpIdx)]);
    plots = {plots;eqPlot;refPlot};
end

%hNRPlotPositions returns plot sizes and positions based on screen resolution
%   [REPOSITIONPLOTS,POSITION] = hNRPlotPositions() determines the
%   support for plot tiling based on the user's display resolution. If plot
%   tiling is supported, REPOSITIONPLOTS is returned as true and
%   POSITION contains an array of vectors of the form [left bottom width
%   height] indicating the X-Y position and size for each plot. The plots
%   are arranged in a 2x3 tiled arrangement. If plot tiling is not
%   supported, REPOSITIONPLOTS is returned as false and POSITION is not
%   valid.
%
%   The plots are arranged in the following manner:
%    2x3 Tile
%            Plot 5    Plot 1    Plot 2
%            Plot 6    Plot 3    Plot 4

function [repositionPlots,position] = hNRPlotPositions()

    su=get(0,'Units');
    set(0,'Units','pixels');
    res = get(0,'ScreenSize');
    set(0,'Units',su);

    % Minimum PC - windows resolution,
    % arbitrarily chosen for plot positioning
    minres = 1280; 
    if (res(3)>minres)
        xpos = fix(res(3)*[1/2; 3/4; 1/2; 3/4; 1/4; 1/4]);
        ypos = fix(res(4)*[1/2; 1/2; 1/16; 1/16; 1/2; 1/16]);
        xsize = (xpos(2) - xpos(1) - 20)*[1; 1; 1; 1; 1; 1];   
        ysize = fix(xsize(1) * 5 / 6)*[1; 1; 1; 1; 1; 1];
        position = [xpos ypos xsize ysize];
        repositionPlots = true;
    else
        position = zeros(6,4);
        repositionPlots = false;
    end
end