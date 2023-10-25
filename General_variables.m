%% General variables
%% Simulation variables

SNRdB = 10;                % SNR in dB

totalNoSlots = 3;          % Number slots containing information simulated

perfectEstimation = false; % Perfect synchronization and channel estimation
rng("default");            % Set default random number generator for repeatability

%Files for read and write operations
fileID = fopen('Bits.bin','r'); %send bit file
fileID2 = fopen('RecBits.bin','w'); %recieve bit file

%% Carrier config

carrier = nrCarrierConfig
carrier.NCellID = 1;
carrier.SubcarrierSpacing = 15;
carrier.CyclicPrefix = 'normal';
carrier.NSizeGrid = 
carrier.NSlot =
carrier.NFrame =
carrier.IntraCellGuardBands =

%Read-Only variables
%SymbolsPerSlot: Number of OFDM symbols per slot bases on CyclePrefix
%SlotsPerSubFrame: Number of Slots per 1ms frame, based of the subcarrierspacing
%SlotsPerFrame: Number of Slots per 10ms frame, based of the subcarrierspacing


pdsch = nrPDSCHConfig;
pdsch.Modulation = "16QAM";
pdsch.NumLayers = 1;
pdsch.PRBSet = 0:carrier.NSizeGrid-1;     % Full band allocation

pdsch
