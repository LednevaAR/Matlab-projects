function SimParams = plutoradioqpsktransmitter_init
%   Copyright 2017 The MathWorks, Inc.

%% General simulation parameters
SimParams.Rsym = 0.2e6;             % Symbol rate in Hertz
SimParams.ModulationOrder = 4;      % QPSK alphabet size
SimParams.Interpolation = 2;        % Interpolation factor
SimParams.Decimation = 1;           % Decimation factor
SimParams.Tsym = 1/SimParams.Rsym;  % Symbol time in sec
SimParams.Fs   = SimParams.Rsym * SimParams.Interpolation; % Sample rate

%% Frame Specifications
% [BarkerCode*2 | 'Hello world 000\n' | 'Hello world 001\n' ...];
SimParams.BarkerCode      = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];     % Bipolar Barker Code
SimParams.BarkerLength    = length(SimParams.BarkerCode);
SimParams.HeaderLength    = SimParams.BarkerLength * 2;                   % Duplicate 2 Barker codes to be as a header
SimParams.Message         = 'Hello world';
%SimParams.MessageLength   = length(SimParams.Message) + 5;                % 'Hello world 000\n'...
%SimParams.NumberOfMessage = 100;                                       % Number of messages in a frame
SimParams.MessageLength = 12000;
SimParams.PayloadLength = SimParams.MessageLength;
SimParams.FrameSize = (SimParams.HeaderLength + SimParams.PayloadLength) ...
    / log2(SimParams.ModulationOrder);
%SimParams.PayloadLength   = SimParams.NumberOfMessage * SimParams.MessageLength * 7; % 7 bits per characters
%SimParams.FrameSize       = (SimParams.HeaderLength + SimParams.PayloadLength) ...
%    / log2(SimParams.ModulationOrder);                                    % Frame size in symbols
SimParams.FrameTime       = SimParams.Tsym*SimParams.FrameSize;

%% Tx parameters
SimParams.RolloffFactor     = 0.5;                                        % Rolloff Factor of Raised Cosine Filter
SimParams.ScramblerBase     = 2;
SimParams.ScramblerPolynomial           = [1 1 1 0 1];
SimParams.ScramblerInitialConditions    = [0 0 0 0];
SimParams.RaisedCosineFilterSpan = 10; % Filter span of Raised Cosine Tx Rx filters (in symbols)

%% Message generation
% Pluto transmitter parameters
SimParams.PlutoCenterFrequency      = 915e6;
SimParams.PlutoGain                 = 0;
SimParams.PlutoFrontEndSampleRate   = SimParams.Fs;
SimParams.PlutoFrameLength          = SimParams.Interpolation * SimParams.FrameSize;

% Simulation Parameters
SimParams.FrameTime = SimParams.PlutoFrameLength/SimParams.PlutoFrontEndSampleRate;
SimParams.StopTime  = 1000;
grass = double(imresize(imread('mount', 'jpg'), [50, 50]));
singleLineGarry = reshape(grass,1,[]);
singleLineGarryBIN = reshape(dec2bin(singleLineGarry, 8)',1,[])';
SimParams.array = [ones(1, 12000) bin2dec(singleLineGarryBIN)'];