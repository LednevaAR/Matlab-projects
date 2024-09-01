%% FM Broadcast Receiver in Simulink
% This example shows how to build an FM mono or stereo receiver using
% Simulink(R) and Communications Toolbox(TM). You can either use captured
% signals, or receive signals in real time using the RTL-SDR, ADALM-PLUTO
% or USRP radio. 

% Copyright 2013-2021 The MathWorks, Inc.

%% Required Hardware and Software
% To run this example using captured signals, you need the following
% software:
%
% * <https://www.mathworks.com/products/simulink.html _Simulink_>
% * <https://www.mathworks.com/products/communications.html _Communications
% Toolbox(TM)_>
%
% To receive signals in real time, you also need one of the following
% hardware:
%
% * RTL-SDR radio and the corresponding software
% <https://www.mathworks.com/hardware-support/rtl-sdr.html _Communications
% Toolbox Support Package for RTL-SDR Radio_>
%
% * ADALM-PLUTO radio and the corresponding software
% <https://www.mathworks.com/hardware-support/adalm-pluto-radio.html
% _Communications Toolbox Support Package for ADALM-PLUTO Radio_>
% 
% * USRP radio and the corresponding software
% <https://in.mathworks.com/hardware-support/usrp.html
% _Communications Toolbox Support Package for USRP Radio_>
% 
% For a full list of Communications Toolbox supported SDR platforms, refer
% to Supported Hardware section of
% <https://www.mathworks.com/discovery/sdr.html _Software Defined Radio
% (SDR)_>.

%% Introduction
% For an introduction to the FM broadcasting technology and demodulation of
% these signals, refer to the <docid:comm_ug#example-FMReceiverExample FM
% Broadcast Receiver> example.

%% Running the Example
% To run the example using captured signals, select the |FM Broadcast
% Captured Signal| block as the source using the |Signal Source Selector|
% block. Then click the run button.
% 
% To run the example using the RTL-SDR radio, ADALM-PLUTO radio or USRP
% radio as the source, select the corresponding |RTL-SDR Receiver|,
% |ADALM-PLUTO Radio Receiver|  or |USRP Radio Receiver| block as the
% source using the |Signal Source Selector| block. Double-click the |Center
% Frequency (MHz)| block and select the value to the center frequency to a
% broadcast FM radio station near you.
%
% If you hear some dropouts or delay in the sound, run the model in
% Accelerator mode.  From the model menu, select Simulation->Accelerator,
% then click the run button.  If you still experience dropouts or delay in
% Accelerator mode, try running the model in Rapid Accelerator mode.

%% Receiver Structure
% The receiver structure in this block diagram represents the
% FMReceiverSimulinkExample.slx model used in this example. The processing
% has three main parts: signal source, FM broadcast demodulation, and audio
% output. 
%
modelName = 'FMReceiverSimulinkExample';
open_system(modelName);
set_param(modelName,'SimulationCommand','update');

%%
% *Signal Source*
%
% This example can use three signal sources:
%
% # ''Captured Signal'': Over-the-air signals written to a file and sourced
% using a <docid:comm_ref#bvb0lb6 Baseband File Reader> block at 228e3
% samples per second.
% # ''RTL-SDR Radio'': RTL-SDR radio running at 200e3 samples per second.
% Set the center frequency to a broadcast FM radio station near you.
% # ''ADALM-PLUTO Radio Receiver'': ADALM-PLUTO radio running at 200e3
% samples per second. Set the center frequency to a broadcast FM radio
% station near you.
% # ''USRP Radio Receiver'': USRP radio running at 200e3 samples per
% second. Set the center frequency to a broadcast FM radio station near
% you. 

%%
% *FM Broadcast Demodulation*
%
% The baseband samples received from the signal source are processed by the
% |FM Broadcast Demodulation Baseband| block. This block converts the input
% sampling rate of 228 kHz to 45.6 kHz, the sampling rate for your host
% computer's audio device. According to the FM broadcast standard in the
% United States, the de-emphasis lowpass filter time constant is set to 75
% microseconds. This example processes received mono signals. The
% demodulator can also process stereo signals.
%
% To perform stereo decoding, the FM Broadcast Demodulator Baseband object
% uses a peaking filter which picks out the 19 kHz pilot tone from which
% the 38 kHz carrier is created. Using the resulting carrier signal, the FM
% Broadcast Demodulator Baseband block downconverts the L-R signal,
% centered at 38 kHz, to baseband. Afterwards, the L-R and L+R signals pass
% through a 75 microsecond de-emphasis filter. The FM Broadcast Demodulator
% Baseband block separates the L and R signals and converts them to the
% 45.6 kHz audio signal.

%%
% *Audio Device Writer*
%
% Play the demodulated audio signals through your computer's speakers using
% the |Audio Device Writer| block.

%% Further Exploration
% To further explore the example, you can vary the center frequency of the
% RTL-SDR radio, ADALM-PLUTO radio or USRP radio and listen to other radio
% stations using the |Center Frequency (MHz)| block.
%
% You can set the Stereo property of the |FM Broadcast Demodulator
% Baseband| block to true to process the signals in stereo fashion and
% compare the sound quality.

%% Selected Bibliography 
% 
% https://en.wikipedia.org/wiki/FM_broadcasting

close_system(modelName,0)
