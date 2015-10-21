function f= midi2hz(midi)
%
% convert MIDI note number to frequency
%

f = power(2,(midi-69)/12)*440;
