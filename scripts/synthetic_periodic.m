% out = synthetic_periodic (period, duration, duty)
%
% Generate a sampled periodic waveform, whose repeating pattern is given
% by a vector  duty,  or a square wave with duty-cycle  duty.
% The period of the underlying waveform is  period  sample times.
% The peak-peak amplitude of the output is  ampl,
% and if duty=0.5, the mean is ampl/2.
function out = synthetic_periodic (period, duration, duty, ampl, phase)
  if nargin < 5
    phase = 0;
    if nargin < 4
      ampl = 1;
      if nargin < 3
        duty = 0.5;
      end
    end
  end
  if numel (duty) == 1
    [num1, den1] = rat (   duty  * period, 1/duration);
    [num2, den2] = rat ((1-duty) * period, 1/duration);
    len = max (den1, den2);
    d = round (duty * len);
    duty = [ones(d,1); zeros(len - d, 1)];
  else
    duty = duty(:);               % force to column vector
    warning ('explicit duty not yet implemented');
    [n, d] = rat (length (duty) / period, 1/duration);
  end

  [scale_num, scale_den] = rat (length(duty)/period);
  len = duration * scale_num;
  replicas = ceil (len / (length (duty) * scale_den)) + 1;
  orig = repmat (interp (duty, scale_den), [replicas, 1]);
  offset = round (phase * length(duty) * scale_den);
  orig = orig (offset + 1: offset + len);
  out = sum (reshape (orig, [scale_num, duration]))';
  out = out * ampl / max (abs (out));
end

