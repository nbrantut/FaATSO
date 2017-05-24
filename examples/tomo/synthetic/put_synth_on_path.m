%PUT_SYNTH_ON_PATH
%
%This routine allows to access the synthetic programs (to generate 
%synthetic data for use in the tomography examples) from
%within Matlab by adding the bin/ folder containing the executable onto the
%PATH environment variable used by Matlab.

here = pwd;

%add to environment:
pth = getenv('PATH');
if isempty(strfind(here,pth))
    setenv('PATH', [pth ':' here '/bin/']);
end