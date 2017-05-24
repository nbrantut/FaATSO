%PUT_FAATSO_ON_PATH
%
%This routine allows to access the faatso programs (tomo and fmm) from
%within Matlab by adding the bin/ folder of the faatso project onto the
%PATH environment variable used by Matlab.
%
%At the moment you HAVE TO MODIFY this routine manually to set Faatso's
%root directory.

%set the root directory of Faatso here (modify the path according to where
%Faatso is installed):
faatso_root_dir = '/Volumes/Data/Projects/Tomography/Faatso/';


%add to environment:
pth = getenv('PATH');
if isempty(strfind(pth, faatso_root_dir))
    setenv('PATH', [pth ':' faatso_root_dir 'bin/']);
end