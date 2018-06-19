%load_AxIS_file [[DEPRECATED]] reads raw recording files created by AxIS
%
%  Legal forms:
%     data = load_AxIS_file(filename);
%     data = load_AxIS_file(filename, well);
%     data = load_AxIS_file(filename, electrode);
%     data = load_AxIS_file(filename, well, electrode);
%     data = load_AxIS_file(filename, timespan);
%     data = load_AxIS_file(filename, well, timespan);
%     data = load_AxIS_file(filename, electrode, timespan);
%     data = load_AxIS_file(filename, well, electrode, timespan);
%
%  Required arguments:
%    filename    Pathname of the file to load
%
%  Optional arguments:
%    well        String listing which wells (in a multiwell file) to load.
%                Format is a comma-delimited string with whitespace ignored, e.g.
%                'A1, B2,C3' limits the data loaded to wells A1, B2, and C3.
%                Also acceptable: 'all' to load all wells.
%                If this parameter is omitted, all wells are loaded.
%                For a single-well file, this parameter is ignored.
%
%    electrode   Which electrodes to load.  Format is either a comma-delimited string
%                with whitespace ignored (e.g. '11, 22,33') or a single channel number;
%                that is, a number, not part of a string.
%                Also acceptable: 'all' to load all channels and 'none', '-1', or -1
%                to load no data (returns only header information).
%                If this parameter is omitted, all channels are loaded.
%
%    timespan    Span of time, in seconds, over which to load data.  Format is a two-element
%                array, [t0 t1], where t0 is the start time and t1 is the end time and both
%                are in seconds after the first sample in the file.  Samples returned are ones
%                that were taken at time >= t0 and <= t1.  The beginning of the file
%                is at 0 seconds.
%                If this parameter is omitted, the data is not filtered based on time.
%

function aData = load_AxIS_file(aFileName, varargin)
    aData = AxisFile(aFileName).DataSets(1).load_as_legacy_struct(varargin);
end

