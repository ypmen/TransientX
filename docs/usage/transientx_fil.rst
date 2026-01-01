==================
transientx_fil
==================
The options are described below:

- Read file options:
  	- ``--jump`` or ``-j``: the skip time at the beginning and end of the data file, in seconds. For example, `-j 10 10` skips the first and last 10 seconds of data.
	- ``--seglen`` or ``-l``: the length of each data chunck, in seconds. 
	- ``--wts``: apply weights in the psrfits file.
	- ``--scloffs``: apply scale and offset in the psrfits file.
	- ``--zero_off``: remove the offset in the psrfits file.
	- ``--cont``: skip checking data continuity between multiple input files.
	- ``--psrfits``: flag to indicate that the input file is in psrfits format. If not provided, the input file is assumed to be in sigproc filterbank format.
	- ``--input`` or ``-f``: path to the input data files. Multiple files can be provided, separated by space.

- Pipeline options:
	- ``--config`` or ``-c``: path to a json configuration file specifying the pipeline parameters. If this option is provided, other pipeline step options will be ignored.

- Downsample options (:doc:`preprocess/downsample <../pipeline/preprocess>`):
  	- ``--td``: time downsampling factor.
  	- ``--fd``: frequency downsampling factor.

- Baseline removal options (:doc:`preprocess/baseline_removal <../pipeline/preprocess>`):
	- ``--baseline w1 w2``: window sizes for the moving average filters used in baseline removal, in number of seconds. The first window size ``w1`` is not used anymore.

- RFI mitigation options (:doc:`rfi_mitigation <../pipeline/rfi_mitigation>`):
  	- ``--zapthre``: threshold for the skewness–kurtosis filter, expressed in interquartile-range threshold.
	- ``--rfi`` or ``-z``: RFI mitigation method. Options include 'zap flow fhigh' (zap channels manually), 'zero' (zero-DM filter), 'zdot' (zero-DM matched filter), 'kadaneF td fd' (kadaneF filter), 'mask td fd' (mask filter). The options can be combined, e.g., `-z zdot kadaneF 4 8 zap 1000 1100 zap 1100 1200`.
	- ``--widthlimit``: window width limit for the kadaneF filter, in number of seconds, below which the data will not be zapped.
	- ``--threMask``: threshold for the mask filter, expressed in interquartile-range threshold.
	- ``--threKadaneF``: threshold for the kadaneF filter, expressed in S/N.
	- ``--fill``: flag to indicate whether to fill the zapped data with random noise or mean.

- Dedispersion options (:doc:`dedispersion <../pipeline/dedispersion>`):
	- ``--dms``: starting DM for the DM-search range, in pc cm^-3.
	- ``--ddm``: DM step for the DM-search range, in pc cm^-3.
	- ``--ndm``: number of DMs to search.
	- ``--overlap``: ratio of data overlap between two adjacent data blocks during dedispersion, expressed as a fraction of the data block length.
	- ``--ddplan``: path to a file specifying a custom DM plan. If this option is provided, the options ``--dms``, ``--ddm``, and ``--ndm`` will be ignored.
	- ``--mean``: mean value to dump the dedispersed time series to disk. 
	- ``--std``: standard deviation value to dump the dedispersed time series to disk.
	- ``--nbits``: number of bits to use when dumping the dedispersed time series to disk.
	- ``--savetim``: flag to indicate whether to save the dedispersed time series to disk.
	- ``--format``: format to use when saving the dedispersed time series to disk. Options include 'sigproc', 'presto' and 'pulsarx'.

- Matched filtering options (:doc:`matched_filtering <../pipeline/matched_filtering>`):
	- ``--minw``: minimum boxcar width for matched filtering, in number of seconds.
	- ``--maxw``: maximum boxcar width for matched filtering, in number of seconds.
	- ``--snrloss``: maximum tolerated S/N loss due to boxcar width step.
	- ``--iqr``: when set, uses interquartile range instead of standard deviation to estimate noise level for S/N calculation. This method is slower but more robust when strong signals are present.

- Clustering options (:doc:`clustering <../pipeline/clustering>`):
	- ``--thre``: S/N threshold for candidate selection after matched filtering.
	- ``--radius`` or ``-r``: clustering radius in the DM-time space, in units of milliseconds.
	- ``--neighbors`` or ``-k``: minimum number of neighbors within the clustering radius to form a cluster.
	- ``--maxncand``: maximum number of candidates to keep, ranked by S/N, in each data segment.
	- ``--minpts``: minimum number of points to form a cluster.
	- ``--drop``: flag to indicate whether to drop candidates with boxcar width equal to the maximum search width.

- Candidate plotting options (:doc:`candidate_plotting <../pipeline/candidate_plotting>`):
	- ``--ra``: right ascension of the observation target, in hh:mm:ss format. If not provided, the RA value in the data header will be used.
	- ``--dec``: declination of the observation target, in dd:mm:ss format. If not provided, the DEC value in the data header will be used.
	- ``--ibeam`` or ``-i``: beam of the observation. If not provided, the period value in the data header will be used.
	- ``--telescope``: name of the telescope. If not provided, the telescope value in the data header will be used.
	- ``--incoherent``: flag indicating that the beam is an incoherent beam (ifbf) or coherent beam (cfbf) for the radio array.
	- ``--source_name``: name of the observation target. If not provided, the source name in the data header will be used.
	- ``--rootname``: root name for the output files.
	- ``--saveimage``: flag to indicate whether to save the candidate plots in fits format.
