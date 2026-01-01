# TransientX

TransientX is a one command line high performance transient search software (filtterbank and psrfits (--psrfits) are both supported). A detailed user guide can be found in https://transientx.readthedocs.io/en/latest/.

## Dependencies

- boost > 1.56
- PlotX (https://github.com/ypmen/PlotX)
- XLibs (https://github.com/ypmen/XLibs)
- sofa (https://www.iausofa.org/2020_0721_C/sofa_c-20200721.tar.gz) or erfa (https://github.com/liberfa/erfa)

## Installation
- ./bootstrap
- ./configure --prefix=[install_path] LDFLAGS="-L/path_to_sofa" CPPFLAGS="-I/path_to_sofa"
- make
- make install

## Docker
```
docker pull ypmen/pulsarx
```

## Example usage (**please let me know if you meet any problems!**)
- export YMW16_DIR=/install_path/src/ymw16
- transientx_fil -h for help

### transientx_fil: search for single pulse
```
transientx_fil -v -o test -t 4 --zapthre 3.0 --fd 1 --overlap 0.1 --ddplan ddplan.txt --thre 7 --maxw 0.1 --snrloss 0.1 -l 2.0 --drop -z kadaneF 8 4 zdot -f *.fil
```
or using a json file
```
transientx_fil -v -o test -t 4 -l 2.0 -c examples/transientx_config.json -f *.fil

```
- --saveimage (save image for AI classifier)
### replot_fil: filter out the bad candidates
```
replot_fil -v -t 4 --zapthre 3.0 --td 1 --fd 1 --dmcutoff 3 --widthcutoff 0.1 --snrcutoff 7 --snrloss 0.1 --zap --zdot --kadane 8 4 7 --candfile test.cands --clean -f *.fil
```

### example ddplan
[ddplan.txt](examples/ddplan.txt)
### example candidate

![exmaple](examples/example.png)

### Skills
- When searching in high-resolution data (e.g., <50us), please apply multiple search strategies to detect bursts with various timescales: (1) For microsecond bursts, use a smaller data block size (-l), a smaller DBSCAN radius (-r), a smaller maximum search pulse width (--maxw), etc;(2) For millisecond bursts, you can use larger values for the above parameters. This approach can help avoid slowing down the clustering process, as discussed in the Benchmark Section of the TransientX paper [TransientX](https://arxiv.org/abs/2401.13834).
  
- use the option --iqr when there are many strong bursts, which can improve the estimation of baseline and noise level.

## Citation

Please cite [https://arxiv.org/abs/2401.13834](https://arxiv.org/abs/2401.13834), if you are using TransientX.

## Acknowledgement
Thanks for very helpful suggestions and bug reports from Guðjón Henning Hilmarsson, Emma Carli and PKU pulsar group.

