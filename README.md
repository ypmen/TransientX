# TransientX

TransientX is a one command line high performance transient search software (filtterbank and psrfits (--psrfits) are both supported).

## Dependencies

- boost > 1.56
- PlotX (https://github.com/ypmen/PlotX)
- sofa (https://www.iausofa.org/2020_0721_C/sofa_c-20200721.tar.gz)

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
### replot_fil: filter out the bad candidates
```
replot_fil -v -t 4 --zapthre 3.0 --td 1 --fd 1 --dmcutoff 3 --widthcutoff 0.1 --snrcutoff 7 --snrloss 0.1 --zap --zdot --kadane 8 4 7 --candfile test.cands --clean -f *.fil
```
### example ddplan
[ddplan.txt](examples/ddplan.txt)
### example candidate

![exmaple](examples/example.png)

## Acknowledgement
Thanks for very helpful suggestions and bug reports from Guðjón Henning Hilmarsson, Emma Carli and PKU pulsar group.

