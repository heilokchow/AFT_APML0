# AFT_APML0 Model

## Download repository
```
git clone "https://github.com/heilokchow/AFT_APML0/"
```

## Build illustration used program in unix system (n = 100, p = 200)
```
make test
```

## Run illustration used program in background (20 seconds)
```
./test > log.txt &
```

## Output file (AL_path_lasso.txt, AL_path_apml0.txt)
The output file is comma delimited which can be inported to excel. The column 2-8 are insample CV-error, tunning parameter,
number of paramters selected, number of True paramters selected, outsample c-index, insample c-index, outsample CV-error.
Later columns are estimated parameters with corresponding model.
