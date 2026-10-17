# Running the analyses

The [run-code directory](https://github.com/byuzbasi/AdaReFun/tree/main/replication/code) contains the simulation and
real-data scripts, setup instructions, and checks. Numerical settings and
seeds are included with the code.

From the repository root:

```bash
cd replication/code
bash run.sh verify
```

Follow the [run guide](https://github.com/byuzbasi/AdaReFun/blob/main/replication/code/README.md) to install dependencies,
prepare a run, retrieve raw data, and execute the simulations or data analyses.
It also gives progress, resume, reporting, and output-validation commands.

The run scripts retrieve raw data from Zenodo and CRAN. The original preparation
and comparison scripts produce the analysis inputs; existing inputs can instead
be imported after checking their hashes in `inputs.csv`. No data or prediction
files are bundled. Data sources and licenses are listed in
[DATA_SOURCES.md](https://github.com/byuzbasi/AdaReFun/blob/main/replication/code/DATA_SOURCES.md).
