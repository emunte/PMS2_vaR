# PMS2_vaR

PMS2_vaR is a framework to assist the inclusion of PMS2 mutational analysis in routine NGS diagnostic pipelines.

### Prerequisites ###

The following tools should be properly installed.
- samtools (v1.10)
- picard (v.2.26.4.jar)
- bwa (0.7.17)
- vardictJava: https://github.com/AstraZeneca-NGS/VarDictJava


### How to use ###

1. Get Code

```
git clone https://github.com/emunte/PMS2_vaR.git

```

2. **Configure tools.yaml**
3. **Configure vardicjavaParams.yaml**


4. Launch PMS2_vaR
```
cd PMS2_vaR
Rscript runPMS2_vaR.R [-t tools_file] [-b bam_txt_files] [-r reference_genome] [-v vardictjava]
```

### Output ###


