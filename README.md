# PMS2_vaR

PMS2_vaR is a framework to assist the inclusion of PMS2 mutational analysis in routine NGS diagnostic pipelines.

### Prerequisites ###

The following tools should be properly installed.
- samtools (v1.10)
- picard (v.2.26.4.jar)
- bwa (0.7.17)
- vardictJava: https://github.com/AstraZeneca-NGS/VarDictJava

Also, R/Bioconductor should be installed with at least these packages: Biostrings, dplyr, httr, jsonlite, optparse, purrr, Rsamtools, stringr, universalmotif, VariantAnnotation, vcfR

### How to use ###

1. **Get Code**

```
git clone https://github.com/emunte/PMS2_vaR.git
```

2. **Configure tools.yaml**
3. **Obtain a modified reference genome file without** ***PMS2CL*** **sequence** . This step only needs to be done once.
```
cd PMS2_vaR
Rscript modify_reference.R [-r --reference reference_genome_path] [-p --pms2CLfasta  full_path_to_PMS2CL_fasta_file] [-o --outputdir path_to_output_directory]
```

4. **Configure vardicjavaParams.yaml**


5. **Launch PMS2_vaR**
```
cd PMS2_vaR
Rscript run_PMS2_vaR.R [-t tools_file] [-b bam.txt_path] [-r modified_reference_genome_path][-g genome_assembly] [-v vardictjava_path] [-n --samplesname samplesname] [-o --outputdir path_to_output_directory]
```

### Output ###
An excel file is obtained in the output directory (-o). 
An algorithm was designed to recommend if a (likely) pathogenic variant would need confirmation by LR-PCR and, if located in PMS2, should be reported
![See Algorithm](https://github.com/emunte/PMS2_vaR/blob/main/WF/Figure2_08012024.png)




