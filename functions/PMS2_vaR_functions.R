
pms2_realignment <- function(args, tools, vardict, output.samples.dir){
bams <- read.table(args$bamTxt)
dir.create(output.samples.dir, showWarnings = FALSE)
for (b in seq_len(nrow(bams))){
  #Get sample name
  sample.name <-bams[b,1] %>% as.character %>%
    basename() %>%
    stringr::str_replace(".bam", "")
  full <- bams[b,1] %>% as.character()

  #Create results dir per sample
  results.dir.sample <- file.path(output.samples.dir, sample.name)
  dir.create (results.dir.sample, showWarnings = FALSE)

  #Create folder to store first approach
  results.dir.1 <- file.path(results.dir.sample, "Approach1")
  results.dir.2 <- file.path(results.dir.sample, "Approach2")
  dir.create(results.dir.1)
  dir.create(results.dir.2)

  #Location for different files
  ## APPROACH 1
  vcfvardict <- file.path (results.dir.1, paste0(sample.name, "freq", vardict$freq, "_VARDICT.vcf"))
  filtered.bam <- paste0(results.dir.1, "/filtered", sample.name, ".bam" )
  r1 <- paste0(results.dir.1,"/r1_", sample.name, ".fastq")
  r2 <- paste0(results.dir.1,"/r2_", sample.name, ".fastq")
  sam <- paste0(results.dir.1, "/segon_alineament", sample.name, ".sam")
  second.bam<- paste0(results.dir.1, "/second_alignment_", sample.name, ".bam")
  sorted.bam <- paste0(results.dir.1, "/second_alignment_", sample.name, "_sorted")


  ## APPROACH 2
  vcfvardict.E11 <- file.path (results.dir.2, paste0(sample.name, "freq", vardict$freq, "_VARDICT.vcf"))
  filtered.bam.E11 <- paste0(results.dir.2, "/filtrar", sample.name, ".bam" )
  filtered.bam.E11.1 <- paste0(results.dir.2, "/filtrar_resta_", sample.name, ".bam" )
  filtered.bam.E11.2 <- paste0(results.dir.2, "/filtrar_E11_", sample.name, ".bam" )
  filtered.bam.E11.3 <- paste0(results.dir.2, "/filtrar_E11_picard", sample.name, ".bam" )
  read.names <- paste0(results.dir.2, "/read_names", sample.name, ".txt" )
  r1.E11 <- paste0(results.dir.2,"/r1_", sample.name, ".fastq")
  r2.E11 <- paste0(results.dir.2,"/r2_", sample.name, ".fastq")
  sam.E11 <- paste0(results.dir.2, "/segon_alineament", sample.name, ".sam")
  second.bam.E11 <- paste0(results.dir.2, "/second_aligment_", sample.name, ".bam")
  third.bam.E11 <- paste0(results.dir.2, "/third_alignment", sample.name, ".bam")
  sorted.bam.E11 <- paste0(results.dir.2, "/third_alignment", sample.name, "sorted")


  #PIPELINE

  #APROACH 1 ----

  ### 1. Filter gen + pseudogen regions with samtools
  cmd1 <-paste(tools$samtools, "view -b -h", full, "7:6012830-6049024 7:6774736-6791432 >", filtered.bam)
  print(cmd1); system(cmd1)

  ### 2. Go back to fastq with picard
  cmd2 <- paste("java -jar", tools$picardJar, "SamToFastq -I", filtered.bam, "-F", r1 , "-F2", r2)
  print(cmd2); system(cmd2)

  ### 3. Realign with bwa mem
  cmd3 <- paste(tools$bwa,  "mem -t 1", fasta, r1, r2, ">", sam)
  print(cmd3); system(cmd3)

  ### 4. Convert sam to bam
  cmd4 <- paste(tools$samtools, "view -b -S", sam, ">", second.bam)
  print(cmd4); system(cmd4)

  ### 5. Sort the bam file amb index the bam file
  cmd5 <- paste (tools$samtools , "sort", second.bam, sorted.bam )
  print(cmd5); system(cmd5)
  cmd6 <- paste (tools$samtools , "index", paste0(sorted.bam, ".bam"))
  print(cmd6); system(cmd6)

  #delete intermediate files
  system(paste("rm", filtered.bam))
  system(paste("rm", r1))
  system(paste("rm", r2))
  system(paste("rm", sam))
  system(paste("rm", second.bam))

  # 6. Variant calling with VardictJava

  cmd7 <- paste (paste(paste0(tools$vardict,  "/bin/VarDict"),
                   "-G", fasta,
                   "-X", vardict$X,
                   "-q", vardict$phred_score,
                   "-m", vardict$missmatches,
                   "-f", (vardict$freq/100),
                   "-N", sample.name,
                   "-b ", paste0(sorted.bam, ".bam"),
                   "-R chr7:6012350-6049257:PMS2 | ",
                   paste0(tools$vardict , "/bin/teststrandbias.R | "),
                   paste0(tools$vardict, "/bin/var2vcf_valid.pl"),
                   "-N", sample.name,
                   "-E -A ",
                   "-d 8",
                   "-q", vardict$phred_score,
                   "-f", (vardict$freq/100), " > ",
                   vcfvardict ))
  print(cmd7); system(cmd7)

  # APPROACH 2 ----

  ## 1. Filter regions:

  ### 1.1 Outside E11 -> Filter with samtools the regions of interst
  cmd1 <-paste(tools$samtools, "view -b -h", full, "7:6012830-6022822 7:6029231-6049024 7:6774736-6775220 7:6781116-6791432 >", filtered.bam.E11.1)
  print(cmd1); system(cmd1)

  ### 1.2 Inside E11
  #### 1.2.1 Filter considering the invariable variants detailed in Gould et al
  cmd1.2 <- paste(tools$samtools, "view -b -h", full, "7:6026364-6026364 7:6026598-6026598 7:6026601-6026601 7:6026625-6026625 7:6026636-6026636 7:6026707-6026707 7:6027017-6027017 7:6027474-6027474>", filtered.bam.E11.2)
  print(cmd1.2); system(cmd1.2)
  #### 1.2.2 Get reads name that overlap with this 7 positions (without duplicates)
  cmd1.3 <- paste( tools$samtools, "view", filtered.bam.E11.2 , " | cut -f1 | sort | uniq >", read.names)
  print(cmd1.3); system(cmd1.3)

  #### 1.2.3 Filter reads from the original bam file taking into consideration the name (picard).
  cmd1.4 <- paste("java -jar", tools$picard, "FilterSamReads -I", full, "-O",  filtered.bam.E11.3, "-FILTER includeReadList -READ_LIST_FILE", read.names)
  print(cmd1.4); system(cmd1.4)


  ## 2.  Go back to fastq with picard:
  ### 2.1 Outside E11
  cmd2 <- paste("java -jar", tools$picardJar, "SamToFastq -I", filtered.bam.E11.1, "-F", r1.E11 , "-F2", r2.E11)
  print(cmd2); system(cmd2)

  ### 2.2 Inside E11 -> not necessary

  ## 3. Realign with bwa mem
  ### 3.1 Outside E11
  cmd3 <- paste(tools$bwa,  "mem -t 1", fasta, r1.E11, r2.E11, ">", sam.E11)
  print(cmd3); system(cmd3)

  ### 3.2 Inside E11 -> not necessary

  ## 4. Convert from sam to bam
  ### 4.1 Outside E11
  cmd4 <- paste(tools$samtools, "view -b -S", sam.E11, ">", second.bam.E11)
  print(cmd4); system(cmd4)

  ### 4.2 Inside E11 -> not necessary

  ## 5. Merge the two bams
  cmd4.1 <- paste("java -jar", tools$picard, "MergeSamFiles -I", second.bam.E11, "-I", filtered.bam.E11.3, "-O", third.bam.E11)
  print(cmd4.1); system(cmd4.1)

  ## 6. Sort + index
  cmd5 <- paste (tools$samtools , "sort", third.bam.E11, sorted.bam.E11 )
  print(cmd5); system(cmd5)
  cmd6 <- paste (samtools , "index", paste0(sorted.bam.E11, ".bam"))
  print(cmd6); system(cmd6)

  #remove intermediate files
  system(paste("rm", filtered.bam.E11.1))
  system(paste("rm", filtered.bam.E11.2))
  system(paste("rm", filtered.bam.E11.3))
  system(paste("rm", r1.E11))
  system(paste("rm", r2.E11))
  system(paste("rm", sam.E11))
  system(paste("rm", second.bam.E11))
  system(paste("rm", third.bam.E11))

  ## 7. Variant calling with VardictJava

  cmd7 <- paste (paste(paste0(tools$vardict,  "/bin/VarDict"),
                       "-G", fasta,
                       "-X", vardict$X,
                       "-q", vardict$phred_score,
                       "-m", vardict$missmatches,
                       "-f", (vardict$freq/100),
                       "-N", sample.name,
                       "-b ", paste0(sorted.bam.E11, ".bam"),
                       "-R chr7:6012350-6049257:PMS2 | ",
                       paste0(tools$vardict , "/bin/teststrandbias.R | "),
                       paste0(tools$vardict, "/bin/var2vcf_valid.pl"),
                       "-N", sample.name,
                       "-E -A ",
                       "-d 8",
                       "-q", vardict$phred_score,
                       "-f", (vardict$freq/100), " > ",
                       vcfvardict.E11 ))
  print(cmd7); system(cmd7)



  }
}


zip_tabix <- function(resultsDir, patro){
  directoris <- list.dirs(resultsDir, recursive = TRUE)
  lapply(directoris, function(x){
    message(paste("Tabix index created for", x))
    file <-list.files (x , pattern=patro , recursive= FALSE)
    message(length(file))
    if (length(file)>0){
      f <- file.path (x , file)
      cmd9 <- paste("bgzip", f)
      print(cmd9); system(cmd9)
      cmd10 <- paste ("tabix" , paste0(f, ".gz"))
      print(cmd10); system(cmd10)
    }

  })
}

convert_vcf_txt <- function(resultsDir, vcf.pattern, rng, args){
  files <- list.files (resultsDir, pattern= vcf.pattern, recursive= TRUE, full.names = TRUE)
  samples.name <- bams[,1] %>%
    as.character() %>%
    basename() %>%
    stringr::str_replace(".bam", "") %>%
    rep(each=2)

  for (f in seq_len(length(files))){ #mirem cada fitxer dins de files
    sample.name <- samples.name[f]
    tab <- TabixFile(files[f])
    #Read the complete vcf
    vcf <- read.vcfR(tab$path)
    # Read vcf only with inROI
    vcf.rng <- readVcf(tab, "hg19", param=rng)
    vcf.rng.exons <- readVcf(tab, "hg19", param=bedGR)
    exons.df <- data.frame(ID= names(vcf.rng.exons@rowRanges),
                           exon = vcf.rng.exons@rowRanges$paramRangeID %>% as.character())

    writeVcf(vcf.rng, file.path(dirname(files[f]), "rng_vardict.vcf"))
    vcf.rng <- read.vcfR(file.path(dirname(files[f]), "rng_vardict.vcf"))

    #tidy vcfs and make a dataframe
    #a) all vcf
    vcf.tidy <- vcfR2tidy(vcf, single_frame = TRUE)
    vcf.df <- as.data.frame(vcf.tidy$dat) %>%
      dplyr::mutate(ID= paste0(CHROM, ":", POS, "_", REF, "/", ALT))
    #b) in roi vcf
    vcf.tidy.rng <- vcfR2tidy(vcf.rng, single_frame = TRUE)
    vcf.df.rng <- as.data.frame(vcf.tidy.rng$dat)



    #Merge two dataframes to mark in ROIS
    sel <- vcf.df$POS %in% vcf.df.rng$POS
    all.variants <- merge(vcf.df, exons.df, by="ID", all= TRUE) %>%
      mutate(INROI= sel,
              ID_sample=all_of(sample.name)) %>%
      relocate(CHROM, POS, ID, REF, ALT,INROI, exon)

    write.table(all.variants, file.path(dirname(files[f]), "all.variants.txt") , sep="\t", row.names = F)
  }
}

merge_pipelines <- function(resultsDir){
  samples <- list.dirs(resultsDir, recursive=FALSE)
  for ( i in seq_len(length(samples))){
    message(paste("Annotating", basename(samples[i]), "variants"))
    all.vars <- list.files(samples[i], pattern= "all.variants.txt", recursive=TRUE, full.names = TRUE)
    file1 <- read.delim(all.vars[1], sep="\t", stringsAsFactors = FALSE)
    file2 <- read.delim(all.vars[2], sep="\t", stringsAsFactors = FALSE)
    both <- merge(file1, file2[,c(3, 8, 14)], by="ID", all=TRUE )
    all.variants.final <- cdnav2(dataset = both) %>%
      dplyr::relocate (cdna, prot, class_paralogous, paralogous_above60, class)
    all.variants.final.2 <- filter_variants(all.variants.final)
    final_results <- file.path(resultsDir, "final_results")
    dir.create(final_results)
    write.table(all.variants.final.2, file.path(final_results, paste0(samples[i], "PMS2_variants.txt")), row.names=FALSE, col.names=TRUE, quote=FALSE)
  }
}

cdnav2 <- function (dataset){
  final.file <- dataset %>% mutate(variant = ifelse(str_length(REF)==1,
                                                       paste0("NC_000007.13%3Ag.", POS, REF,">", ALT),
                                                       paste0("NC_000007.13%3Ag.",as.numeric(POS+1),"_",as.numeric(POS+str_length(REF)-1), "del", str_sub(REF,2, -1)))) %>%
    mutate (variant = ifelse(str_length(ALT)==1,
                             variant,
                             paste0("NC_000007.13%3Ag.",as.numeric(POS), "_",as.numeric(POS+1), "ins", str_sub(ALT,2, -1)))) %>%
    mutate(variant = ifelse(str_length(ALT)>1 & str_length(REF)>1,
                            paste0("NC_000007.13%3Ag.", as.numeric(POS),"_", as.numeric(POS+str_length(REF)-1),"delins", ALT ),
                            variant))
  final.file$cdna <- ":"
  final.file$prot <- ":"
  server_mutalyzerv3 <- "https://v3.mutalyzer.nl/api/normalize/"
  for (i in seq_len(nrow(final.file))){
    if (final.file$CHROM[i]==7){
      if(!(final.file$TYPE[i] %in% c("INV", "<DEL>", "<DUP>")) && stringr::str_length(final.file$variant[i])<100){
        mut <- NA
        while(is.na(mut)[1]| class(mut)=="try-error"){
          mut <- try(read_json(paste0(server_mutalyzerv3,final.file$variant[i])))
        }
        if(!is.null(mut$equivalent_descriptions)){
          number <- lapply(mut$equivalent_descriptions, function(x){
            stringr::str_detect(x, "NM_000535.7")
            }) %>%
            unlist() %>%
            which()
          if(length(number)>0){
            final.file$cdna[i] <-  unlist(mut$equivalent_descriptions$c[[number]][1])
            final.file$prot[i] <- unlist(mut$equivalent_descriptions$c[[number]][2])}

        }else{

        }
        print (i)
      }}
  }

  #ADD extra information
  ## Classification
  ## Present or not in both pipelines
  ## Split cdna and prot into different variables
  ## Determine if paralogous and >60
  final.file <- final.file %>%
    dplyr::rowwise() %>%
    dplyr::mutate(class_paralogous= ifelse(ID %in% variants.pseudogen, "paralogous",""),
                  class=ifelse(ID %in% ben, "BEN", "Not found"),
                  class=ifelse(ID %in% vus, "VUS", class),
                  class=ifelse(ID%in% pat, "PAT", class),
                  class=ifelse(ID%in% lpat, "lPAT", class),
                  class=ifelse(ID%in% lben, "lBEN", class),
                  present_pipelines = ifelse(is.na(AF.y),
                                                   "general",
                                                   ifelse(is.na(AF.x),
                                                          "E11",
                                                          "2 pipelines" )),
                  NM= "NM_000535.7",
                  cdna = stringr::str_split(cdna, ":") %>% purrr::map(2) %>% unlist,
                  NP = "NP_000526.2",
                  prot = stringr::str_split(prot, ":") %>% purrr::map(2) %>% unlist,
                  paralogous_above60 = class_paralogous == "paralogous" && AF.x >= 0.6) %>%
    dplyr::relocate(cdna, .after=ID) %>%
    dplyr::relocate(prot, .after=cdna)%>%
    dplyr::relocate(paralogous_above60, .before=class_paralogous) %>%
    dplyr::relocate(present_pipelines, .before=paralogous_above60) %>%
    dplyr::relocate(INROI, .after=class) %>%
    dplyr::relocate(FILTER.x, .after=INROI) %>%
    dplyr::relocate(FILTER.y, .after=FILTER.x)

  return(final.file)
}


filter_variants <- function(all.variants.final){
  all.variants.final2 <- all.variants.final %>%
    dplyr::rowwise() %>%
    dplyr::mutate(what_to_do = ifelse(!FILTER.x %in% c("PASS", "NM5.25", "InIns", "NM5.25;InIns", "InIns;NM5.25;" ),
                                        "QUALITY FILTER NOT PASSED",
                                        ifelse(INROI == "FALSE", "Out of ROI",
                                               ifelse(paralogous_above60 =="",
                                                      "Do not classify, paralogous variant <60",
                                                      ifelse(class %in% c("Not found") ,
                                                             "Classify variant",
                                                             ifelse(class %in% c("BEN", "lBEN"),
                                                                    "Do not report, BEN/lBEN variant",
                                                                    "CHEK E11 pipeline")))))) %>%
    dplyr::mutate(what_to_do = ifelse(what_to_do != "CHEK E11 pipeline",
                                      what_to_do,
                                      ifelse(class_paralogous == "",
                                             ifelse(class %in% c("PAT", "lPAT"),
                                                    ifelse(present_pipelines == "2 pipelines",
                                                           "Perform LR-PCR",
                                                           "Only perform LR-PCR if IHC PMS2-"),
                                                    ifelse(class %in% c("VUS") && present_pipelines == "2 pipelines",
                                                           "Report VUS",
                                                           "Do not report VUS")
                                                          ),
                                                    ifelse(class %in% c("PAT", "lPAT"),
                                                           ifelse(present_pipelines == "2 pipelines",
                                                                  "Perform LR-PCR",
                                                                  "Do not perform LR-PCR"),
                                                           ifelse(present_pipelines == "2 pipelines",
                                                           "Report VUS",
                                                           "Do not report VUS")))))



}
