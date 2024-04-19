### Read variant classifications

vars.class <- read.csv("./data/variants_classified.csv")

# Classify variants into groups
ben <- vars.class %>%
  dplyr::filter(classification %in% c("BEN")) %>%
  dplyr::mutate(ID.vars = stringr::str_replace(ID.vars, " ", ""),
                ID.vars.hg38= stringr::str_replace(ID.vars.hg38, " ", ""))

lben <- vars.class %>%
  dplyr::filter(classification %in% c("LBEN")) %>%
  dplyr::mutate(ID.vars = stringr::str_replace(ID.vars, " ", ""),
                ID.vars.hg38= stringr::str_replace(ID.vars.hg38, " ", ""))

pat <- vars.class %>%
  dplyr::filter(classification %in% c("PAT")) %>%
  dplyr::mutate(ID.vars = stringr::str_replace(ID.vars, " ", ""),
                ID.vars.hg38= stringr::str_replace(ID.vars.hg38, " ", ""))
lpat <- vars.class %>%
  dplyr::filter(classification %in% c("LPAT")) %>%
  dplyr::mutate(ID.vars = stringr::str_replace(ID.vars, " ", ""),
                ID.vars.hg38= stringr::str_replace(ID.vars.hg38, " ", ""))

vus <- vars.class %>%
  dplyr::filter(classification %in% c("VUS")) %>%
  dplyr::mutate(ID.vars = stringr::str_replace(ID.vars, " ", ""),
                ID.vars.hg38= stringr::str_replace(ID.vars.hg38, " ", ""))


classification <- list(ben = ben,
                       lben = lben,
                       vus = vus,
                       lpat = lpat,
                       pat = pat)


# Paralogous variant list
vars.paralogous <- read.csv("./data/variants_pms2CL.csv", header = TRUE, sep="\t")

# BEd
bed.file <-  read.csv2(ifelse(args$genome == "hg19", "./data/PMS2_bed_file.bed", "./data/PMS2_bed_file_hg38.bed"), sep="\t", header = TRUE)
bedGR <- regioneR::toGRanges(bed.file)
names(bedGR) <- paste0("E", 15:1)


