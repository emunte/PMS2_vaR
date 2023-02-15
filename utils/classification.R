### Read variant classifications

vars.class <- read.csv("./data/variants_classificated.csv")

# Classify variants into groups
ben <- vars.class %>%
  dplyr::filter(classification %in% c("BEN")) %>%
  dplyr::mutate(ID.vars = stringr::str_replace(ID.vars, " ", ""))

lben <- vars.class %>%
  dplyr::filter(classification %in% c("LBEN")) %>%
  dplyr::mutate(ID.vars = stringr::str_replace(ID.vars, " ", ""))

pat <- vars.class %>%
  dplyr::filter(classification %in% c("PAT")) %>%
  dplyr::mutate(ID.vars = stringr::str_replace(ID.vars, " ", ""))
lpat <- vars.class %>%
  dplyr::filter(classification %in% c("LPAT")) %>%
  dplyr::mutate(ID.vars = stringr::str_replace(ID.vars, " ", ""))

vus <- vars.class %>%
  dplyr::filter(classification %in% c("VUS")) %>%
  dplyr::mutate(ID.vars = stringr::str_replace(ID.vars, " ", ""))


# Paralogous variant list
vars.paralogous <- read.csv("./data/variants_pms2CL.csv") %>% as.matrix()

# BEd
bed.file <-  read.csv2("./data/PMS2_bed_file.bed", sep="\t", header = TRUE)
bedGR <- regioneR::toGRanges(bed.file)
names(bedGR) <- paste0("E", 15:1)


