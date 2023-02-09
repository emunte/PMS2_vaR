# Description: Run the PMS2_vaR pipeline
# USAGE: Rscript runPMS2_vaR.R [-t tools_file] [-d datasets_file]

#libs


# Check runPMS2_vaR.R is being called from PMS2_vaR folder

if (length(list.files(pattern = "runPMS2_vaR.R"))== 0){
  cat("Sorry, runPMS2_vaR.R .R should be called from PMS2_vaR.R  folder\n")
  quit()
}


