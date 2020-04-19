# filename = "CurOvGradeKEGGnets.Rdata"
cmdArgs <- commandArgs(trailingOnly = FALSE)
filename = commandArgs(trailingOnly=TRUE)
# Rscript --vanilla read_cur_data.R CurOvGradeKEGGnets.RData

# must install optparse
# library("optparse")
# 
# option_list = list(
#   make_option(c("-f", "--file"), type="character", default=NULL, 
#               help="dataset file name", metavar="character"),
# ); 
# 
# filename = opt$file
# 
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);

needle <- "--file="
match <- grep(needle, cmdArgs)
path = ""
if (length(match) > 0) {
  # Rscript
  path = normalizePath(sub(needle, "", cmdArgs[match]))
} else {
  # 'source'd via R console
  path = normalizePath(sys.frames()[[1]]$ofile)
}

path = substr(path, 1, nchar(path) - 15)
path_out = paste(path, "files/", sep="")
data = environment()
val = get(load(filename, data))
count = 0
for (i in ls(val)){
  a <- get(i, val)
  expr <- get("expr", a)
  write.csv(expr, file = paste(path_out, "expr", count, ".csv", sep = ""))
  grade <- get("grade", a)
  write.csv(grade, file = paste(path_out, "grade", count, ".csv", sep = ""))
  count = count + 1
}


