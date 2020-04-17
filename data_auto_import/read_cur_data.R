# filename = "CurOvGradeKEGGnets.Rdata"
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

data = environment()
val = get(load(filename, data))
count = 0
for (i in ls(val)){
  a <- get(i, val)
  expr <- get("expr", a)
  write.csv(expr, file = paste("expr", count, ".csv", sep = ""))
  grade <- get("grade", a)
  write.csv(grade, file = paste("grade", count, ".csv", sep = ""))
  count = count + 1
}


