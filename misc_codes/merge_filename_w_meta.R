## usage:
## input 1: sample_filename_map.txt (a file name mapping file)
## input 2: rk1.txt_fixed.txt (experimental meta data)

args <- commandArgs(TRUE)

fn_map <- read.delim(args[1], header=F)
names(fn_map) <- c("original_sample_names", "SAMPLES")

meta <- read.delim(args[2])

dim(fn_map)
dim(meta)

meta.new <- merge(meta, fn_map, by.x="CIBNOR_id", by.y="original_sample_names")
