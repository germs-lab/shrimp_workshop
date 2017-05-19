## usage:
## input 1: sample_filename_map_suffix_rmed.txt (a file name mapping file)
## input 2: rk1.txt_fixed.txt (experimental meta data)
## input 3: output file location and name (eg., rk1_meta_sample_data.txt)

args <- commandArgs(TRUE)

fn_map <- read.delim(args[1], header=F)
names(fn_map) <- c("original_sample_names", "SAMPLES")

meta <- read.delim(args[2])
dim(meta)

meta.new <- merge(meta, fn_map, by.x="CIBNOR_id", by.y="original_sample_names")
dim(meta.new)

write.table(meta.new, args[3], sep="\t", row.names=F, quote=F)
