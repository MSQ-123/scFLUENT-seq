
library(data.table)
Args <- commandArgs(TRUE)

R1_bed12 <- fread(Args[1])

# split into R1 and R2

R1_bed12$id <- unlist(lapply(strsplit(R1_bed12$V4, split = '/'), FUN = function(x)x[[1]]))

R1 <- R1_bed12[grep(pattern = '/1', x = R1_bed12$V4), ]
R2 <- R1_bed12[grep(pattern = '/2', x = R1_bed12$V4), ]

R2 <- R2[match(R1$id, R2$id),]
print(table(R2$V6 == R1$V6))
print(table(R2$V1 == R1$V1))

# pseudo-long reads # for junction reads it will record more than 1 blocks, so replace the block num to 1

idx <- R2$V6 == '+'
R2$start <- R2$V2
R2$end <- R2$V3
R2[idx, ]$end <- R1[idx, ]$V3

idx <- R2$V6 == '-'
R2[idx, ]$start <- R1[idx, ]$V2

R2$V2 <- R2$start
R2$V3 <- R2$end

R2$start <- NULL
R2$id <- NULL
R2$end <- NULL

R2$V7 <- R2$V2
R2$V8 <- R2$V3
R2$V11 <- R2$V3 - R2$V2


#replace the block num to 1, ignore splicing
R2$V10 <- 1
R2$V12 <- 0


write.table(R2, col.names = F, row.names = F, quote = F, sep = "\t",
            file = Args[2])
