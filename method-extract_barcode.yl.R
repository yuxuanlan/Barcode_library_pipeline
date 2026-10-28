# modified version of cellbarcode to use in Barcode_library_pipeline
#  old name: method-extract_barcode.yl1.R
library(CellBarcode)
library(magrittr)
library(stringr)
# install.packages("data.table") # for testing cellbarcode
library(data.table)
library(Biostrings)

ifnullelse <- function(cand1, other) {
  if (is.null(cand1)) {
    return(other)
  } else {
    return(cand1)
  }
}

bc_process_sample_name <- function(sample_name, metadata, input_names) {
  # input metadata,
  # input samle_name,
  # names of input data or the file name
  # metadata = sample_name > names or file name
  if (is.null(sample_name)) {
    if (!is.null(metadata)) {
      sample_name <- rownames(metadata)
    } else {
      sample_name <- input_names
    }
  } else {
    if (length(sample_name) != length(input_names)) {
      stop("Sample name length does not meet sample number.")
    }
  }
  
  if (!is.null(metadata)) {
    if (!all(sample_name == rownames(metadata))) {
      stop("Sample name does not match row name of metadata.")
    }
  }
  sample_name
}

bc_process_metadata <- function(sample_name, old_metadata, new_metadata) {
  res <- merge(old_metadata, new_metadata, by = 0, all=TRUE)
  res <- as.data.frame(res)
  rownames(res) <- res$Row.names
  res$Row.names <- NULL
  res[sample_name, ]
}

bc_extract_metadata <- function(x, sample_name) {
  d <- vapply(x, function(x_i) {
    res <- c(attr(x_i, "raw_read_count"), attr(x_i, "barcode_read_count"))
    res
  }, c("raw_read_count" = 1, "barcode_read_count" = 1))
  d <- as.data.frame(t(d))
  rownames(d) <- sample_name
  d
}

bc_extract_df2 <- function(x, 
                           pattern="", 
                           maxLDist=0,
                           costs = list(sub = 1, ins = 99, del = 99), 
                           ordered = TRUE) {
  # function by Yuxuan Lan based on CellBarcode::bc_extract()
  # input barcode_df, with column freq and seq
  # input pattern for matching, allowing regexp
  # output match_result:
  #
  # matched_list.fwd_gfp$d[1:2,]
  #               barcode_seq count
  #                    <char> <int>
  # 1: CGGCATGGACGAGCTGTACCAG    17
  # 2: CGGCATGGACGAGCTGTACCAT     6
  #   matched_list.fwd_gfp$reads[1:2,]
  #   seq
  # 1 GTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCG...
  # 2 GCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTC...
  #   freq                     bc
  # 1    3 CGGCATGGACGAGCTGTACCAG
  # 2    2 CGGCATGGACGAGCTGTACCAT
  
  # sequence frequecy
  reads_freq <- x$freq
  # sequence
  reads_seq <- x$seq
  
  match_result <- list()
  # if maxLDist > 0, use utils::aregexec to match the barcode else use
  # stringr::str_match to perform the barcode match. The stringr::str_match is
  # faster
  if (maxLDist > 0) {
    # CURRENTLY NOT WORKING ON THIS SECTION SO CODE RELATED TO UMI COULD BE WRONG
    # regular expression match
    m <- utils::aregexec(
      pattern, 
      reads_seq, 
      fixed = FALSE, 
      max.distance = maxLDist, 
      costs = costs
    )
    
    # captured sequence
    seq_l <- regmatches(reads_seq, m)
    
    # matched sequence
    seq_v <- plyr::laply(seq_l, function(x) { x[1] })
    # captured barcode
    barcode_v <- plyr::laply(seq_l, function(x) { x[1] })
    
    # dataframe listing reads matching the pattern
    match_idx <- which(!is.na(seq_v))
    df.matched_reads <- data.frame(seq=x$seq[match_idx], freq=x$freq[match_idx], bc=seq_v[match_idx])
    
  } else {
    # regular expression match
    m <- stringr::str_match(reads_seq, pattern)
    #      [,1]                                                
    # [1,] "GCTCTCCCAGATAGTTAAGACTACCCTCCCTAGGAACGACTATGAGCGGT"
    # [2,] "GCTCGAACAGTGGGTTTGAACTTGACTCACAAGGGAGGACCTAGAGCGGT"
    # [3,] NA                                                  
    #      [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8] 
    # [1,] "TCC" "ATA" "AAG" "ACC" "CCT" "AAC" "TAT"
    # [2,] "GAA" "TGG" "TGA" "TGA" "ACA" "GAG" "CTA"
    # [3,] NA    NA    NA    NA    NA    NA    NA   
    # capture sequence
    seq_v <- m[, 1]
    
    # matched barcode
    # barcode_v <- m[, pattern_type['barcode'] + 1] # original
    barcode_v <- m[, 1]
    # [1] "GCTCTCCCAGATAGTTAAGACTACCCTCCCTAGGAACGACTATGAGCGGT"
    # [2] "GCTCGAACAGTGGGTTTGAACTTGACTCACAAGGGAGGACCTAGAGCGGT"
    # [3] NA 
    
    # dataframe listing reads matching the pattern
    match_idx <- which(!is.na(m[,1]))
    df.matched_reads <- data.frame(seq=x$seq[match_idx], freq=x$freq[match_idx], bc=m[match_idx,1])
  }
  # use this when UMI is not in pattern_type
  d <- data.table(
    barcode_seq = barcode_v, 
    count = reads_freq
  )[, .(count = sum(count)), by = barcode_seq]
  
  # order the data by counts
  if (ordered) {
    d <- d[order(count, decreasing = TRUE)]
  } 
  #                                           barcode_seq count
  #                                                <char> <num>
  # 1: GCTCTCCCAGATAGTTAAGACTACCCTCCCTAGGAACGACTATGAGCGGT    20
  # 2: GCTCGAACAGTGGGTTTGAACTTGACTCACAAGGGAGGACCTAGAGCGGT    20
  # 3:                                               <NA>    20
  d <- stats::na.omit(d)
  
  # prepare the row_read_count & barcode_read_count
  attr(d, "raw_read_count") <- sum(reads_freq) # 60
  attr(d, "barcode_read_count") <- sum(d$count) # 40
  
  # d
  #                                           barcode_seq count
  #                                                <char> <num>
  # 1: GCTCTCCCAGATAGTTAAGACTACCCTCCCTAGGAACGACTATGAGCGGT    20
  # 2: GCTCGAACAGTGGGTTTGAACTTGACTCACAAGGGAGGACCTAGAGCGGT    20
  match_result$d <- d
  match_result$reads <- df.matched_reads
  
  return(match_result)
}

bc_extract_df <- function(x, 
                          pattern = "", 
                          sample_name = NULL,
                          metadata = NULL, 
                          maxLDist = 0, 
                          pattern_type = c(barcode = 1), 
                          costs = list(sub = 1, ins = 99, del = 99), 
                          ordered = TRUE) {
  
  # sequence frequecy
  reads_freq <- x$freq
  # sequence
  reads_seq <- x$seq
  
  # if maxLDist > 0, use utils::aregexec to match the barcode else use
  # stringr::str_match to perform the barcode match. The stringr::str_match is
  # faster
  if (maxLDist > 0) {
    # CURRENTLY NOT WORKING ON THIS SECTION SO CODE RELATED TO UMI COULD BE WRONG
    # regular expression match
    m <- utils::aregexec(
      pattern, 
      reads_seq, 
      fixed = FALSE, 
      max.distance = maxLDist, 
      costs = costs
    )
    
    # captured sequence
    seq_l <- regmatches(reads_seq, m)
    
    # matched sequence
    seq_v <- plyr::laply(seq_l, function(x) { x[1] })
    # captured barcode
    barcode_v <- plyr::laply(seq_l, function(x) { x[1] })
  } else {
    # regular expression match
    m <- stringr::str_match(reads_seq, pattern)
    #      [,1]                                                
    # [1,] "GCTCTCCCAGATAGTTAAGACTACCCTCCCTAGGAACGACTATGAGCGGT"
    # [2,] "GCTCGAACAGTGGGTTTGAACTTGACTCACAAGGGAGGACCTAGAGCGGT"
    # [3,] NA                                                  
    #      [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8] 
    # [1,] "TCC" "ATA" "AAG" "ACC" "CCT" "AAC" "TAT"
    # [2,] "GAA" "TGG" "TGA" "TGA" "ACA" "GAG" "CTA"
    # [3,] NA    NA    NA    NA    NA    NA    NA   
    # capture sequence
    seq_v <- m[, 1]
    
    # matched barcode
    # barcode_v <- m[, pattern_type['barcode'] + 1] # original
    barcode_v <- m[, 1]
    # [1] "GCTCTCCCAGATAGTTAAGACTACCCTCCCTAGGAACGACTATGAGCGGT"
    # [2] "GCTCGAACAGTGGGTTTGAACTTGACTCACAAGGGAGGACCTAGAGCGGT"
    # [3] NA 
  }
  # use this when UMI is not in pattern_type
  d <- data.table(
    barcode_seq = barcode_v, 
    count = reads_freq
  )[, .(count = sum(count)), by = barcode_seq]
  
  # order the data by counts
  if (ordered) {
    d <- d[order(count, decreasing = TRUE)]
  } 
  #                                           barcode_seq count
  #                                                <char> <num>
  # 1: GCTCTCCCAGATAGTTAAGACTACCCTCCCTAGGAACGACTATGAGCGGT    20
  # 2: GCTCGAACAGTGGGTTTGAACTTGACTCACAAGGGAGGACCTAGAGCGGT    20
  # 3:                                               <NA>    20
  d <- stats::na.omit(d)
  
  # prepare the row_read_count & barcode_read_count
  attr(d, "raw_read_count") <- sum(reads_freq) # 60
  attr(d, "barcode_read_count") <- sum(d$count) # 40
  
  d
  #                                           barcode_seq count
  #                                                <char> <num>
  # 1: GCTCTCCCAGATAGTTAAGACTACCCTCCCTAGGAACGACTATGAGCGGT    20
  # 2: GCTCGAACAGTGGGTTTGAACTTGACTCACAAGGGAGGACCTAGAGCGGT    20
  return(d)
}
check_sample_name <- function(object) {
  metaname <- bc_names(object)
  
  # metadata v.s. messyBc
  if (all(metaname != names(object@messyBc)))
    stop("The messyBc sample names are not consistent with metadata")
  
  # metadata v.s. cleanBc
  if (!is.null(object@cleanBc)) {
    if (all(metaname != names(object@cleanBc)))
      stop("The cleanBc sample names are not consistent with metadata")
  }
}
BarcodeObj <- setClass("BarcodeObj", 
                       slots=list(
                         metadata="data.frame",
                         messyBc="list",
                         cleanBc="listOrNull"
                       ),
                       prototype=list(
                         metadata=data.frame(),
                         messyBc=list(),
                         cleanBc=NULL
                       ),
                       validity=check_sample_name
)

is_pure_dataframe <- function(obj) {
  inherits(obj, "data.frame") && !inherits(obj, "tbl") && !inherits(obj, "data.table")
}
read_fastq <- function(in_file_path) {
  .Call('_CellBarcode_read_fastq', PACKAGE = 'CellBarcode', in_file_path)
}
read_fastq_gz <- function(in_file_path) {
  .Call('_CellBarcode_read_fastq_gz', PACKAGE = 'CellBarcode', in_file_path)
}

# function from : setMethod("bc_extract", "character", function
bc_extract_fastq <- function(
    x,
    direction='fwd', # YL
    pattern,
    sample_name=NULL,
    maxLDist = 0,
    metadata = NULL,
    pattern_type=c(barcode=1), 
    ordered = TRUE){
  if (!is.null(metadata) & !is_pure_dataframe(metadata)) {
    stop("metadata should be a data.frame.")
  }
  
  ## check input file
  if (!all(file.exists(x))) {
    stop(paste0("Input file(s) do not exist."))
  }
  
  input_names <- basename(x)
  # [1] "plate1_merged.fastq.gz"
  
  sample_name <- bc_process_sample_name(sample_name, metadata, input_names)
  # [1] "plate1_merged.fastq.gz"
  
  pattern = rep(pattern, length.out = length(x))
  maxLDist = rep(maxLDist, length.out = length(x))
  
  # how many reads in fastq.gz?
  # barcode_df = read_fastq_gz(x[1])
  # dim(barcode_df)
  # [1] 599983      2
  messyBc <- lapply(seq_along(x), function(i) {
    if (grepl(".fq$", x[i]) | grepl(".fastq$", x[i])) {
      barcode_df = read_fastq(x[i])
    } else if (grepl(".fq.gz$", x[i]) | grepl(".fastq.gz$", x[i])) {
      barcode_df = read_fastq_gz(x[i])
    } else {
      stop("The input is not Fastq file. Please check the input.")
    }
    # YL:
    if (direction=='rc'){
      rc <- as.character(reverseComplement(DNAStringSet(x=barcode_df$seq)))
      barcode_df$seq= rc
    }
    bc_extract_df(
      barcode_df,
      sample_name = sample_name[i], 
      pattern = pattern[i], 
      maxLDist = maxLDist[i], 
      pattern_type = pattern_type, 
      costs = costs, 
      ordered = ordered)
  })
  # class(messyBc)
  # [1] "list"
  
  # dim(messyBc[[1]])
  # [1] 23858     2
  
  # messyBc
  # [[1]]
  #                                               barcode_seq count
  #                                                    <char> <int>
  #     1: GCTCCGGCAGAGCGTTTGCACTGACCTCTTGAGGATGGACAGGGAGCGGT 20855
  #     2: GCTCGGACAGTTAGTTGACACTGTGCTCAATAGGGTTGACGTTGAGCGGT 19228
  #     3: GCTCTTACAGCCGGTTGAAACTAACCTCAAAAGGATGGACGGCGAGCGGT 18202
  #     4: GCTCATACAGTAAGTTGCTACTTGACTCCCTAGGTAAGACTACGAGCGGT 18153
  #     5: GCTCTAGCAGCTTGTTTCTACTTCACTCCATAGGTAAGACCAAGAGCGGT 17623
  #     
  new_metadata <- bc_extract_metadata(messyBc, sample_name)
  metadata <- bc_process_metadata(sample_name, metadata, new_metadata)
  
  names(messyBc) <- sample_name
  # output <- list(messyBc = messyBc, metadata = metadata)
  # class(output) <- "BarcodeObj"
  
  # check if no barcodes found
  if (all(sapply(messyBc, nrow) == 0)) {
    message("No barcode found. Please check the input file and pattern.")
  }
  output <- BarcodeObj(metadata=metadata, messyBc=messyBc)
  # > output@metadata
  #                        raw_read_count barcode_read_count
  # plate1_merged.fastq.gz        5603546            2518821
  
  # output@messyBc
  # $plate1_merged.fastq.gz
  #                                               barcode_seq count
  #                                                    <char> <int>
  #     1: GCTCCGGCAGAGCGTTTGCACTGACCTCTTGAGGATGGACAGGGAGCGGT 20855
  #     2: GCTCGGACAGTTAGTTGACACTGTGCTCAATAGGGTTGACGTTGAGCGGT 19228
  #     3: GCTCTTACAGCCGGTTGAAACTAACCTCAAAAGGATGGACGGCGAGCGGT 18202
  #     4: GCTCATACAGTAAGTTGCTACTTGACTCCCTAGGTAAGACTACGAGCGGT 18153
  #     5: GCTCTAGCAGCTTGTTTCTACTTCACTCCATAGGTAAGACCAAGAGCGGT 17623
  return(output)
}


