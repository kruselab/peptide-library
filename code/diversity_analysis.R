##### GENERATE CLONES ##########################################################

# define NNK codons
AA_vec = c("C", "D", "E", "F", "H", "I", "K", "M", "N", "Q", "W", "Y", "A", "A",
           "G", "G", "P", "P", "T", "T", "V", "V", "L", "L", "L", "R", "R", "R",
           "S", "S", "S", "*")
AA_table = setNames(AA_vec, as.character(seq_along(AA_vec)))


# generate a library of size n and length L; save in multiple files such that
# each file contains the unique sequences from a smaller number, chunk_size, of
# generated clones.
build_single_length_library = function(n, 
                                       L = 8, 
                                       chunk_size = 1e6, 
                                       output_dir = "chunks") {
  if (!dir.exists(output_dir)) {dir.create(output_dir)}
  total_chunks = ceiling(n / chunk_size)
  set.seed(0)
  for (i in 1:total_chunks) {
    cat("Generating chunk", i, "of", total_chunks, "\n")
    current_chunk_size = min(chunk_size, n - (i - 1) * chunk_size)
    codons = sample(x = 31, size = current_chunk_size * L, replace = TRUE)
    AA_matrix = matrix(AA_table[codons], nrow = current_chunk_size, ncol = L)
    sequences = do.call(paste0, as.data.frame(AA_matrix))
    unique_sequences = sort(unique(sequences))
    output_file = file.path(output_dir, paste0("chunk_", i, "_unique_sequences.txt"))
    data.table::fwrite(list(unique_sequences), file = output_file)
  }
  cat("All chunks generated.\n")
  return(invisible())
}


##### MERGE SORT AND COUNT #####################################################

# get the target file based on the first three characters
get_target_file = function(sequence, sort_dir = "sorted_files") {
  first_three_chars = substr(sequence, 1, 3)
  target_file = file.path(sort_dir, paste0("sorted_", first_three_chars, ".txt"))
  return(target_file)
}


# append sequences to appropriate three-character files in bulk
append_sequences_to_files = function(sequences, sort_dir = "sorted_files") {
  if (!dir.exists(sort_dir)) {dir.create(sort_dir)}
  sequence_list = split(sequences, substr(sequences, 1, 3))
  for (first_three_chars in names(sequence_list)) {
    target_file = file.path(sort_dir, paste0("sorted_", first_three_chars, ".txt"))
    data.table::fwrite(list(sequence_list[[first_three_chars]]), 
                       file = target_file, 
                       sep = "\n", 
                       append = TRUE)
  }
  return(invisible())
}


# process sequence chunk into three-character files
process_sequence_chunk = function(chunk_file, sort_dir = "sorted_files") {
  sequences = data.table::fread(file = chunk_file, header = FALSE)[[1]]
  append_sequences_to_files(sequences, sort_dir)
  return(invisible())
}


# process all chunks in dir
process_chunk_dir = function(dir = "chunks", sort_dir = "sorted_files") {
  files = list.files(dir, pattern = "unique_sequences.txt$", full.names = TRUE)
  for (i in 1:length(files)) {
    cat("Processing file", i, "of", length(files), ":", basename(files[i]), "\n")
    process_sequence_chunk(files[i])
  }
  cat("All chunk files processed!\n")
  return(invisible())
}


# count unique sequences across all three-character files
count_total_uniques <- function(dir = "sorted_files") {
  files = list.files(dir, pattern = "sorted_.*\\.txt$", full.names = TRUE)
  total_unique_count = 0
  for (i in 1:length(files)) {
    cat("Processing file", i, "of", length(files), ":", basename(files[i]), "\n")
    sequences = data.table::fread(file = files[i], header = FALSE)[[1]]
    total_unique_count = total_unique_count + length(unique(sequences))
    cat(length(unique(sequences)), "uniques in file,", total_unique_count, "total.\n")
  }
  return(total_unique_count)
}


##### RUN ######################################################################

# # example run
# build_single_length_library(1e8)
# process_chunk_dir()
# count_total_uniques()
