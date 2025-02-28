##### STANDARD ANALYSIS FUNCTIONS ##############################################

# packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggseqlogo))
suppressPackageStartupMessages(library(stringdist))


# generate list of common directories and make them if they don't already exist
setup_wd = function(wd) {
  for (i in 1:length(wd)) {
    dir.create(path = wd[[i]])
  }
  return(invisible())
}


# load merged reads fastq file
load_merged_reads = function(file, format = "fastq") {
  reads = readDNAStringSet(filepath = file, use.names = F, format = format)
  return(reads)
}


# from a DNAStringSet, return those that match F and R primer sequences
filter_reads_by_primers = function(reads, library, F_seq = "", R_seq = "", remove_flanks = T) {

  # choose primer strings based on library
  if (library == "peptide") {
    F_primer = "CAAGGAGTTTTT"
    R_primer = "TACACATACTTT"
  } else {
    F_primer = F_seq
    R_primer = R_seq
  }
  
  # find reads that match primers (either forward or reverse complement)
  matches = c(reads[vcountPattern(F_primer, reads) > 0 & vcountPattern(R_primer, reads) > 0],
              reverseComplement(reads)[vcountPattern(F_primer, reverseComplement(reads)) > 0 & vcountPattern(R_primer, reverseComplement(reads)) > 0])
  
  # remove any flanking characters outside primers
  if (remove_flanks) {
    matches = DNAStringSet(sub(paste0(".*(", F_primer, ".*", R_primer, ").*"), "\\1", matches))
  }
  return(matches)
}


# translate reads
translate_reads = function(reads) {
  translations = translate(reads, no.init.codon = T, if.fuzzy.codon = "solve")
  return(translations)
}


# from a AAStringSet, return those that match library-specific regex pattern
filter_translations_by_library_design = function(translations, library, pattern = "") {
  
  # choose regex pattern based on library
  if (library == "peptide") {
    pattern = "QGVFEYYKSVTFVSNCGSHPSTTSKGSPINTQYVFKLLQA[GS]*YPYDVPDYAGGGGS\\w{0,8}\\*"
  } else {
    pattern = pattern
  }

  # filter by regex
  translations = translations[grepl(pattern, translations)]
  return(translations)
}


# convert StringSet to df
convert_StringSet_to_df = function(stringset, combine_duplicates = T) {
  sequence_df = as_tibble(as.character(stringset))
  names(sequence_df) = "sequence"
  if (combine_duplicates) {
    sequence_df = sequence_df %>%
      group_by(sequence) %>%
      summarise(read_count = n()) %>%
      ungroup() %>%
      arrange(desc(read_count), sequence)
  } else {
    sequence_df = sequence_df %>%
      mutate(read_count = 1)
  }
  return(sequence_df)
}


# extract peptide sequences (remove Aga2p/linkers/etc.)
extract_peptides = function(sequence_df, pattern = ".*QGVFEYYKSVTFVSNCGSHPSTTSKGSPINTQYVFKLLQA[GS]*YPYDVPDYAGGGGS(\\w{0,8})\\*.*") {
  sequence_df$peptide_sequence = gsub(pattern, "\\1", sequence_df$sequence)
  return(sequence_df)
}


# normalize abundance to reads per million
normalize_depth = function(data) {
  data = data %>%
    mutate(reads_per_million = round(read_count * (1e6 / sum(read_count)), 3))
  return(data)
}


# plot sequence logos
plot_peptide_sequence_logos = function(data, returned_alignment = "N") {
  data_processed = data %>% 
    mutate(NT_aligned_peptide = paste0(peptide_sequence, strrep("-", 8 - nchar(peptide_sequence)))) %>%
    mutate(CT_aligned_peptide = paste0(strrep("-", 8 - nchar(peptide_sequence)), peptide_sequence)) %>%
    uncount(read_count)
  NT_sequence_logo = ggplot() + 
    geom_logo(data_processed$NT_aligned_peptide, 
              seq_type = "aa", 
              method = "prob",
              col_scheme = make_col_scheme(chars = c("X"), values = 1)) +
    scale_x_continuous("Position (N-aligned)", breaks = 1:8) +
    scale_y_continuous(breaks = seq(0,1,.2)) +
    theme_logo() + 
    theme(axis.ticks.y = element_line())
  CT_sequence_logo = ggplot() + 
    geom_logo(data_processed$CT_aligned_peptide, 
              seq_type = "aa", 
              method = "prob",
              col_scheme = make_col_scheme(chars = c("X"), values = 1)) +
    scale_x_continuous("Position (C-aligned)", breaks = 1:8) +
    scale_y_continuous(breaks = seq(0,1,.2)) +
    theme_logo() + 
    theme(axis.ticks.y = element_line())
  ggsave("output/NT_aligned_sequence_logo.pdf", NT_sequence_logo, height = 2, width = 4)
  ggsave("output/CT_aligned_sequence_logo.pdf", CT_sequence_logo, height = 2, width = 4)
  if (returned_alignment == "N") {return(NT_sequence_logo)}
  else if (returned_alignment == "C") {return(CT_sequence_logo)}
  else {
    warning("Invalid returned_alignment argument. Returning NT alignment by default.", call. = F)
    return(NT_sequence_logo)
  }
}


# run standard analysis pipeline; generate peptide list and sequence logos
run_pipeline = function(merged_reads_file) {
  setup_wd(list(output = "output/"))
  all_reads = load_merged_reads(merged_reads_file)
  matching_reads = all_reads %>% filter_reads_by_primers(library = "peptide")
  peptides = matching_reads %>%
    translate_reads() %>%
    filter_translations_by_library_design(library = "peptide") %>%
    convert_StringSet_to_df() %>%
    extract_peptides() %>%
    group_by(peptide_sequence) %>%
    summarise(read_count = sum(read_count)) %>%
    arrange(desc(read_count)) %>%
    normalize_depth()
  write_csv(peptides, "output/peptides.csv")
  plot_peptide_sequence_logos(peptides)
  return(invisible())
}


##### ADDITIONAL ANALYSIS FUNCTIONS ############################################

# calculate and plot Levenshtein edit distances from a given peptide
plot_edit_distances = function(data, peptide, save_plot = TRUE, save_csv = TRUE, prefix = "") {
  data_processed = data %>%
    mutate(edit_distance = stringdist(peptide, peptide_sequence, method = "lv")) %>%
    uncount(read_count) %>%
    group_by(edit_distance) %>%
    summarize(count = n()) %>%
    merge(expand_grid(edit_distance = 0:8, count0 = 0), all.y = TRUE) %>%
    mutate(count = ifelse(!is.na(count), count, count0)) %>%
    select(-count0) %>%
    mutate(frequency = round(count / sum(count), 5))
  plot_edit_dist = ggplot(data_processed) +
    geom_line(aes(x = edit_distance, y = frequency)) +
    scale_x_continuous(paste0("Levenshtein edit distance from \"", peptide, "\""), limits = c(0,8), breaks = 0:8) +
    scale_y_continuous("Frequency", limits = c(0,1), breaks = seq(0, 1, 0.2)) +
    theme_cowplot()
  if (save_plot) {
    ggsave(paste0(prefix, peptide, "_edit_distance_plot.pdf"), plot_edit_dist, width = 5, height = 4)
  }
  if (save_csv) {
    write_csv(data_processed, paste0(prefix, peptide, "_edit_distances.csv"))
  }
  return(plot_edit_dist)
}