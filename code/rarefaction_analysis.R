# packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(immunarch))
suppressPackageStartupMessages(library(cowplot))


# converts NGS pipeline output csv file to immunarch-readable file
convert_NGS_output_to_immunarch = function(input_file,
                                           output_file = gsub("(.*)\\.csv", 
                                                              "\\1_immunarch.csv", 
                                                              input_file)) {
  input_data = read_csv(input_file,
                        na = character(), # ensures "" and "NA" peptides aren't interpretted as missing or converted to NAs
                        show_col_types = FALSE) %>%
    mutate(peptide_sequence = ifelse(peptide_sequence == "", # handles specific "" edge case
                                     "EMPTYPEPTIDE",
                                     peptide_sequence)) %>%
    mutate(peptide_sequence = ifelse(peptide_sequence == "NA", # handles specific "NA" edge case
                                     "NAPEPTIDE", 
                                     peptide_sequence))
  output_data = data.frame(CDR3.aa = input_data$peptide_sequence,
                           Clones = input_data$read_count,
                           Proportion = NA, 
                           CDR3.nt = NA, 
                           V.name = NA, 
                           D.name = NA, 
                           J.name = NA, 
                           V.end = NA, 
                           D.start = NA, 
                           D.end = NA, 
                           J.start = NA, 
                           VJ.ins = NA, 
                           VD.ins = NA, 
                           DJ.in = NA, 
                           Sequence = NA)
  write_csv(output_data, output_file)
  return(invisible())
}


# convert all files in dir to immunarch-readable format in output_dir
convert_NGS_outputs_in_dir = function(input_dir, output_dir) {
  dir.create(path = output_dir)
  files = list.files(input_dir, full.names = TRUE)
  for (file in files) {
    output_path = paste0(output_dir, basename(file))
    convert_NGS_output_to_immunarch(file, output_file = output_path)
  }
  return(invisible())
}


# calculate reference samples (i.e., total NGS reads and number of unique clones for each sample)
calc_ref_samples = function(data) {
  ref_samples = data.frame(sample = as.vector(data$meta),
                           n = NA,
                           uniques = NA)
  for (s in ref_samples$Sample) {
    d = as_tibble(data$data[[s]])
    ref_samples[ref_samples$Sample == s, "n"] = sum(d$Clones)
    ref_samples[ref_samples$Sample == s, "uniques"] = nrow(d)
  }
  return(ref_samples)
}


# calculate Chao1 nonparametric asymptotic estimators
calc_chao1 = function(data) {
  chao1 = as.data.frame(repDiversity(data$data, 
                                     .do.norm = TRUE, 
                                     .method = "chao1"))
  chao1 = mutate(chao1, Sample = rownames(chao1))
  return(chao1)
}


# calculate rarefaction-extrapolation curves
calc_rarefaction_curves = function(data, ref_samples) {
  rarefaction_curves = repDiversity(data$data,
                                    .do.norm = TRUE,
                                    .method = "raref", 
                                    .norm = FALSE, 
                                    .extrapolation = 1e6, 
                                    .verbose = FALSE) %>%
    left_join(ref_samples, by = "Sample") %>%
    filter(Size <= 2 * n)
  return(rarefaction_curves)
}


# plot rarefaction curves with Chao1 estimators
plot_rarefaction_curves = function(chao1_estimators, 
                                   rarefaction_curves, 
                                   ref_samples,
                                   save_plot = TRUE,
                                   prefix = "") {
  rarefaction_plot = ggplot() +
    geom_hline(data = chao1_estimators, 
               aes(yintercept = Estimator, color = Sample),
               linetype = "dotted") +
    geom_text(data = chao1_estimators, 
              aes(y = Estimator, label = round(Estimator, 0)),
              x = max(rarefaction_curves$Size), hjust = 1, vjust = -0.2) +
    geom_line(data = rarefaction_curves %>% filter(Type == "interpolation"),
              aes(x = Size, y = Mean, group = Sample, color = Sample)) +
    geom_line(data = rarefaction_curves %>% filter(Type == "extrapolation"),
              aes(x = Size, y = Mean, group = Sample, color = Sample),
              linetype = "longdash") +
    geom_point(data = ref_samples, 
               aes(x = n, y = uniques, color = Sample),
               size = 2, stroke = 0) +
    scale_x_continuous("NGS reads") +
    scale_y_continuous("Unique clones") +
    scale_color_brewer(palette = "Blues") +
    guides(label = "none") +
    theme_cowplot()
  if (save_plot) {
    ggsave(paste0(prefix, "rarefaction_plot.pdf"),
           rarefaction_plot,
           height = 6, width = 7)
  }
  return(rarefaction_plot)
}


# run full rarefaction analysis
run_rarefaction_analysis = function(NGS_outputs_dir, 
                                    immunarch_data_dir = "immunarch/", 
                                    prefix = "") {
  convert_NGS_outputs_in_dir(NGS_outputs_dir, immunarch_data_dir)
  data = repLoad(immunarch_data_dir)
  ref_samples = calc_ref_samples(data)
  chao1_estimators = calc_chao1(data)
  rarefaction_curves = calc_rarefaction_curves(data, ref_samples)
  rarefaction_plot = plot_rarefaction_curves(chao1_estimators,
                                             rarefaction_curves,
                                             ref_samples,
                                             prefix = prefix)
  return(list(ref_samples, 
              chao1_estimators,
              rarefaction_curves,
              rarefaction_plot))
}
