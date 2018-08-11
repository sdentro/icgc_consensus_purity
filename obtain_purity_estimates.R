library(readr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(MASS)

args = commandArgs(T)
basedir = args[1]
outdir = args[2]

make_figures = F
write_intermediary_files = F

# basedir = "/Users/dentros/old_disk/Projects/icgc/consensus_subclonal_copynumber/final_run_testing/final_purity_ploidy_patch/data_bundle/"
# outdir = "~/Documents/Projects/icgc/consensusCopyNumberFinal_purity_w_sd/"

purity_file = file.path(basedir, "consensus.20170119.purity.ploidy.annotated.txt")
phylo_file = file.path(basedir, "summary_table.consensus.sample-5000.pwgs.run1.20170208.txt")
ccube_file = file.path(basedir, "summary_table_consensus_2017-02-01_run1_ccube.txt")
clonehd_file = file.path(basedir, "cloneHD.run1.summary.table.txt")
ctpsingle_file = file.path(basedir, "summary_table_run_1_CTPsingle.txt")
clip_file = file.path(basedir, "CliP_run1_Jan30_summarytable_CliP.txt")

power_file = file.path(basedir, "20170207_dpclust_run4_finalConsensusCopynum_allSegments_noReplace.txt")
aberration_summary_file = file.path(basedir, "aberration_summary.txt")
wgd_status_20170217 = file.path(basedir, "WGD-info_20170217.txt")
wgd_status_20170216 = file.path(basedir, "WGD-info_20170216.txt")

purity_cna = readr::read_tsv(purity_file)
consensus = purity_cna[,c(1,2)]
purity_cna = purity_cna[,c(-2,-3)]
phylo = readr::read_tsv(phylo_file)
ccube = readr::read_tsv(ccube_file)
clonehd = readr::read_tsv(clonehd_file)
ctpsingle = readr::read_tsv(ctpsingle_file)
clip = readr::read_tsv(clip_file)

purity_snv = data.frame(samplename=purity_cna$samplename)
purity_snv$phylowgs = phylo$purity[match(purity_snv$samplename, phylo$samplename)]
purity_snv$ccube = ccube$purity[match(purity_snv$samplename, ccube$samplename)]
purity_snv$clonehd = clonehd$purity[match(purity_snv$samplename, clonehd$samplename)]
purity_snv$ctpsingle = ctpsingle$purity[match(purity_snv$samplename, ctpsingle$samplename)]
purity_snv$clip = clip$purity[match(purity_snv$samplename, clip$samplename)]

#####################################################################
# Apply prefiltering of the calls - detect outliers
#####################################################################

#' Calculate pair-wise distances between the purity calls
#' If a single method disagrees more than the distance_threshold it is overruled
filter_purities = function(purity_table, methods_indices, distance_threshold) {
  overrulings = rep(list(list()), length(methods_indices)) 
  for (samplename in purity_table$samplename) {
    sel2 = purity_table$samplename==samplename
    sel2[is.na(sel2)] = F
    values = dist(as.numeric(purity_table[sel2, methods_indices]))
    for (i in 1:nrow(as.matrix(values))) {
      d = as.matrix(values)[i,-i]
      if (all(d > distance_threshold, na.rm=T) & !all(is.na(d))) {
        overrulings[[i]] = append(overrulings[[i]], samplename)
      }
    }
  }
  return(overrulings)
}
snv_overrulings = filter_purities(purity_snv, 2:ncol(purity_snv), 0.1)

# Do not filter out cases where multiple methods disagree with eachother - we expect some kind of consensus
r = table(unlist(lapply(snv_overrulings, unlist)))
avail_purity_counts = apply(purity_snv[,2:6], 1, function(x) { sum(!is.na(x)) })
avail_purity_counts = data.frame(samplename=purity_snv$samplename, avail_purity_counts)

#' Make sure we're not removing a critical mass. When small number of calls we don't always 
do_not_overrule = c()
for (samplename in names(r)) {
  if (r[samplename] > 2 | avail_purity_counts$avail_purity_counts[avail_purity_counts$samplename==samplename] < 3) {
    do_not_overrule = c(do_not_overrule, samplename)
  }
}
  

for (i in 1:length(snv_overrulings)) {
  overruled_samples = unlist(snv_overrulings[[i]])
  overruled_samples = overruled_samples[!overruled_samples %in% do_not_overrule]
  purity_snv[match(overruled_samples, purity_snv$samplename), i+1] = NA
}

cna_overrulings = filter_purities(purity_cna, 2:ncol(purity_cna), 0.2)
for (i in 1:length(cna_overrulings)) {
  overruled_samples = unlist(cna_overrulings[[i]])
  purity_cna[match(overruled_samples, purity_cna$samplename), i+1] = NA
}

#####################################################################
# Get densities
#####################################################################
plot_purity = function(purity_m, xlabel, consensus=NULL) {
  
  theme_plot = function(p) {
    p = p + xlim(0,1) + xlab(xlabel) + ylab("Density") +
            theme_bw() + theme(axis.title.x=element_text(colour="black",size=24,face="plain"),
                               axis.text.x=element_text(colour="black", size=22, face="plain"),
                               axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.title.y = element_text(colour="black",size=24,face="plain"),
                               strip.text.x = element_text(colour="black",size=24,face="plain"),
                               plot.title = element_text(colour="black",size=24,face="plain"),
                               legend.position="bottom")
    return(p)
  }
  
  if (all(is.na(purity_m$value))) {
    p_cna = ggplot(purity_m)
    p_cna = theme_plot(p_cna)
    
    purity_max_density = NA
    
  } else {
    
    p_cna = ggplot(purity_m, aes(x=value)) + stat_density(geom="line", size=1)
    p_cna = theme_plot(p_cna)
    
    pp_cna = ggplot_build(p_cna)
    bin_max_density = which.max(pp_cna$data[[1]]$ymax)
    max_val = pp_cna$data[[1]]$ymax[bin_max_density]
    purity_max_density = pp_cna$data[[1]]$x[bin_max_density]
    
    purity_m$y_value = max_val*0.6
    p_cna = p_cna + geom_segment(data=purity_m, mapping=aes(x=value, xend=value, y=0, yend=y_value, colour=variable))
    
    if (!is.null(consensus)) {
      consensus$y_value = max_val*1.2
      p_cna = p_cna + geom_segment(data=consensus, mapping=aes(x=value, xend=value, y=0, yend=y_value), colour="black", linetype=3)
    }
  }
  
  return(list(plot=p_cna, purity_max_density=purity_max_density))
}

# Init data
purity_cna_m = melt(purity_cna)
purity_snv_m = melt(purity_snv)
purity_consensus_m = melt(consensus)

purity_cna_m$y_value = 0.2
purity_snv_m$y_value = 0.2
purity_consensus_m$my_value = 0.2

# samplename = "00508f2b-36bf-44fc-b66b-97e1f3e40bfa"
# samplename = "0176cf1d-0760-4769-a493-277f4bb7585e"

all_consensus_purities = data.frame()

for (samplename in purity_cna$samplename) {
  outfile = file.path(outdir, paste0(samplename, "_purity.png"))
  # if (!file.exists(outfile)) {
  if (! samplename %in% all_consensus_purities$samplename) {
    cna_temp = purity_cna_m[purity_cna_m$samplename==samplename,]
    snv_temp = purity_snv_m[purity_snv_m$samplename==samplename,]
    all_temp = rbind(cna_temp, snv_temp)
    cons_temp = purity_consensus_m[purity_consensus_m$samplename==samplename,]
    
    res = plot_purity(cna_temp, "CNA purity", consensus=cons_temp)
    p_cna = res$plot
    cna_purity_max_density = res$purity_max_density
    
    res = plot_purity(snv_temp, "SNV purity")
    snv_purity_max_density = res$purity_max_density
    
    res = plot_purity(all_temp, "All purity")
    p_all = res$plot
    all_purity_max_density = res$purity_max_density
    
    data_table = data.frame(type=c("Current cons", "CNA peak", "SNV peak", "All peak", "Purity SD"), 
               purity=c(round(purity_consensus_m[purity_consensus_m$samplename==samplename,]$value, 3), 
                        round(cna_purity_max_density, 3), 
                        round(snv_purity_max_density, 3),
                        round(all_purity_max_density, 3),
                        round(sd(all_temp$value, na.rm=T), 3)))
    
    cna_closest_method_index = which.min(abs(cna_temp$value - data_table$purity[2]))
    snv_closest_method_index = which.min(abs(snv_temp$value - data_table$purity[3]))
    all_closest_method_index = which.min(abs(all_temp$value - data_table$purity[4]))
    
    
    data_table_2 = data.frame(type=c("CNA", "SNV", "All"), 
                              method=c(ifelse(length(cna_closest_method_index) > 0, as.character(cna_temp$variable[cna_closest_method_index]), "NA"), 
                                       ifelse(length(snv_closest_method_index) > 0, as.character(snv_temp$variable[snv_closest_method_index]), "NA"), 
                                       as.character(all_temp$variable[all_closest_method_index])), 
                              purity=c(ifelse(length(cna_closest_method_index) > 0, round(cna_temp$value[cna_closest_method_index], 3), NA), 
                                       ifelse(length(snv_closest_method_index) > 0, round(snv_temp$value[snv_closest_method_index], 3), NA), 
                                       round(all_temp$value[all_closest_method_index], 3)))
    
    if (length(cna_closest_method_index)==0) {
      data_table_2$method[1] = NA
      data_table_2$purity[1] = NA
    }
    
    if (make_figures) {
      png(outfile, height=350, width=1200)
      grid.arrange(arrangeGrob(p_cna, p_all, 
                               arrangeGrob(tableGrob(data_table, rows=NULL), tableGrob(data_table_2, rows=NULL), ncol=1), 
                               ncol=3,
                               widths=c(2,2,1)))
      dev.off()
    }
    
    all_consensus_purities = rbind(all_consensus_purities,
                                   data.frame(samplename=samplename,
                                              current_cons=data_table$purity[1],
                                              cna_peak=data_table$purity[2],
                                              snv_peak=data_table$purity[3],
                                              all_peak=data_table$purity[4],
                                              cna_closest_method=data_table_2$purity[1],
                                              snv_closest_method=data_table_2$purity[2],
                                              all_closest_method=data_table_2$purity[3]))
  }
}
if (write_intermediary_files)
  write.table(all_consensus_purities, file=file.path(basedir, "combined_consensus_purities.txt"), row.names=F, sep="\t", quote=F)

#####################################################################
# Plot complete distribution for all methods
#####################################################################
if (make_figures) {
  all_consensus_purities = read.table(file.path(basedir, "combined_consensus_purities.txt"), header=T, stringsAsFactors=F)
  m = match(all_consensus_purities$samplename, purity_snv$samplename)
  all_consensus_purities = cbind(all_consensus_purities, purity_snv[m,2:7])
  
  m = match(all_consensus_purities$samplename, purity_cna$samplename)
  all_consensus_purities = cbind(all_consensus_purities, purity_cna[m,2:7])
  
  colnames(all_consensus_purities)[9] = "clonehd_snv"
  ymax = 100
  all = ggplot() + geom_histogram(data=all_consensus_purities, mapping=aes(x=current_cons, y=..count..), binwidth=0.01, colour="black", fill="grey") + xlim(0, 1) + ylim(0, ymax)
  p1 = ggplot() + geom_histogram(data=all_consensus_purities, mapping=aes(x=phylowgs, y=..count..), binwidth=0.01, colour="black", fill="grey") + xlim(0, 1) + ylim(0, ymax)
  p2 = ggplot() + geom_histogram(data=all_consensus_purities, mapping=aes(x=ccube, y=..count..), binwidth=0.01, colour="black", fill="grey") + xlim(0, 1) + ylim(0, ymax)
  p3 = ggplot() + geom_histogram(data=all_consensus_purities, mapping=aes(x=clonehd_snv, y=..count..), binwidth=0.01, colour="black", fill="grey") + xlim(0, 1) + ylim(0, ymax)
  p4 = ggplot() + geom_histogram(data=all_consensus_purities, mapping=aes(x=ctpsingle, y=..count..), binwidth=0.01, colour="black", fill="grey") + xlim(0, 1) + ylim(0, ymax)
  p5 = ggplot() + geom_histogram(data=all_consensus_purities, mapping=aes(x=clip, y=..count..), binwidth=0.01, colour="black", fill="grey") + xlim(0, 1) + ylim(0, ymax) 
  png(file.path(basedir, "purity_histogram_snv_methods.png"), width=1000, height=600)
  grid.arrange(grid.arrange(all, p1, p2, ncol=1),
               grid.arrange(p3, p4, p5, ncol=1),
               ncol=2)
  dev.off()
  
  p1 = ggplot() + geom_histogram(data=all_consensus_purities, mapping=aes(x=aceseq, y=..count..), binwidth=0.01, colour="black", fill="grey") + xlim(0, 1) + ylim(0, ymax)
  p2 = ggplot() + geom_histogram(data=all_consensus_purities, mapping=aes(x=absolute, y=..count..), binwidth=0.01, colour="black", fill="grey") + xlim(0, 1) + ylim(0, ymax)
  p3 = ggplot() + geom_histogram(data=all_consensus_purities, mapping=aes(x=clonehd, y=..count..), binwidth=0.01, colour="black", fill="grey") + xlim(0, 1) + ylim(0, ymax)
  p4 = ggplot() + geom_histogram(data=all_consensus_purities, mapping=aes(x=sclust, y=..count..), binwidth=0.01, colour="black", fill="grey") + xlim(0, 1) + ylim(0, ymax)
  p5 = ggplot() + geom_histogram(data=all_consensus_purities, mapping=aes(x=jabba, y=..count..), binwidth=0.01, colour="black", fill="grey") + xlim(0, 1) + ylim(0, ymax) 
  p6 = ggplot() + geom_histogram(data=all_consensus_purities, mapping=aes(x=battenberg, y=..count..), binwidth=0.01, colour="black", fill="grey") + xlim(0, 1) + ylim(0, ymax)
  png(file.path(basedir, "purity_histogram_cna_methods.png"), width=1000, height=600)
  grid.arrange(grid.arrange(p1, p2, p3, ncol=1),
               grid.arrange(p4, p5, p6, ncol=1),
               ncol=2)
  dev.off()
}

#####################################################################
# Learn a standard deviation using the given purities
#####################################################################
all_consensus_purities = read.table(file.path(basedir, "combined_consensus_purities.txt"), header=T, stringsAsFactors=F)
m = match(all_consensus_purities$samplename, purity_snv$samplename)
all_consensus_purities = cbind(all_consensus_purities, purity_snv[m,2:5])

m = match(all_consensus_purities$samplename, purity_cna$samplename)
all_consensus_purities = cbind(all_consensus_purities, purity_cna[m,2:7])

all_consensus_purities = all_consensus_purities[, !grepl("bayclone", colnames(all_consensus_purities))]

fit_dist_1 = function(all_purities, cons_purity, fitmethod="BFGS") {
  x = as.numeric(all_purities)
  x = x[!is.na(x)]
  sd_est = fitdistr(x=x, 
                    densfun="normal", 
                    mean=cons_purity, 
                    method=fitmethod)
  return(sd_est$estimate["sd"])
}

fit_dist_2 = function(all_purities, cons_purity, start, fitmethod="BFGS") {
  x = as.numeric(all_purities)
  x = x[!is.na(x)]
  sd_est = fitdistr(x=x, 
                    dnorm, 
                    mean=cons_purity, 
                    start=start, 
                    method=fitmethod)
  return(sd_est$estimate)
}

calc_ll = function(all_purities, cons_purity, sd_purity) {
  x = as.numeric(all_purities)
  x = x[!is.na(x)]
  return(sum(dnorm(x=x, mean=cons_purity, sd=sd_purity, log=T)))
}

# Fit in two ways to see what is more reliable
all_consensus_purities$sd_1 = NA
all_consensus_purities$sd_2 = NA
for (i in 1:nrow(all_consensus_purities)) {
  dat = all_consensus_purities[i,7:16]
  all_consensus_purities$sd_1[i] = fit_dist_1(dat, all_consensus_purities$all_peak[i])

  raw_sd = sd(dat, na.rm=T)
  # Check whether the data does not agree too much for fitting to fail
  if (raw_sd > 0.01 & !is.na(raw_sd)) {
    all_consensus_purities$sd_2[i] = fit_dist_2(all_purities=dat, 
                                                cons_purity=all_consensus_purities$all_peak[i],
                                                start=list(sd=raw_sd))
  } else {
    all_consensus_purities$sd_2[i] = 0
  }
}

# Calculate the mad for each sample
all_consensus_purities$mad = apply(all_consensus_purities[, 7:17], 1, mad, na.rm=T)
if (write_intermediary_files)
  write.table(all_consensus_purities, file=file.path(basedir, "combined_consensus_purities_w_sd.txt"), row.names=F, sep="\t", quote=F)

if (make_figures) {
  all_consensus_purities$peak_diff = abs(all_consensus_purities$cna_peak-all_consensus_purities$snv_peak)
  p = ggplot(all_consensus_purities) + aes(x=sd_2, y=peak_diff) + geom_point() + xlim(0,1) + ylim(0,1)
  png(file.path(basedir, "peakDiff_vs_sd.png"), height=500, width=500)
  print(p)
  dev.off()
  
  p = ggplot(all_consensus_purities) + aes(x=mad, y=peak_diff) + geom_point() + xlim(0,1) + ylim(0,1)
  png(file.path(basedir, "peakDiff_vs_mad.png"), height=500, width=500)
  print(p)
  dev.off()
}

# Get list of samples with high 1+1:
# cn_states = read.table(cn_states_file, header=T, stringsAsFactors=F)
# m = match(all_consensus_purities$samplename, cn_states$samplename)
# all_consensus_purities$size_1_1 = cn_states$tot_1_1[m]

# Annotate power per sample
power = read.table(power_file, header=T, stringsAsFactors=F)
power$total_snv = power$num_clonal+power$num_subclonal
m = match(all_consensus_purities$samplename, power$samplename)
all_consensus_purities$power = power$nrpcc[m]
all_consensus_purities$num_snv = power$total_snv[m]

aberration_summary = read.table(aberration_summary_file, header=T, stringsAsFactors=F)
m = match(all_consensus_purities$samplename, aberration_summary$samplename)
all_consensus_purities$frac_aberrated = aberration_summary$frac_aberrated[m]

all_consensus_purities$power_class = NA
all_consensus_purities$power_class[all_consensus_purities$power < 5] = "low"
all_consensus_purities$power_class[all_consensus_purities$power > 5 & all_consensus_purities$power < 10] = "medium"
all_consensus_purities$power_class[all_consensus_purities$power > 10] = "high"
all_consensus_purities$power_class = factor(all_consensus_purities$power_class)
# p = ggplot(all_consensus_purities) + aes(x=frac_aberrated, y=num_snv, colour=power_class) + geom_point() + scale_y_log10()
# print(p)


# # Select the consensus
# quiet = all_consensus_purities[all_consensus_purities$size_1_1 > 2700,]
# quiet[quiet$power_class=="low",]
# 
# quiet = all_consensus_purities[all_consensus_purities$frac_aberrated < 0.05,]
# head(quiet[with(quiet, order(frac_aberrated, decreasing=T)),], 25)
# 
# quiet[quiet$power_class=="low",]


all_consensus_purities$aberrated_category = NA
all_consensus_purities$aberrated_category[all_consensus_purities$frac_aberrated > 0.1] = "> 0.1"
all_consensus_purities$aberrated_category[all_consensus_purities$frac_aberrated > 0.09 & all_consensus_purities$frac_aberrated < 0.1] = "> 0.09, < 0.1"
all_consensus_purities$aberrated_category[all_consensus_purities$frac_aberrated > 0.08 & all_consensus_purities$frac_aberrated < 0.09] = "> 0.08, < 0.09"
all_consensus_purities$aberrated_category[all_consensus_purities$frac_aberrated > 0.07 & all_consensus_purities$frac_aberrated < 0.08] = "> 0.07, < 0.08"
all_consensus_purities$aberrated_category[all_consensus_purities$frac_aberrated > 0.06 & all_consensus_purities$frac_aberrated < 0.07] = "> 0.06, < 0.07"
all_consensus_purities$aberrated_category[all_consensus_purities$frac_aberrated > 0.05 & all_consensus_purities$frac_aberrated < 0.06] = "> 0.05, < 0.06"
all_consensus_purities$aberrated_category[all_consensus_purities$frac_aberrated > 0.04 & all_consensus_purities$frac_aberrated < 0.05] = "> 0.04, < 0.05"
all_consensus_purities$aberrated_category[all_consensus_purities$frac_aberrated > 0.03 & all_consensus_purities$frac_aberrated < 0.04] = "> 0.03, < 0.04"
all_consensus_purities$aberrated_category[all_consensus_purities$frac_aberrated > 0.02 & all_consensus_purities$frac_aberrated < 0.03] = "> 0.02, < 0.03"
all_consensus_purities$aberrated_category[all_consensus_purities$frac_aberrated > 0.01 & all_consensus_purities$frac_aberrated < 0.02] = "> 0.01, < 0.02"
all_consensus_purities$aberrated_category[all_consensus_purities$frac_aberrated < 0.01] = "< 0.01"
all_consensus_purities$aberrated_category = factor(all_consensus_purities$aberrated_category, 
                                                   levels=c("> 0.1", "> 0.09, < 0.1", "> 0.08, < 0.09", "> 0.07, < 0.08", "> 0.06, < 0.07", "> 0.05, < 0.06", "> 0.04, < 0.05", "> 0.03, < 0.04", "> 0.02, < 0.03", "> 0.01, < 0.02", "< 0.01"))
all_consensus_purities$peak_diff = abs(all_consensus_purities$snv_peak - all_consensus_purities$all_peak)

if (make_figures) {
  p = ggplot(all_consensus_purities[, c("peak_diff", "aberrated_category")]) + aes(x=aberrated_category, y=peak_diff) + geom_boxplot()
  png(file.path(basedir, "peak_diff_vs_aberrated.png"), height=500, width=900)
  print(p)
  dev.off()
  
  p = ggplot(all_consensus_purities) + aes(x=frac_aberrated, y=num_snv, colour=power_class) + geom_point() + scale_y_log10() + geom_vline(xintercept = 0.08)
  png(file.path(basedir, "frac_abberated_vs_num_snv.png"), height=500, width=900)
  print(p)
  dev.off()
}

# temp = all_consensus_purities[all_consensus_purities$aberrated_category=="> 0.07, < 0.08",]
# head(temp[with(temp, order(peak_diff, decreasing=T)),], 25)
# 
# quiet_low = all_consensus_purities[all_consensus_purities$frac_aberrated < 0.08 & all_consensus_purities$power_class=="low",]


#####################################################################
# Choose the purity - using the threshold
#####################################################################
all_consensus_purities$new_consensus = all_consensus_purities$all_closest_method
all_consensus_purities$consensus_choice = "all_peak"
all_consensus_purities$new_consensus[all_consensus_purities$frac_aberrated < 0.08] = all_consensus_purities$snv_closest_method[all_consensus_purities$frac_aberrated < 0.08]
all_consensus_purities$consensus_choice[all_consensus_purities$frac_aberrated < 0.08] = "snv_peak"
if (write_intermediary_files)
  write.table(all_consensus_purities, file=file.path(basedir, "combined_consensus_purities_w_sd_allAnno.txt"), row.names=F, sep="\t", quote=F)

#####################################################################
# Create final formatting
#####################################################################

wgd_status_16 = read.table(wgd_status_20170216, header=T, stringsAsFactors=F, sep="\t")
wgd_status_17 = read.table(wgd_status_20170217, header=T, stringsAsFactors=F, sep="\t")


wgd_status2_sub = wgd_status_17[row.names(wgd_status_17) %in% row.names(wgd_status_16),]
wgd_status2_sub = rbind(wgd_status2_sub, wgd_status_16[!row.names(wgd_status_16) %in% row.names(wgd_status_17),])
wgd_status2_sub$samplename = row.names(wgd_status2_sub)
wgd_status2_sub = wgd_status2_sub[with(wgd_status2_sub, order(samplename)),]


consensus_output = all_consensus_purities[, c("samplename", "new_consensus", "mad")]

consensus_output = data.frame(samplename=all_consensus_purities[, c("samplename")],
                              purity=all_consensus_purities[, c("new_consensus")],
                              ploidy=power$ploidy,
                              purity_conf_mad=round(all_consensus_purities[, c("mad")], 3), stringsAsFactors=F)
consensus_output$wgd_status = "no_wgd"

# present=red, possible=blue, absent=black, blue ones non-WGD
consensus_output$wgd_status[wgd_status2_sub$WGD_call=="present"] = "wgd"
if (write_intermediary_files)
  write.table(consensus_output, file=file.path(basedir, "combined_consensus_purities_w_sd_allAnno_finalSlimmed.txt"), row.names=F, sep="\t", quote=F)

current_date = gsub("-", "", Sys.Date())
write.table(consensus_output, file=file.path(basedir, paste0("consensus.", current_date, ".purity.ploidy.txt")), row.names=F, sep="\t", quote=F)
