##uses the tama_observations_bed.py file to produce the csv##

data <- read.csv("Gac_white_70hpf_ens_as_1_distributions.csv")

hist(data$size_dif)
hist(data$start_diff)
hist(data$end_diff)

mean(data$size_dif)
mean(data$start_diff)
mean(data$end_diff)


max(data$size_dif)
max(data$start_diff)
max(data$end_diff)

median(data$size_dif)
median(data$start_diff)
median(data$end_diff)

library(ggplot2)
library(RColorBrewer)

ggplot(data, aes(x=size_dif)) + geom_histogram(fill="aquamarine4", bins=50) + xlim(c(-50000, 50000)) + scale_y_continuous(name="Counts", limits=c(0, 50000)) +
  theme(panel.background = element_blank(), panel.grid.major.y = element_line( size=.2, color="grey"),
        axis.text.x=element_text(angle=45, hjust = .9, size=10), axis.text.y=element_text(size=10), axis.title.y = element_text(size=12), axis.title.x = element_text(size=12), axis.ticks = element_blank(), legend.position="none") +
  xlab("Size difference (Iso-Seq minus Ensembl)") +
  geom_vline(xintercept=median(data$end_diff), lwd=.5, linetype=2, color="black")

ggplot(data, aes(x=start_diff)) + geom_histogram(fill="mediumorchid4", bins=50) + xlim(c(-50000, 50000)) + scale_y_continuous(name="Counts", limits=c(0, 50000)) +
  theme(panel.background = element_blank(), panel.grid.major.y = element_line( size=.2, color="grey"),
        axis.text.x=element_text(angle=45, hjust = .9, size=10), axis.text.y=element_text(size=10), axis.title.y = element_text(size=12), axis.title.x = element_text(size=12), axis.ticks = element_blank(), legend.position="none") +
  xlab("Starting position difference (Iso-Seq minus Ensembl)") +
  geom_vline(xintercept=median(data$end_diff), lwd=.5, linetype=2, color="black")

ggplot(data, aes(x=end_diff)) + geom_histogram(fill="hotpink1", bins=50) + xlim(c(-50000, 50000)) + scale_y_continuous(name="Counts", limits=c(0, 50000))  +
  theme(panel.background = element_blank(), panel.grid.major.y = element_line( size=.2, color="grey"),
        axis.text.x=element_text(angle=45, hjust = .9, size=10), axis.text.y=element_text(size=10), axis.title.y = element_text(size=12), axis.title.x = element_text(size=12), axis.ticks = element_blank(), legend.position="none") +
  xlab("Ending position difference (Iso-Seq minus Ensembl)") +
  geom_vline(xintercept=median(data$end_diff), lwd=.5, linetype=2, color="black") 


