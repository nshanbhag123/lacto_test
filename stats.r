library(tidyverse)
library(dplyr)
library(ggplot2)
library(readr)
library(EnvStats)
library(argparse)

parser <- ArgumentParser(description='Visualize ANI Distribution and Compute Empirical P-Values')
parser$add_argument('--simu_anis', help='specify path to collate ani file')
parser$add_argument('--test_anis', help = 'specify path to test ani file')
parser$add_argument('--outpath', help = 'specify path to output folder')
args <- parser$parse_args()


# wk.dir <- setwd("/home/nshanbhag/lacto_compare/")
# ani <- paste(wk.dir, "/", "out", "/", "coallate_ani.csv", sep = "")
df <- read.csv(args$simu_anis)
df <- pivot_longer(df, cols = -X, names_to = "strain", values_to = "ani")

###############################
stats_output <- paste(args$outpath, "/", "stats.txt", sep = "")
sink(stats_output)

ani <- df$ani
strain <- df$strain

model <- lm(ani ~ strain)

shapiro <- shapiro.test(model$residuals)
bartlett <- bartlett.test(model$residuals, strain)
anova <- anova(model)
tukey <- TukeyHSD(aov(ani~strain), conf.level=0.95)

if (shapiro$p.value & bartlett$p.value < 0.05) {
    print("Normality and Homoscedasticity Conditions were Met - ANOVA test printed below")
    print(anova)
} else {
    print("Normality and Homoscedasticity Conditions were not Met - non-parametric Kruskal test printed below")
    kruskal <- kruskal.test(ani~strain)
    print(kruskal)
}

if (anova$"Pr(>F)"[1] < 0.05){
    print("Because the ANOVA test was significant, perform post-hoc test")
    print("Tukey's Test to Determine the Pairwise Difference of Means Between Groups")
    print(tukey)
}

sink()
##################################
boot_strap <-sample(ani,size=10000,replace = TRUE) 
boot_strap_df <- data.frame(boot_strap)
gg <- boot_strap_df %>% ggplot( aes(x=boot_strap)) + geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +ggtitle("Density plot of ANI values") + xlab("ANI values(bootstrapped)")
gg_output <- paste(args$outpath, "/", "density.png", sep = "") 
ggsave(gg_output, plot = gg, device = "png")


real_ani <- read_tsv(args$test_anis, col_names = FALSE)
col_names <- c("strain1", "strain2", "ani", "seq_fragments", "aligned")
names(real_ani) <- col_names

emp<-rep(0,nrow(real_ani))
for (i in 1:nrow(real_ani)){
    {emp[i]<-pemp(real_ani$ani[i],boot_strap,discrete = FALSE)}
}
output <- cbind(real_ani, emp)
outpath <- paste(outpath, "/", "empirical_ps.csv" ,sep = "")

write.csv(output, outpath)