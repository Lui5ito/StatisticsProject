library(Seurat)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(plotly)
library(cowplot)

rm(list=ls())

#On récupère les données fournies
#data <- Read10X_h5("raw_feature_bc_matrix.D7.h5")

#Les données sont sous formes d'une "liste" de matrices. On les extraits. 
#gene_expression <- data$`Gene Expression`
#antibody_capture <- data$`Antibody Capture`
#save(gene_expression, file = "gene_expression.RData")
#save(antibody_capture, file = "antibody_capture.RData")

load("gene_expression.RData")
load("antibody_capture.RData")

#Mais ce sont encore des objet particulier. On ne peut pas les trasnformer en dataframe car les matrices sont trop grosses.
#Ici ça explique comment gérer un objetdgCMatrix : https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/

as.data.frame(gene_expression[c(1:3), c(1:3)])
summary(gene_expression)
str(gene_expression)

#Liste de toutes les non-zéros valeures.
gene_expression@x
summary(gene_expression@x)
str(gene_expression@x)

#Renvoie le numéro de la ligne de chaque valeure non-zéro. (Correspondance avec x)
gene_expression@i
summary(gene_expression@i)

#Liste de taille le nombre de colonne, et renvoie le nombre de non-zéro par colonne (en faisant la différence, cf le site)
gene_expression@p
summary(gene_expression@p)

#Nombre de lignes et nombre de colonnes
gene_expression@Dim

#Donne le nom des lignes s'ils existent et le nom des colonnes s'ils existent
gene_expression@Dimnames

# A quoi ressemble les lignes de genes_expression
as.data.frame(gene_expression@Dimnames[[1]]) %>% slice(1:50) %>% View()
nrow_genes <- gene_expression@Dimnames[[1]] %>% length() %>% as.numeric()
# A quoi ressemble les colonnes de genes_expression
as.data.frame(gene_expression@Dimnames[[2]]) %>% slice(1:50) %>% View()
ncol_genes <- gene_expression@Dimnames[[2]] %>% length() %>% as.numeric()


#Pour l'instant la matrice qui nous intéresse c'est gene_expression. C'est la qu'on trouve l'information des goutelettes.


#Ici on veut calculer le nombre de transcrit par goutelette.
sum_gene_df <- as.data.frame(colSums(gene_expression)) %>% 
  rename("nb_transcrit" =`colSums(gene_expression)`) %>% 
  filter(nb_transcrit != 0)


#plot(sum_gene_df$nb_transcrit)

#plot(sum_gene_df$nb_transcrit, log='y')

#plot(sum_gene_df$nb_transcrit, pch='.', log='y')
#abline(h = 2000, col = "red")
#abline(h = 12000, col = "red")


ggplot(sum_gene_df, aes(x = seq_along(nb_transcrit), y = nb_transcrit)) +
  geom_point(shape = '.') +
  scale_y_continuous(trans = 'log10') +
  geom_hline(yintercept = 2000, linetype = "dashed", color = 'red', linewidth = 0.6) +
  geom_hline(yintercept = 12000, linetype = "dashed", color = 'red', linewidth = 0.6) +
  theme_bw() +
  ggtitle("Scatter plot,\n and a first determination of the groups") +
  xlab("Index") + ylab("Number of transcript")

sum_gene_df$classe <- case_when(sum_gene_df$nb_transcrit < 2000 ~ "1",
                                sum_gene_df$nb_transcrit >= 2000 & sum_gene_df$nb_transcrit < 12000 ~ "2",
                                sum_gene_df$nb_transcrit >= 12000 ~ "3")

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

sum_gene_df %>% 
  filter(classe == i) %>% 
  summary()

vp1 <- ggplot(sum_gene_df %>% 
         filter(classe == 1) %>% 
           filter(nb_transcrit > 30), aes(x = as.factor(classe), y = nb_transcrit)) +
  geom_violin(trim = FALSE, fill = "lightblue", color = "#66CCFF") +
  stat_summary(fun.data=data_summary, color="red") +
  theme_bw() + 
  ggtitle("Violin plot for group 1") +
  xlab("Group 1") + ylab("Number of transcript")
bp1 <- ggplot(sum_gene_df %>% 
         filter(classe == 1) %>% 
           filter(nb_transcrit > 30), aes(x = nb_transcrit, y = after_stat(density))) +
         geom_histogram(aes(x = nb_transcrit), alpha = 0.8, color = "#66CCFF", fill = "lightblue") +
  theme_bw() +
  geom_density(aes(x = nb_transcrit), color = "red", linewidth = 0.666) + 
  ggtitle("Histogramm plot and \nnp-estimated density for group 1") +
  xlab("Number of transcript") + ylab("Number of goutellette")

data_cdf_poisson = rpois(1000, 1)
cdf1 <- ggplot(data.frame(data_cdf_poisson), aes(x = data_cdf_poisson)) + 
  stat_ecdf(geom = "point",  color = "red") +
  labs(title="Cumulative Density Function \nfor Poisson",
       y = "F(x)", x="x") +
  theme_bw()

ecdf1 <- ggplot(sum_gene_df %>% 
         filter(classe == 1)%>% 
         filter(nb_transcrit > 30), aes(nb_transcrit)) + 
  stat_ecdf(geom = "point",  color = "#66CCFF") +
  labs(title="Empirical Cumulative \nDensity Function",
       y = "Nombre de transcrit", x="Index de la goutelette") +
  theme_bw()

gp1 <- plot_grid(vp1, bp1, cdf1, ecdf1, ncol = 2, nrow = 2)



vp2 <- ggplot(sum_gene_df %>% 
                filter(classe == 2), aes(x = as.factor(classe), y = nb_transcrit)) +
  geom_violin(trim = FALSE, fill = "lightblue", color = "#66CCFF") +
  stat_summary(fun.data=data_summary, color="red") +
  theme_bw() + 
  ggtitle("Violin plot for group 2") +
  xlab("Group 2") + ylab("Number of transcript")
bp2 <- ggplot(sum_gene_df %>% 
                filter(classe == 2), aes(x = nb_transcrit, y = after_stat(density))) +
  geom_histogram(aes(x = nb_transcrit), alpha = 0.8, color = "#66CCFF", fill = "lightblue") +
  theme_bw() +
  geom_density(aes(x = nb_transcrit), color = "red", linewidth = 0.666) + 
  ggtitle("Histogramm plot and \nnp-estimated density for group 2") +
  xlab("Number of transcript") + ylab("Number of goutellette")

data_cdf_ln1 = rlnorm(1000, 0, 1)
cdf2 <- ggplot(data.frame(data_cdf_ln1), aes(x = data_cdf_ln1)) + 
  stat_ecdf(geom = "point",  color = "red") +
  labs(title="Cumulative Density Function \nfor log-Normal",
       y = "F(x)", x="x") +
  theme_bw()

ecdf2 <- ggplot(sum_gene_df %>% 
                  filter(classe == 2), aes(nb_transcrit)) + 
  stat_ecdf(geom = "point",  color = "#66CCFF") +
  labs(title="Empirical Cumulative \nDensity Function",
       y = "Nombre de transcrit", x="Index de la goutelette") +
  theme_bw()
gp2 <- plot_grid(vp2, bp2, cdf2, ecdf2, ncol = 2, nrow = 2)



#Group 3
vp3 <- ggplot(sum_gene_df %>% 
                filter(classe == 3), aes(x = as.factor(classe), y = nb_transcrit)) +
  geom_violin(trim = FALSE, fill = "lightblue", color = "#66CCFF") +
  stat_summary(fun.data=data_summary, color="red") +
  theme_bw() + 
  ggtitle("Violin plot for group 3") +
  xlab("Group 3") + ylab("Number of transcript")
bp3 <- ggplot(sum_gene_df %>% 
                filter(classe == 3), aes(x = nb_transcrit, y = after_stat(density))) +
  geom_histogram(aes(x = nb_transcrit), alpha = 0.8, color = "#66CCFF", fill = "lightblue") +
  theme_bw() +
  geom_density(aes(x = nb_transcrit), color = "red", linewidth = 0.666) + 
  ggtitle("Histogramm plot and \nnp-estimated density for group 3") +
  xlab("Number of transcript") + ylab("Number of goutellette")

data_cdf_ln2 = rlnorm(1000, 100, 1)
cdf3 <- ggplot(data.frame(data_cdf_ln2), aes(x = data_cdf_ln2)) + 
  stat_ecdf(geom = "point",  color = "red") +
  labs(title="Cumulative Density Function \nfor log-Normal",
       y = "F(x)", x="x") +
  theme_bw()

ecdf3 <- ggplot(sum_gene_df %>% 
                  filter(classe == 3), aes(nb_transcrit)) + 
  stat_ecdf(geom = "point",  color = "#66CCFF") +
  labs(title="Empirical Cumulative \nDensity Function",
       y = "Nombre de transcrit", x="Index de la goutelette") +
  theme_bw()
gp3 <- plot_grid(vp3, bp3, cdf3, ecdf3, ncol = 2, nrow = 2)


fp <- plot_grid(gp1, gp2, gp3, labels = c("Group 1", "Group 2", "Group 3"), label_x = c(0.45, 0.45, 0.45), label_y = c(0.15, 0.15, 0.15)  ,ncol = 1, nrow = 3)


gp1
gp2
gp3








library(KernSmooth)
plot(bkde((sum_gene_df %>% 
             filter(classe == 1) %>% 
             filter(nb_transcrit > 25))$nb_transcrit))
##Adding the gamma distribution
#define x-values
x <- seq(0, 100, by=0.1)   
#calculate gamma density for each x-value
y <- dgamma(x, shape = 10, rate = 0.5) 
plot(y)
#create cumulative density plot

lines(y)

?density

#logL(lambda) = sum(i=1:n) (-lambda+sum_gene_df$nb_transcrit[i]*log(lambda) - log(sum_gene_df$nb_transcrit[i]!))

logL <- function(lambda, x) {
  sum(-lambda+x*log(lambda) - lfactorial(x))
}
enattendant <- sum_gene_df%>% 
  filter(classe == 1) %>% 
  filter(nb_transcrit > 25)

estimation <- optim(par = 1, logL, x = enattendant$nb_transcrit, method = "Brent", lower = 0, upper = 10)
lambda <- estimation$par
lambda

x_seq <- seq(0, max(enattendant$nb_transcrit), length.out = 100)
y_pois <- dpois(x_seq, lambda)
x <- rpois(1000000, lambda)
hist(x)
lines(x_seq, y_pois, col = "red", lwd = 2)


sum_gene_df %>% 
  filter(classe == 2) %>% 
  summary()

ggplot(sum_gene_df %>% 
                        filter(classe == 2), aes(x = nb_transcrit, y = after_stat(density))) +
  geom_histogram(aes(x = nb_transcrit), alpha = 0.8, color = "#66CCFF", fill = "lightblue") +
  theme_bw() +
  geom_density(aes(x = nb_transcrit), color = "red", linewidth = 0.666)

sum_gene_df %>% 
  filter(classe == 3) %>% 
  summary()

ggplot(sum_gene_df %>% 
                        filter(classe == 3), aes(x = nb_transcrit, y = after_stat(density))) +
  geom_histogram(aes(x = nb_transcrit), alpha = 0.8, color = "#66CCFF", fill = "lightblue") +
  theme_bw() +
  geom_density(aes(x = nb_transcrit), color = "red", linewidth = 0.666)



str(sum_gene_df)

sum_gene_df_sans1 <- subset(sum_gene_df, nb_transcrit != 1)
sum_gene_df_sans12 <- subset(sum_gene_df_sans1, nb_transcrit != 2)
sum_gene_df_sans123 <- subset(sum_gene_df_sans12, nb_transcrit != 3)
plot(sum_gene_df_sans123$nb_transcrit, log='y')

str(sum_gene_df)
summary(sum_gene_df)

ggplot(sum_gene_df, aes(x = nb_transcrit, y = c(1:327395))) +
  geom_point() +
  scale_y_continuous(trans='log10') +
  geom_jitter()
  

nb_transcrit <- as.data.frame(gene_expression@x)

group_nb_transcrit <- as.data.frame(table(as.character(nb_transcrit$`gene_expression@x`)))
group_nb_transcrit$Var1 <- as.numeric(group_nb_transcrit$Var1)
group_nb_transcrit <- group_nb_transcrit[order(group_nb_transcrit$Var1), ]
group_nb_transcrit <- group_nb_transcrit %>% slice(-1)
plot(group_nb_transcrit)


plot(nb_transcrit)

group_gene_df <- table(nb_transcrit)

plot(gene_expression@x)

summary(sum_gene_df)
library(pastecs)

stat.desc(sum_gene_df$nb_transcrit, basic=TRUE, desc=TRUE, norm=FALSE, p=0.95)


#On veut récuperer les fréquences des goutelettes qui ont x transcrits
freq_table <- as.data.frame(ftable(as.data.frame(colSums(gene_expression))))
colnames(freq_table)[colnames(freq_table) == 'colSums.gene_expression.'] <- 'Nombre_de_goutelettes_qui_ont_le_même_nombre_de_transcrits'
colnames(freq_table)[colnames(freq_table) == 'Freq'] <- 'Nombre_de_transcrits'

#On plot sur une échelle logarithmique le nombre de goutelettes qui ont le même nombre de transcrit.
#Logarithmique car le nombre de transcrit commence à 0 et peut aller à plusieurs milliers.
ggplot(freq_table, aes(x = Nombre_de_transcrits, y = as.numeric(Nombre_de_goutelettes_qui_ont_le_même_nombre_de_transcrits))) + 
  geom_point() +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  geom_jitter() +
  labs(title = 'Nombre de transcrits par goutelettes en fonctions des goutelettes qui ont le même nombre de transcrits',subtitle = 'Echelle logarithmique en abscisse', x = 'Nombre de goutelettes qui ont le même nombre de transcrit.', y = 'Nombre de transcrits dans ces goutelettes')

sum_gene_df$classes <-case_when(log(sum_gene_df$nb_transcrit) > 1e+04 ~ 1,
                                log(sum_gene_df$nb_transcrit) > 1e+01 & log(sum_gene_df$nb_transcrit) < 1e+04 ~ 2,
                                log(sum_gene_df$nb_transcrit) <= 1e+01,
                                .default = 0)

summary(sum_gene_df$classes)
