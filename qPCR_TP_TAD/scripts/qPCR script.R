library(tidyr)
library(ggplot2)
library(dplyr)
library(janitor)
library(finalfit)
library(qpcR)
library(readr)
library(here)

qtable<-TP_TAD_qPCR_groupe_jeudi_vendredi_donne_es_nifH_septembre_2019

write.csv2(qtable,here::here("raw","qPCR_table.csv"),row.names =FALSE)

str(qtable)

#Nettoyer
qtable<-qPCR_table
clean_names(qtable)

#Renommer
#qtable <- qtable %>% 
 # rename(code_postal=code_po_stal)


#1 Tableaux 1 : echantillons

echantillons<- qtable %>%
  filter(Content=="maïs")

echantillons_2<- qtable %>%
  filter(Content=="colza")

echantillons_3<- qtable %>%
  filter(Content=="sol nu")

echantillons<-full_join(echantillons,echantillons_2)
echantillons<-full_join(echantillons,echantillons_3)

sample<- echantillons|>select(Content, Ct)
sample<-drop_na(sample,Ct)

write.csv2(sample, here::here("data","sample.csv"),row.names = FALSE)

#2 Tableaux 2 : Standards

ST<- qtable %>%
  filter(Content=="Std 1")

STd<- qtable %>%
  filter(Content=="Std 2")

std<-full_join(ST,STd)
 
std<- std|>select(Content, Ct, `Starting Quantity (SQ)`,`Log Starting Quantity`,`SQ Mean`)

write.csv2(std, here::here("data","std.csv"),row.names = FALSE)

#3 NTC = Ctrl neg
Ctrl<- qtable %>%
  filter(Content=="NTC")

#1 Eliminer les colonnes superflues (well, biological set name, set point, well note)




head(qpcr_data)
primer_key <- data.frame(row = c("A", "B", "C", "D", "E", "F", "G", "H"),
                         primers = c(rep("HK", 4), rep("test", 4)))

tidy_data <- separate(qpcr_data, Well, into = c("row", "column"), 
                      sep = 1, convert = TRUE)
head(tidy_data)
tidy_data <- left_join(tidy_data, primer_key, by = "row")

ggplot(tidy_data, aes(x = column, y = row, fill = primers, label = Sample.Name)) +
  geom_tile(colour = "black") +
  geom_text()

#Challenge
ggplot(tidy_data, aes(x = column, y = row, fill = primers, label = Sample.Name)) +
  geom_tile(colour = "black") +
  geom_text() +
  scale_y_discrete(limits =  c("H", "G", "F", "E", "D", "C", "B", "A")) +
  scale_x_continuous(breaks = 1:12)

#Calculating average Ct values
tidy_data <- filter(tidy_data, Sample.Name != "NTC", Sample.Name != "")
head(tidy_data$Ct)
class(tidy_data$Ct)

summarised_data <- tidy_data %>%
  mutate(Ct = as.numeric(Ct)) %>%
  group_by(Sample.Name, primers) %>%
  summarise(mean_Ct = mean(Ct))

head(summarised_data)

#Challenge
summarised_data <- summarised_data %>%
  separate(Sample.Name,into = c("RNAi", "replicate"), sep = "-")

ggplot(summarised_data, aes(x = RNAi, y = mean_Ct, colour = primers)) +
  geom_point()

#DDCt methode
test_data <- summarised_data %>%
  filter(primers == "test")

ref_data <- summarised_data %>%
  filter(primers == "HK") %>%
  rename("ref_Ct" = "mean_Ct")

new_data <- left_join(test_data, ref_data, by = c( "RNAi", "replicate"))

new_data <- mutate(new_data, delta_Ct = ref_Ct - mean_Ct)
ggplot(new_data, aes(x = RNAi, y = delta_Ct)) +
  geom_point()

#Mean delta Ct for each treatment.
treatment_summary <- new_data %>%
  group_by(RNAi) %>%
  summarise(mean_delta_Ct = mean(delta_Ct))

mean_control <- filter(treatment_summary, RNAi == "Control") %>% pull(mean_delta_Ct)

# Delta delta Ct
combined_data <- new_data %>% 
  mutate(delta_delta_Ct = mean_control - delta_Ct)

ggplot(combined_data, aes(x = RNAi, y = delta_delta_Ct)) +
  geom_point()

#Calculating relative DNA concentration

combined_data <- combined_data %>%
  mutate(rel_conc = 2^-delta_delta_Ct)

ggplot(combined_data, aes(x = RNAi, y = rel_conc)) +
  geom_point() +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.3))

#Barplot with relative concentration

combined_data <- combined_data %>%
  group_by(RNAi) %>%
  mutate(mean_rel_conc = mean(rel_conc))
http://127.0.0.1:38721/graphics/db282dba-c7c9-4f32-aaa0-39a6cb5afac9.png

ggplot(combined_data, aes(x = primers.x, y = mean_rel_conc, fill = RNAi, group = RNAi)) +
  geom_col(position ="dodge") +
  geom_point(aes(y = rel_conc), position = position_dodge(width = 0.9)) +
  scale_y_continuous("Relative concentration", labels = scales::percent) +
  scale_fill_manual("", values = c("Control" = "grey", "RNAi1" = "lightblue", "RNAi2" = "steelblue4")) +
  scale_x_discrete("") +
  ggtitle("qPCR analysis of knockdown for Test") +
  theme_classic(base_size = 12)

# Enregistrer les graphiques dans un fichier pdf : format vectoriel = meilleure qualitée

pdf("qPCR.pdf")
ggplot(new_data, aes(x = RNAi, y = delta_Ct)) +
  geom_point()

ggplot(combined_data, aes(x = RNAi, y = delta_delta_Ct)) +
  geom_point()

ggplot(combined_data, aes(x = RNAi, y = rel_conc)) +
  geom_point() +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.3))

ggplot(combined_data, aes(x = primers.x, y = mean_rel_conc, fill = RNAi, group = RNAi)) +
  geom_col(position ="dodge") +
  geom_point(aes(y = rel_conc), position = position_dodge(width = 0.9)) +
  scale_y_continuous("Relative concentration", labels = scales::percent) +
  scale_fill_manual("", values = c("Control" = "grey", "RNAi1" = "lightblue", "RNAi2" = "steelblue4")) +
  scale_x_discrete("") +
  ggtitle("qPCR analysis of knockdown for Test") +
  theme_classic(base_size = 12)

dev.off() 
