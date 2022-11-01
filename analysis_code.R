#Preamble #######
packages <- c("esc",
              "data.table",
              "lme4",
              "lmerTest",
              "vroom",
              "meta",
              "metafor",
              "clubSandwich",
              "ggplot2",
              "tidyverse",
              "hablar",
              "janitor",
              "tidylog",
              "dmetar",
              "metaviz",
              "orchaRd",
              "patchwork",
              "RColorBrewer",
              "magrittr",
              "broom",
              "broom.mixed",
              "ggeasy",
              "ggstance",
              "metaDigitise",
              "klippy",
              "emmeans",
              "gt",
              "ggbeeswarm",
              "MuMIn",
              "vroom",
              "orchaRd",
              "metaAidR",
              "shinyDigitise",
              "devtools",
              "ggtext")

lapply(packages, require, character.only = TRUE)

#Pre1 - Extract data with metadigitise or shiny#########

data <- shinyDigitise(dir = "Desktop")
data <- metaDigitise(dir = "Desktop")

#Pre2 - Summary data########

#summary data for raw
Linddat<-read_csv(here::here("ExtractData", "Raw", "LindRaw.csv")) %>%
  janitor::clean_names() %>%
  group_by(selection_treatment,sex) %>%
  summarise(count = n())

Angelldat<-read_csv(here::here("ExtractData", "Raw", "AngellRaw.csv")) %>%
  janitor::clean_names() %>%
  filter(sex == "M") %>%
  group_by(male_parent) %>%
  summarise(count = n())

Dowlingdat<-read_csv(here::here("ExtractData", "Raw", "DowlingRaw.csv")) %>%
  janitor::clean_names() %>%
  group_by(treatment, sex) %>%
  summarise(count = n())

Bouwhuisdat<-read_csv(here::here("ExtractData", "Raw", "BouwhuisRaw.csv")) %>%
  janitor::clean_names() %>%
  group_by(sex) %>%
  summarise(count = n())

samp <- dat_tot[,3] %>%
  rename(count = Sample_Size)

samp_tot<-data.frame(a = bind_rows(samp, Linddat[,3], Angelldat[,2], Dowlingdat[,3], Bouwhuisdat[,2]))
median(samp_tot$count)


#Pre3 - Beta ########

#optim function
int = .1
get.slope <- function(dat){
  
  WSS.optim <- function(par){
    SSy = 0
    y.hat = par[1] + par[2]*dat$x
    for(i in 1:nrow(dat)){
      range = seq(-6 *dat$y.sd[i] + dat$y.mean[i], 6 *dat$y.sd[i] +dat$y.mean[i], by = int)
      for (j in range){
        SSy = SSy + dat$n[i]*dnorm(x = j, mean = dat$y.mean[i], sd=dat$y.sd[i])*int*(j-y.hat[i])^2}
      
    }
    SSy}
  
  result <- optim(par = c(1,2), fn = WSS.optim)
  
  xbar <- weighted.mean(x = dat$x, w = dat$n)
  SSx = 0
  for(i in 1:nrow(dat)){SSx = SSx + dat$n[i]*(dat$x[i] - xbar)^2}
  
  num = result$value/(sum(dat$n)-2)
  SE = sqrt(num/SSx)
  
  data.frame(estimate = result$par[2], SE = SE, study = dat$studytreat[1])}


#create list of all data for optim for humans

file_list <- list.files(path = here::here("ExtractData", "Non-human"),pattern = ".csv")

setwd(here::here("ExtractData", "Non-human"))

dat_tot<-vroom(file_list, id = "Study") %>%
  mutate(StudyTreat = paste(Study, "-", Treat))

datout <- split(dat_tot, dat_tot$StudyTreat)

setwd(here::here())

#create list of all data for optim for non-humans

file_listhum <- list.files(path = here::here("ExtractData", "Human"),pattern = ".csv")

setwd(here::here("ExtractData", "Human"))

dat_tothum<-vroom(file_listhum, id = "Study") %>%
  mutate(StudyTreat = paste(Study, "-", Treat))

datouthum <- split(dat_tothum, dat_tothum$StudyTreat)

setwd(here::here())

#optim for data number 3 - set columns names
colnames= c("file", "x","n", "y.mean","y.sd", "y.se", "treat", "studytreat")
datout<-lapply(datout, setNames, colnames) 

colnames= c("file", "x","n", "y.mean","y.sd", "y.se", "treat", "studytreat")
datouthum<-lapply(datouthum, setNames, colnames) 

#function to create terminal ages for humans and non-human
extract.age <- function(dat){
  x<-nrow(dat)
  y<-x-1
  
  dat <- dat[c(y:x),]
}


extract.humans <- function(dat){
  
  dat <- dat %>% filter(x > 26.3)
  
}

datoutterm <- lapply(datout, extract.age)
datouttermhum <- lapply(datouthum, extract.humans)

#get beta values
#for non-humans and humans
mat_age<-lapply(datout, get.slope)
mat_agehum<-lapply(datouthum, get.slope)

#terminal
mat_ageterm<-lapply(datoutterm, get.slope)
mat_agetermhum<-lapply(datouttermhum, get.slope)

#bind and save
mat_age2<-bind_rows(mat_age)
mat_agehum2<-bind_rows(mat_agehum)
mat_age2term<-bind_rows(mat_ageterm)
mat_age2termhum<-bind_rows(mat_agetermhum)

write.csv(mat_age2, "mat_age.csv", row.names = F)
write.csv(mat_agehum2, "mat_agehum.csv", row.names = F)
write.csv(mat_age2term, "mat_ageterm.csv", row.names = F)
write.csv(mat_age2termhum, "mat_agetermhum.csv", row.names = F)

#for raw data - split by treatment/subset fro terminal

#lind
Lind<-read_csv(here::here("ExtractData", "Raw", "LindRaw.csv")) %>%
  janitor::clean_names() %>%
  convert(fct(selection_treatment, line, maternal_id, sex)) %>%
  mutate(treat = paste(selection_treatment, "-",sex)) %>%
  mutate(mat_age2 = maternal_age^2) %>%
  #terminal?
  #filter(maternal_age == 3| maternal_age == 4) %>%
  group_split(treat, keep = T)
  
mod_outL <- matrix(nrow = length(Lind), ncol = 3)

for (i in 1:length(Lind)){
  mod<-lmer(lifespan ~ maternal_age + (1|maternal_id), data = Lind[[i]])
  summary(mod)
  print(i)
  
  mod_outL[i,1:3] <- c(summary(mod)$coef[2,1],
                      summary(mod)$coef[2,2],
                       Lind[[i]]$treat[1])
  
}

write.csv(mod_outL, "mat_ageLindTerm.csv", row.names = F)

#Angell (Two)
Angell<-read_csv(here::here("ExtractData", "Raw", "AngellRaw.csv")) %>%
  janitor::clean_names() %>%
  filter(sex == "M") %>%
  mutate(mat_age = case_when(female_parent == "YF" ~ 2,
                             female_parent == "OF" ~ 12)) %>%
  convert(fct(male_parent, block, sex), num(mat_age)) %>%
  group_split(male_parent, keep = T)

mod_outA <- matrix(nrow = length(Angell), ncol = 3)

for (i in 1:length(Angell)){
  mod<-lmer(lifespan ~ mat_age + (1|block), data = Angell[[i]])
  summary(mod)
  print(i)
  
  mod_outA[i,1:3] <- c(summary(mod)$coef[2,1],
                       summary(mod)$coef[2,2],
                       as.character(Angell[[i]]$male_parent[1]))
  
}

write.csv(mod_outA, "mat_ageAngell.csv", row.names = F)

#Dowling (Two)
Dowling<-read_csv(here::here("ExtractData", "Raw", "DowlingRaw.csv")) %>%
  janitor::clean_names() %>%
  mutate(mat_age = case_when(day == "D5" ~ 5,
                             day == "D10" ~ 10)) %>%
  mutate(treat = paste(treatment, "-", sex)) %>%
  convert(fct(treat, vial, sex), num(mat_age, days)) %>%
  group_split(treat, keep = T)

mod_outD <- matrix(nrow = length(Dowling), ncol = 3)

for (i in 1:length(Dowling)){
  mod<-lmer(days ~ mat_age + (1|vial) + (1|matrep), data = Dowling[[i]])
  summary(mod)
  print(i)
  
  mod_outD[i,1:3] <- c(summary(mod)$coef[2,1],
                       summary(mod)$coef[2,2],
                       as.character(Dowling[[i]]$treat[1]))
  
}

write.csv(mod_outD, "mat_ageDowl.csv", row.names = F)

#Bouwhuis
Bouwhuis<-read_csv(here::here("ExtractData", "Raw", "BouwhuisRaw.csv")) %>%
  janitor::clean_names() %>%
  convert(fct(sex)) %>%
  mutate(mat_age2 = maternal_age^2) %>%
  #terminal?
  #filter(maternal_age > 10) %>%
  group_split(sex, keep = T) 

mod_outB <- matrix(nrow = length(Bouwhuis), ncol = 3)

for (i in 1:length(Bouwhuis)){
  mod<-lmer(lifespan ~ d_maternal_age + av_maternal_age + (1|mum_id), data = Bouwhuis[[i]])
  summary(mod)
  print(i)
  
  mod_outB[i,1:3] <- c(summary(mod)$coef[2,1],
                       summary(mod)$coef[2,2],
                       as.character(Bouwhuis[[i]]$sex[1]))
  
}
# 2 = female
# 1 = male

write.csv(mod_outB, "mat_ageBouTerm.csv", row.names = F)

#Kroeger (cant as standardised)

Kroeger<-read_csv(here::here("ExtractData", "Raw", "KroegerRaw.csv")) %>%
  janitor::clean_names() %>%
  convert(fct(valley)) %>%
  #filter(maternal_age > 10) %>%
  group_split(valley, .keep = T) 

mod_outK <- matrix(nrow = length(Kroeger), ncol = 3)

i <- 1
for (i in 1:length(Kroeger)){
  mod<-lmer(lifespan ~ mab + (1|mother_id) + (1|yrborn), data = Kroeger[[i]])
  summary(mod)
  print(i)
  
  mod_outK[i,1:3] <- c(summary(mod)$coef[2,1],
                       summary(mod)$coef[2,2],
                       as.character(Kroeger[[i]]$valley[1]))
  
}

write.csv(mod_outK, "mat_ageKro.csv", row.names = F)

#IC 
IC1<-read_csv("ExtractData/Raw/ICRaw.csv") %>%
  janitor::clean_names() %>%
  #teerminal
  filter(care_egg == 2) %>%
  filter(care_age == 8|care_age==11) %>%
  convert(fct(age_factor, block, dam_name), num(care_age, death_age))

m1 <- lmer(death_age ~ age_mated + carcass_weight + age_at_death + (1|block/dam_name), data = IC1)
summary(m1)

IC2<-read_csv("ExtractData/Raw/ICRaw.csv") %>%
  janitor::clean_names() %>%
  filter(care_age == 2) %>%
  filter(care_egg == 8|care_egg==11)%>%
  convert(fct(age_factor, block, dam_name), num(care_age, death_age))

m1 <- lmer(death_age ~ age_egg + carcass_weight + age_factor + (1|block/dam_name), data = IC2)
summary(m1)

#Import data######

articleMeanB <- 
  read_csv("metaData.csv") %>%
  janitor::clean_names() %>%
  mutate(bvar = bse^2) %>%
  mutate(termvar = post_se^2) %>%
  mutate(wi = 1/bvar) %>%
  mutate(prec = 1/bse) %>%
  mutate(teffect = b/sqrt(bvar)) %>%
  as.data.frame() %>%
  mutate(year.c = scale(year, scale = F))  #robustness with age classes?
  #filter(number_of_age_classes > 2)

#######MODELS##########
#1. Random effect model######
#all ages
model_age <- rma.mv(yi = b,
                    V = bvar,
                    mods = ~1,
                    random = ~1 | species/study/replicate,
                    test = "t",
                    method = "REML",
                    data=articleMeanB)

summary(model_age)
tidy(model_age)

#overall orchard plot
orchard_plot(model_age, mod = "1", group = "study", data = articleMeanB, xlab = "coef")

#terminal 
model_ageterm <- rma.mv(yi = post_t,
                     V = termvar,
                     mods = ~1,
                     random = ~1 | species/study/replicate,
                     test = "t",
                     method = "REML",
                     data=articleMeanB)

summary(model_ageterm)
tidy(model_ageterm)

#2. Pub bias ######

#all
model_age_all <- rma.mv(yi = b,
                        V = bvar,
                        mods = ~1 + year.c + bvar,
                        random = ~1 | species/study/replicate,
                        test = "t",
                        data=articleMeanB,
                        method = "REML")

summary(model_age_all)

#terminal
model_ageterm <- rma.mv(yi = post_t,
                        V = termvar,
                        mods =  ~1 + year.c + termvar,
                        random = ~1 | species/study/replicate,
                        test = "t",
                        data=articleMeanB,
                        method = "REML")

summary(model_ageterm)
tidy(model_ageterm)




#3. Moderators #######

#remove rodents from terminal as no effect size
termB <- articleMeanB %>%
  filter(!is.na(post_t))

#models
allage_pub <- rma.mv(yi = b,
                           V = bvar,
                           mods = ~ 1 + group + offspring_sex + pac + year.c + bvar,
                           random = ~1 | species/study/replicate,
                           test = "t",
                           data=articleMeanB,
                     method = "REML")
summary(allage_pub)

termage_pub <- rma.mv(yi = post_t,
                            V = termvar,
                            mods = ~ 1 + group + offspring_sex + pac + year.c + termvar,
                            random = ~1 | species/study/replicate,
                            test = "t",
                            data=termB,
                      method = "REML")
summary(termage_pub)

allage <- rma.mv(yi = b,
                 V = bvar,
                            mods = ~ 1 + group + offspring_sex + pac,
                            random = ~1 | species/study/replicate,
                            test = "t",
                            data=articleMeanB,
                 method = "REML")

summary(allage)


termage <- rma.mv(yi = post_t,
                            V = termvar,
                            mods = ~ 1 + group + offspring_sex + pac ,
                            random = ~1 | species/study/replicate,
                            test = "t",
                            data=termB,
                  method = "REML")

summary(termage)

resa <- qdrg(object = allage_pub, data = articleMeanB, at = list(bvar = 0, year.c = 0)) 

em1a<-emmeans(resa, specs = "group", type = "response", adjust = "none")
em2a<-emmeans(resa, specs = "offspring_sex", type = "response", adjust = "none")
em3a<-emmeans(resa, specs = "pac", type = "response", adjust = "none")

resb <- qdrg(object = termage_pub, data = termB,  at = list(termvar = 0, year.c = 0)) 

em1b<-emmeans(resb, specs = "group", type = "response", adjust = "none")
em2b<-emmeans(resb, specs = "offspring_sex", type = "response", adjust = "none")
em3b<-emmeans(resb, specs = "pac", type = "response", adjust = "none")

resc <- qdrg(object =allage , data = articleMeanB) 

em1c<-emmeans(resc, specs = "group", type = "response", adjust = "none")
em2c<-emmeans(resc, specs = "offspring_sex", type = "response", adjust = "none")
em3c<-emmeans(resc, specs = "pac", type = "response", adjust = "none")

resd <- qdrg(object = termage, data = termB) 

em1d<-emmeans(resd, specs = "group", type = "response", adjust = "none")
em2d<-emmeans(resd, specs = "offspring_sex", type = "response", adjust = "none")
em3d<-emmeans(resd, specs = "pac", type = "response", adjust = "none")

#a = betapub, b = termpub, c = termnopub, d = betanopub
em1a.tidy<-tidy(em1a, conf.int = T)
em1b.tidy<-tidy(em1b, conf.int = T)
em1c.tidy<-tidy(em1c, conf.int = T)
em1d.tidy<-tidy(em1d, conf.int = T)

em1a.tidy$term <- "All Ages"
em1a.tidy$pubadjust <- "Adjusted"
em1b.tidy$term <- "Terminal Ages"
em1b.tidy$pubadjust <- "Adjusted"
em1c.tidy$term <- "All Ages"
em1c.tidy$pubadjust <- "Unadjusted"
em1d.tidy$term <- "Terminal Ages"
em1d.tidy$pubadjust <- "Unadjusted"

valuedat <- rbind(em1a.tidy, em1b.tidy,em1c.tidy,em1d.tidy)
valuedat$grouping <- paste(valuedat$term, "-", valuedat$pubadjust)

valuedat %<>% convert(fct(grouping)) %>%
  mutate(grouping = fct_relevel(grouping, 
                            "All Ages - Unadjusted", "All Ages - Adjusted", "Terminal Ages - Unadjusted", 
                            "Terminal Ages - Adjusted")) %>%
  mutate(grouping = fct_rev(grouping))

plotdat <- 
  as.data.frame(cbind(estimate = articleMeanB$b,
                      group = articleMeanB$group,
                      size = 1/articleMeanB$bse)) %>%
  convert(num(estimate, size), fct(group))

totdat<-plotdat %>% group_by(group) %>%
  summarise(count = n())

join1 <- rep("<br>*n*",8)
join2 <- rep("=", 8)
my.labs1 <- c(as.character(totdat$group[1]),
              as.character(totdat$group[2]),
              as.character(totdat$group[3]),
              as.character(totdat$group[4]),
              as.character(totdat$group[5]),
              as.character(totdat$group[6]),
              as.character(totdat$group[7]),
              as.character(totdat$group[8]),
              as.character(totdat$group[9]))

my.labs1 <- gsub("-|\\s+|^$"," ",my.labs1)

my.labs2 <- c(as.character(totdat$count[1]),
              as.character(totdat$count[2]),
              as.character(totdat$count[3]),
              as.character(totdat$count[4]),
              as.character(totdat$count[5]),
              as.character(totdat$count[6]),
              as.character(totdat$count[7]),
              as.character(totdat$count[8]),
              as.character(totdat$count[9]))

my.labels <- paste(my.labs1, join1, join2, my.labs2)

#moderator graph
ggplot2::ggplot(valuedat, aes(x = estimate, y = group, colour = grouping, linetype = grouping)) +
  ggplot2::geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),  height = 0.3,size = 1.2, position = position_dodgev(height = 0.6)) +
  ggplot2::geom_point(position = position_dodgev(height = 0.6), size = 2) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.4) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.title = element_text(size = 9)) +
  ggplot2::theme(legend.direction="horizontal") +
  ggplot2::theme(legend.background = element_blank()) +
  ggplot2::labs(x = "Effect Size", y = "", colour = "Model grouping", linetype = "Model grouping") +
  guides(linetype = guide_legend(nrow = 4), colour = guide_legend(nrow = 4)) +
  scale_y_discrete(labels = my.labels) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_markdown(hjust=0.5, size=12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.box="vertical",
        legend.key.width = unit(5, "line"),
        legend.position="bottom") +
  scale_colour_manual(values=c("black","black", "grey44","grey44"),
                      breaks = c("All Ages - Unadjusted", "All Ages - Adjusted", "Terminal Ages - Unadjusted", 
                                 "Terminal Ages - Adjusted"),
                      labels = c("All Ages - Unadjusted", "All Ages - Adjusted", "Terminal Ages - Unadjusted", 
                                 "Terminal Ages - Adjusted")) +
  scale_linetype_manual(values=c("solid", "dashed", "solid", "dashed"),
                        breaks = c("All Ages - Unadjusted", "All Ages - Adjusted", "Terminal Ages - Unadjusted", 
                                   "Terminal Ages - Adjusted"),
                        labels = c("All Ages - Unadjusted", "All Ages - Adjusted", "Terminal Ages - Unadjusted", 
                                   "Terminal Ages - Adjusted")) +
  ggplot2::geom_hline(yintercept = seq(from = 0.5, to = 8.5, by = 1), linetype = 1, colour = "black", alpha = 0.6) +
  geom_text(data = valuedat,
             aes(label = sprintf("%0.2f", round(estimate,2))),
             vjust = -0.5,
             position = position_dodgev(height = 0.6),
             size = 3) +
  geom_text(data = valuedat,
            aes(x = conf.high, y = group, label = sprintf("%0.2f", round(conf.high,2))),
            vjust = -0.5,
            position = position_dodgev(height = 0.6),
            size = 3)  +
  geom_text(data = valuedat,
            aes(x = conf.low, y = group, label = sprintf("%0.2f", round(conf.low,2))),
            vjust = -0.5,
            position = position_dodgev(height = 0.6),
            size = 3) 

ggplot2::ggsave("mod_fig.tiff", width = 15, height = 12, device = "tiff", dpi= 300)
