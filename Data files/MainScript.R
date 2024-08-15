### Aiden R Script for skua data

#libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(lme4)
library(emmeans)
library(DHARMa)

# colour coordination
c1 <- c("Seabird" = "#440154", "Marine mammal" = "#414487", "Fish" = "#2a788e", "Penguin" = "#22a884", "Rodent" = "#7ad151", "Livestock" = "#fde725", "Other" = "grey")
c2 <- c("12S" = "#fe9f6d", "MarVer"="#8c2981")
c3 <- c("Seabird" = "#440154", "Fish" = "#2a788e", "Rabbit" = "#7ad151", "Prey" = "#fde725", "Other" = "grey")
#size 995, 700 figures

### PART 1: BIAS

bias <- read.csv("bias.csv")
bias$primer <- gsub( " .*$", " ", bias$Primer.mix)
bias$primer <- str_trim(bias$primer)
bias$primer <- factor(bias$primer, levels=c("16S-Vert", "12S-Vert", "MarVer", "FISH"))
bias$Group <- factor(bias$Group, levels=c("Seabird", "Marine mammal", "Fish", "Penguin", "Rodent", "Livestock"))

fig2a <- ggplot(bias, aes(x=primer, y=Expected.., fill=Group)) + geom_bar(position="fill", stat="identity") + theme_classic() +
  scale_fill_manual(values = c1) + labs(y = "Proportion reads", x = "", fill = "Species group") + theme(text = element_text(size = 20),axis.text.x = element_text(angle=90, vjust =.5, hjust=1))
fig2b <- ggplot(bias, aes(x=primer, y=Observed.., fill=Group)) + geom_bar(position="fill", stat="identity") + theme_classic() +
  scale_fill_manual(values = c1) + labs(y = "Proportion reads", x = "", fill = "Species group") + theme(text = element_text(size = 20),axis.text.x = element_text(angle=90, vjust =.5, hjust=1))

### PART 2: SWAB SAMPLES

swab_full <- read.csv("Swabs.csv")
swab <- swab_full %>%
  group_by(group, Primer, Pellet_.ID) %>%
  summarise(reads = sum(readnr))
swab <- swab %>%
  drop_na(reads)

swab$group <- factor(swab$group, levels=c("Seabird", "Marine mammal", "Other"))

swabch <- swab %>%
  select(group, Primer, reads, Pellet_.ID)
swabch <- pivot_wider(swabch, names_from = Primer, values_from = reads)

fig3 <- ggplot(swab, aes(x = Pellet_.ID, y = reads, fill = group)) + geom_bar(position="fill", stat="identity") + theme_classic() +
  scale_fill_manual(values = c1) + 
  labs(y = "Proportion reads", x = "", fill = "Species group") + theme(text = element_text(size = 20),axis.text.x = element_text(angle=90, vjust =.5, hjust=1)) + facet_wrap(~Primer)

swabch <- swabch %>%
  group_by(Pellet_.ID) %>%
  mutate(perc16 = MarVer/sum(MarVer), perc12 = `12S`/sum(`12S`))

swabch %>% group_by(group) %>% summarise(mean(perc16), mean(perc12))

chisq.test(swabch$MarVer,swabch$'12S')

swabsp <- swab_full %>%
  group_by(Specific, Primer, Pellet_.ID) %>%
  summarise(reads = sum(readnr))
swabsp <- swabsp %>%
  drop_na(reads)

fig4 <- ggplot(swabsp, aes(x = Specific, y = log(reads), fill = Primer)) + geom_bar(position="dodge", stat = "identity" ) + theme_classic() +
  scale_fill_manual(values = c2) + 
  labs(y = "Logged reads", x = "", fill = "Primer") + theme(text = element_text(size = 20), axis.text.x = element_text(angle=90, vjust =.5, hjust=1)) + facet_wrap(~Pellet_.ID)

### PART 3: PELLETS

pel1 <- read.csv("pellet1.1.csv")
pel2 <- read.csv("pellet2.1.csv")
pel2 <- pel2 %>%
  drop_na(reads)
pel1 <- pel1 %>%
  drop_na(reads)

pel2 <- pel2 %>% select(-c(spc, barcodenr, readnr))
pel1 <- pel1 %>% select(-c(spec,barcode.1))

pels <- rbind(pel1, pel2)
pels$group <- factor(pels$group, levels=c("Seabird", "Fish", "Rabbit", "Prey", "Other"))

pelgroup <- pels %>% group_by(group,Primer,Pellet_.ID) %>%
  summarise(reads = sum(reads), expect = unique(exp_group), observe = unique(obs_group))
pelgroup <- pelgroup %>% mutate(expect = ifelse(group=="Other", 0, expect))

pelspec <- pels %>% group_by(Specific,Primer,Pellet_.ID) %>%
  summarise(reads = sum(reads), observe = sum(obs_perc), expect = unique(exp_ind), observe = unique(exp_ind))

pelloc <- pels %>% group_by(Pellet_.ID) %>%
  summarise(loc = unique(Location))
peltret <- pels %>% group_by(Pellet_.ID) %>%
  summarise(tret = unique(Extract))
groups <- pels %>% group_by(group) %>%
  summarise(specicific = unique(Specific))

pelspec <- left_join(pelspec,pelloc)
pelgroup <- left_join(pelgroup, pelloc)
pelspec <- left_join(pelspec, peltret)
pelgroup <- left_join(pelgroup,peltret)

pelboth <- pelgroup %>% subset(Pellet_.ID=="2-4" | Pellet_.ID=="3-44")
pelboth2 <- pelspec %>% subset(Pellet_.ID=="2-4" | Pellet_.ID=="3-44")

fig5 <- ggplot(pelboth2, aes(x = Specific, y = log(reads), fill = Primer)) + geom_bar(position="dodge", stat = "identity" ) + theme_classic() +
  scale_fill_manual(values = c2) + 
  labs(y = "Logged reads", x = "", fill = "Primer") + theme(text = element_text(size = 20), axis.text.x = element_text(angle=90, vjust =.5, hjust=1)) + facet_wrap(~Pellet_.ID)

pelch <- pelboth %>%
  select(group, Primer, reads, Pellet_.ID)
pelch <- pivot_wider(pelch, names_from = Primer, values_from = reads)
chisq.test(pelch$MarVer, pelch$'12S')

pelgroup <- pelgroup %>% subset(Primer == "MarVer" & Pellet_.ID != "5-10")
pelnoo <- pelgroup %>% subset(group != "Other")
pelnoosp <- pelspec %>% subset(Primer == "MarVer" & Pellet_.ID != "5-10" & Specific != "Bacteria")
pelspec <- pelspec %>% subset(Primer == "MarVer" & Pellet_.ID != "5-10")

fig6 <- ggplot(pelgroup, aes(x = Pellet_.ID, y = reads, fill = group)) + geom_bar(position="fill", stat="identity") + theme_classic() +
  scale_fill_manual(values=c3) + labs(y = "Proportion reads", x = "", fill = "Species group") + theme(text = element_text(size = 20),axis.text.x = element_text(angle=90, vjust =.5, hjust=1))

fig7 <- ggplot(pelnoosp, aes(x = Specific, y = observe, fill = Pellet_.ID)) + geom_bar(position="dodge", stat = "identity" ) + theme_classic() +
  scale_fill_viridis_d(option="magma") + labs(y = "Logged reads", x = "", fill = "Pellet") + theme(text = element_text(size = 20),axis.text.x = element_text(angle=90, vjust =.5, hjust=1))

fig11 <- ggplot(pelnoo, aes(x = loc, y = reads, fill = group, )) + geom_bar(position="fill", stat="identity") + theme_classic() +
  scale_fill_manual(values = c3) + labs(y = "Proportion reads", x = "", fill = "Species group") + theme(text = element_text(size = 20),axis.text.x = element_text(angle=90, vjust =.5, hjust=1))

peloc <- pelnoo %>% group_by(loc, group) %>% summarise(reads = sum(reads))
peloc <- pivot_wider(peloc, names_from = loc, values_from = reads)
peloc <- peloc %>% remove_rownames %>% column_to_rownames(var = "group")
peloc[is.na(peloc)] <- 0
chisq.test(peloc)

### MODEL CREATION

epsilon <- 1e-8
pelglm <- pelgroup %>% 
  mutate(value = (observe-expect)/(expect+epsilon))
pelglm2 <- pelglm %>% subset(group != "Other")
pelnglm <- pelgroup %>%
  subset(group == "Other")

glm_morph2 <- lmer(value ~ group+(1|loc)+(1|tret), data = pelglm2)
glm_loc2 <- lmer(value ~ group*loc+group+loc + (1|tret), data = pelglm2)
glm_tret <- lmer(observe~tret+(1|loc), data = pelglm)

sim3 <- simulateResiduals(fittedModel = glm_tret)
sim4 <- simulateResiduals(fittedModel = glm_morph2)
sim5 <- simulateResiduals(fittedModel = glm_loc2)

plotQQunif(sim3)
plotQQunif(sim4)
plotQQunif(sim5)

summary(glm_tret)
summary(glm_morph2)
summary(glm_loc2)

ems1 <- emmeans(glm_morph2, ~group)
ems2 <- emmeans(glm_loc2, ~group+loc)
ems3 <- emmeans(glm_tret, ~tret)

pairs(ems1)
pairs(ems2)
pairs(ems3)

### MODEL FITTING

fix1 <- tidy(glm_morph2, effects = "fixed")
fix1$name <- c("Seabird", "Fish","Prey")
fix1$order <- factor(fix1$name, levels=c("Fish", "Seabird", "Prey","Other"))
ggplot(fix1, aes(x = order, y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0.2) +
  labs(x = "",
       y = "") +
  theme_classic() +
  coord_flip() + theme(text = element_text(size = 20))

fix2 <- tidy(glm_loc2, effects = "fixed")
fix2$name <- c("Seabird @ Foula", "Fish @ Foula", "Prey @ Foula", "Seabird @ Noss", "Seabird @ St Kilda", "Fish @ Noss", "Fish @ St Kilda", "Prey @ St Kilda")
ggplot(fix2, aes(x = name, y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0.2) +
  labs(x = "",
       y = "") +
  theme_classic() +
  coord_flip() + theme(text = element_text(size = 20))

fix3 <- tidy(glm_tret, effects = "fixed")
fix3$name <- c("Bleach and PBS", "No wash", "PBS twice")
ggplot(fix3, aes(x = name, y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0.2) +
  labs(x = "",
       y = "Proportion reads") +
  theme_classic() +
  coord_flip()+ theme(text = element_text(size = 20))

### superseded

#contam <- subset(pelgroup, group == "Other")
#contam <- subset(contam, Primer == "MarVer")
#contplot <- ggplot(contam, aes(x = tret, y = observe, fill = tret, )) + geom_boxplot() + theme_classic() +
#  scale_fill_viridis_d() + labs(y = "Proportion reads", x = "", fill = "Extraction") + theme(axis.text.x = element_text(angle=90, vjust =.5, hjust=1))
#consum <- contam %>% group_by(tret) %>% summarise(sum = sum(reads))
#chisq.test(consum$sum, p = c(1/3,1/3,1/3))

#logLik(glm_loc)
#logLik(glm_loc1)
#2*(-313.845--425.6033)
#1-pchisq(223.5166,6)

#logLik(glm_morph)
#2*(-293.3822--420.3738)
#1-pchisq(253.9832,6)

#glm_loc3 <- lmer(value ~ group*loc+tret+group+loc + (1|tret), data = pelglm2)
#glm_loc1 <- lmer(value ~ group+loc + (1|tret), data = MarVer) sigdif from above model, so do not use
#glm_loc <- lmer(value ~ group*loc+group+loc + (1|tret), data = pelglm)
#glm_morph <- lmer(value ~ group+(1|loc)+(1|tret), data = pelglm)

#pelbothfig <- ggplot(pelboth, aes(x = Pellet_.ID, y = reads, fill = group)) + geom_bar(position="fill", stat="identity") + theme_classic() +
#  scale_fill_manual(values=c3) + labs(y = "Proportion reads", x = "Pellet", fill = "Species group") + theme(text = element_text(size = 20),axis.text.x = element_text(angle=90, vjust =.5, hjust=1)) + facet_wrap(~Primer)

#pelbothfige <- ggplot(pelboth, aes(x = Pellet_.ID, y = expect, fill = group)) + geom_bar(position="fill", stat="identity") + theme_classic() +
#  scale_fill_manual(values=c3) + labs(y = "Proportion reads", x = "", fill = "Species group") + theme(axis.text.x = element_text(angle=90, vjust =.5, hjust=1)) + facet_wrap(~Primer)
