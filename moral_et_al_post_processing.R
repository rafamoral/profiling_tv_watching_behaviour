## Profiling Customer Television Watching Behaviour Using A Novel Bayesian
## Hierarchical Joint Model for Time-to-Event and Count Data
## by Moral et al.

## The code in this script file reproduces the following analyses and figures:
## 1. clustering
## 2. correlation plots
## 3. exploratory plots

## loading libraries and helper functions
library(tidyverse)
library(ggplot2)
library(ape)
library(corrplot)
library(dendextend)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## loading model summaries
load("model_summary.RData")

## parameter estimates
delta_summary
gamma_summary
phi_summary
psi_summary
d_summary
g_summary

psi_hat <- psi_summary$mean
phi_hat <- phi_summary$mean[1:40]
d_hat <- d_summary$mean
g_hat <- matrix(g_summary$mean, ncol = 8, byrow = TRUE)

customer_parms <- data.frame(psi_hat, phi_hat, d_hat, g_hat)
names(customer_parms) <- c("psi","phi","d",paste0("b",1:8))

## preparing data for naive approach
load("raw_data_naive_approach.RData")

## clustering
rownames(customer_parms) <- c(paste0("A:", 1:20), paste0("C:", 1:20))
rownames(dc_cc_mat) <- c(paste0("A:", 1:20), paste0("C:", 1:20))

## clustering based on model output
dist_matrix <- dist(scale(customer_parms))

hc <- hclust(dist_matrix, "ward.D2")
clus2 <- cutree(hc, 2)
cols <- gg_color_hue(2)

## clustering based on naive approach
dist_matrix <- dist(scale(dc_cc_mat))

hc_raw <- hclust(dist_matrix, "ward.D2")
clus2_raw <- cutree(hc_raw, 2)
cols <- gg_color_hue(2)

## fan-type dendrograms
par(mfrow = c(1,2))
plot(as.phylo(hc_raw), type = "fan",
     tip.color = cols[rep(1:2, each = 20)],
     font = 2,
     edge.width = .7, main = "(a) Raw data")
plot(as.phylo(hc), type = "fan",
     tip.color = cols[rep(1:2, each = 20)],
     font = 2,
     edge.width = .7, main = "(b) Model parameter estimates")

dend_left <- color_branches(as.dendrogram(hc_raw), k = 2, col = gg_color_hue(3)[c(2,3)])
dend_right <- color_branches(as.dendrogram(hc), k = 2, col = gg_color_hue(3)[c(2,3)])

d <- dendlist(color_labels(dend_left, col = "red", labels = labels(dend_left)[grep("C", labels(dend_left))]),
              color_labels(dend_right, col = "red", labels = labels(dend_right)[grep("C", labels(dend_right))]))

png("clustering.png", w = 10, h = 6, res = 800, units = "in")
tanglegram(d, sort = TRUE,
           common_subtrees_color_lines = FALSE,
           highlight_distinct_edges = FALSE,
           highlight_branches_lwd = FALSE,
           main_left = "Raw data", main_right = "Parameter estimates")
dev.off()

## correlation matrix for all customers
corrplot(cor(customer_parms), method = "square", type = "lower", diag = FALSE)

cor_mixed <- cor1 <- cor(customer_parms[grep("A", rownames(customer_parms)),])
cor2 <- cor(customer_parms[grep("C", rownames(customer_parms)),])

cor_mixed[upper.tri(cor1)] <- cor2[upper.tri(cor2)]

png("corrplots_sample_data.png", res = 800, units = "in", w = 6, h = 7)
corrplot(cor_mixed, method = "square", diag = FALSE,
         tl.srt = 0, tl.offset = 1, tl.col = "white", cl.pos = "b",
         mar = c(0,2,2,0))
text(1:11, 12, expression(psi, varphi, d, b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8]))
text(0, 11:1,  expression(psi, varphi, d, b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8]))
text(-.75, 6, "Active", adj = .5, srt = 90, xpd = NA, font = 2)
text(6, 12.75, "Cancelled", adj = .5, xpd = NA, font = 2)
dev.off()

## exploratory plots
customer_parms$status <- ifelse(substr(rownames(customer_parms), 1, 1) == "A",
                                "active", "cancelled") %>% as.factor
customer_parms2 <- customer_parms
customer_parms2[,1:11] <- scale(customer_parms2[,1:11])

png("estimates.png", res = 800, units = "in", w = 6, h = 5)
customer_parms %>%
  pivot_longer(1:11,
               names_to = "Parameter",
               values_to = "Estimate") %>%
  mutate(Parameter = factor(Parameter, levels = unique(Parameter))) %>%
  ggplot(aes(x = Parameter, y = Estimate, fill = status, col = status)) +
  theme_bw() +
  geom_boxplot(size = .4, outlier.size = .75) +
  labs(fill = "Subscription status", col = "Subscription status") +
  scale_x_discrete(labels = c(expression(psi),expression(varphi),"d",
                              expression(b[1]),expression(b[2]),expression(b[3]),
                              expression(b[4]),expression(b[5]),expression(b[6]),
                              expression(b[7]),expression(b[8]))) +
  scale_fill_manual(values = paste0(gg_color_hue(2), "55")) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom")
dev.off()

## scaled version
customer_parms2 %>%
  pivot_longer(1:11,
               names_to = "Parameter",
               values_to = "Estimate") %>%
  mutate(Parameter = factor(Parameter, levels = unique(Parameter))) %>%
  ggplot(aes(x = Parameter, y = Estimate, fill = status, col = status)) +
  theme_bw() +
  geom_boxplot(size = .4, outlier.size = .75) +
  labs(fill = "Subscription status", col = "Subscription status") +
  scale_x_discrete(labels = c(expression(psi),expression(varphi),"d",
                              expression(b[1]),expression(b[2]),expression(b[3]),
                              expression(b[4]),expression(b[5]),expression(b[6]),
                              expression(b[7]),expression(b[8]))) +
  scale_fill_manual(values = paste0(gg_color_hue(2), "55")) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom") +
  ylab("Scaled estimate")