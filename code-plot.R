library(ggplot2)
library(ggrepel)
library(dplyr)
devtools::source_gist("07db45251922fdf3bdb500ec9493c1bf") #this gets my plotting theme from the internet

params <- c("n", "K", "ymin", "ymax", "fold_change")

for(p in params){
  data <- read.csv(paste("", p, "_morris.csv", sep = ""))
  # data$X <- factor(data$X, levels = data$X)
  data <- arrange(data, desc(mu_star))

  ggplot(data) +
    geom_pointrange(aes(X, y = mu_star,
                        ymin = mu_star - mu_star_conf,
                        ymax = mu_star + mu_star_conf),
                    size = 0.2) +
    scale_x_discrete("", limits = data$X) +
    scale_y_continuous(expression(paste(mu, '*'))) +
    ggtitle(p) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 270, vjust = 0, hjust=0))
  ggsave(paste("sensitivity_", p, "_mu-star.pdf", sep =""), width = 80, height = 60, units = "mm")

  ggplot(data, aes(y = sigma, x = mu)) +
    geom_point(size = 0.5) +
    geom_text_repel(aes(label = X),
                    min.segment.length = 0.2,
                    segment.size = 0.2,
                    force = 2,
                    size = 6 * 0.35) +
    # geom_text(aes(label = X),
    #           check_overlap = T,
    #           hjust=0, vjust=1,
    #           size = 6 * 0.35) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_x_continuous(expression(mu), limits = c(min(min(data$mu), -max(data$mu)),
                                  max(-min(data$mu), max(data$mu)))) +
    # scale_y_continuous(trans = "log10") +
    scale_y_continuous(expression(sigma)) +
    ggtitle(p) +
    theme_ajf(base_size = 8)
  ggsave(paste("sensitivity_", p, "_mu_vs_sigma.pdf", sep =""), width = 60, height = 60, units = "mm")
}

all_data <- c()
for(p in params){
  data <- read.csv(paste("", p, "_morris.csv", sep = ""))

  data$param <- p
  data <- data %>%
    mutate(scaled_mu_star = .$mu_star / sum(data$mu_star))

  all_data <- rbind(all_data, data)
}

ggplot(all_data) +
  geom_tile(aes(X, param, fill = scaled_mu_star)) +
  scale_fill_viridis_c(expression(paste('normalised ', mu, '*'))) +
  scale_x_discrete("Model parameter") +
  scale_y_discrete("") +
  theme_ajf(base_size = 8) +
  theme(axis.text.x = element_text(angle = 270, vjust = 0, hjust=0))
ggsave("sensitivity_mu-star.pdf", width = 100, height = 60, units = "mm")

