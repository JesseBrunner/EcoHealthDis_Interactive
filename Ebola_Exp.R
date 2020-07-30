library(tidyverse)


# Congo data set
Congo <- read_csv("Ebola_Congo_1995.csv"
                  )
Congo <- Congo %>%
  mutate(Cumulative = cumsum(Cases)) %>%
  filter(Day >=7)


# Uganda data set
Uganda <- read_csv("Ebola_Uganda_2000.csv")
Uganda <- Uganda %>%
  mutate(Cumulative = cumsum(Cases))

# Plot data

df <- Congo # choose data set
ylog <- TRUE # choose log-y-axis

cutoff <- 68 # choose cutoff
# True cutoff was 63 for Congo and 56 for Uganda


P <- ggplot(df, aes(x=Day, y=Cumulative)) +
  geom_point() +
  geom_smooth(method="lm",
              data=filter(df, Day <= cutoff),
              fullrange=TRUE, se=F, linetype=2) +
  geom_smooth(method="lm",
              data=filter(df, Day <= cutoff)) +
  scale_x_continuous("Days since first case", breaks = 10*0:16) +

  coord_cartesian(ylim=c(1, max(df$Cumulative)*2))

ifelse(ylog,
       plot(P + scale_y_log10("Cumulative number of cases")),
       plot(P + scale_y_continuous("Cumulative number of cases"))
)

summary(lm(log(Cumulative) ~ Day, data=filter(df, Day <= cutoff)))
