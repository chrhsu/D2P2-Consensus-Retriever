# Remember to change the working directory or this will not work

path = file.choose()



library(ggplot2)
df <- read.csv(path)

df2 = data.frame(df["Organism"], df["Mean.of.Sample.Mean"], df["Core.or.NonCore"], df["Standard.Deviation"])

df2$Organism <- factor(df2$Organism, levels = c("Human", "Cow", "Mouse", "Zebrafish", "Fly", "Worm", "Yeast", "Frog"))

ggplot(data=df2, aes(x=Organism, y=Mean.of.Sample.Mean, fill=Core.or.NonCore)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Mean.of.Sample.Mean-Standard.Deviation, ymax=Mean.of.Sample.Mean+Standard.Deviation), width=.2,position=position_dodge(.9)) +
  ggtitle("Nucleolus Core and NonCore Disorder Protein") +
  theme(plot.title = element_text(hjust = 0.5)) + labs(x="Organism", y="Mean of Sample Means (of Disorder)", fill="Type") + ggsave("hq_map.png", dpi = 300)

