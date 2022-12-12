library(ggplot2)
# Create a scatter plot
a <- ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width)) +
  geom_point(aes(color = Species)) +
  scale_color_viridis_d() +
  theme_minimal()

p <- ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width))+
  geom_point(aes(color = Species, shape = Species))+
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))


multi_panel <- ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width))+
  geom_point(aes(color = Species))+
  geom_smooth(aes(color = Species, fill = Species))+
  facet_wrap(~Species, ncol = 3, nrow = 1)+
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
print(multi_panel)