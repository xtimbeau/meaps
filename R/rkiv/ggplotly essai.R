{r}
#| label: fig-spectreRi
#| fig-scap: "Odd-ratio par commune de résidence fonction de la distance aux communes d'emploi (spectre résidents)"
#| fig-cap: "La figure représente pour les 20 plus grandes communes de l'agglomération de la Rochelle les odd-ratios estimés (configuration 100% des flux) en fonction de la distance entre cette commune et les communes où travaillent les résidents. Les points marqués d'un petit point blancs sont les emplois situés hors du périmètre du SCoT."
library(plotly)
load("output/odds_eff.srda")
data <- odds_eff |> 
  filter(emobpro>450, d<40000) |> 
  select(d, n=nDCLT, o=or_eff, c=config, m=mobpro) |> 
  mutate(c =str_c("c", c |> str_replace_all("é", "e") |> str_replace_all("%",""))) |> 
  pivot_wider(id_cols = c(d,n,m), values_from = o, names_from=c)
gspctrE <- ggplot(data)+
  geom_point(aes(x=d/1000, y = c99, col=m))+
  scale_y_log10(limits=c(0.05, 50))+
  facet_wrap(vars(n))+ 
  theme_ofce(base_size = 4, base_family = "Nunito")
ggplotly(gspctrE) |> 
  layout(
    updatemenus = list(
      list(
        buttons = list(
          list(method = "restyle",
               args = list("y", list(data$diagonale)),  # put it in a list
               label = "Diagonale"),
          list(method = "restyle",
               args = list("y", list(data$c99)),  # put it in a list
               label = "99%")))))