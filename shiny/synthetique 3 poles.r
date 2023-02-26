library(shiny)
library(rmeaps)
library(tidyverse)
library(matrixStats)
library(ofce)
library(progressr)
library(rmeaps)
library(patchwork)
source("../R/radiation functions.r")
binwidth <- 0.075
beta <- 1.6
rh <- 0.25
re <- rh*sqrt(9/10)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("MEAPS simulations synthétiques"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "n_habitants",
                  label = "Nombre d'habitants",
                  min = 1,
                  max = 5000,
                  value = 1000),
      
      sliderInput(inputId = "k_emplois",
                  label = "Nombre d'emplois",
                  min = 1,
                  max = 5000,
                  value = 1000),
      
      sliderInput(inputId = "pc_habitants",
                  label = "Part du centre (habitants)",
                  min = 0,
                  max = 1,
                  value = .7)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotOutput(outputId = "s1map")
      
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  s1 <- reactive({
    n <- input$n_habitants
    k <- input$k_emplois
    pc <- input$pc_habitants
    habc <- cbind(pos_cunif(n=pc*n, rayon = rh, centre = c(1, 1), beta=beta), f=0.1, g = 1)
    habv1 <- cbind(pos_cunif(n=(1-pc)/2*n, rayon = rh*sqrt((1-pc)/2/pc), centre = c(.5, .5), beta=beta), f=0.1, g = 3)
    habv2 <- cbind(pos_cunif(n=(1-pc)/2*n, rayon = rh*sqrt((1-pc)/2/pc), centre = c(1.5, 1.5), beta=beta), f=0.1, g = 2)
    hab <- rbind(habc, habv2, habv1)
    
    empc <- cbind(pos_cunif(n=pc*k, rayon  = re, centre = c(1, 1), beta=beta), p=1, g=1)
    empv1 <- cbind(pos_cunif(n=(1-pc)/2*k, rayon = re*sqrt((1-pc)/2/pc), centre = c(.5, .5), beta=beta), p=1, g=3)
    empv2 <- cbind(pos_cunif(n=(1-pc)/2*k, rayon = re*sqrt((1-pc)/2/pc), centre = c(1.5, 1.5), beta=beta), p=1, g=2)
    emp <- rbind(empc, empv2, empv1)
    
    make_tibs(emp, hab, binwidth)
  })
  
  s2 <- reactive({
     habc <- cbind(pos_cunif(n=pc*n, rayon = rh, centre = c(1, 1), beta=beta), f=0.1, g = 1)
    habv1 <- cbind(pos_cunif(n=(1-pc)/2*n, rayon = rh*sqrt((1-pc)/2/pc), centre = c(.5, .5), beta=beta), f=0.1, g = 3)
    habv12 <- t(t(habv1) + c(-.25, -.25, 0, 0))
    habv2 <- cbind(pos_cunif(n=(1-pc)/2*n, rayon = rh*sqrt((1-pc)/2/pc), centre = c(1.5, 1.5), beta=beta), f=0.1, g = 2)
    habv22 <- t(t(habv2) + c(.25, .25, 0, 0))
    hab <- rbind(habc, habv2, habv1)
    hab2 <- rbind(habc, habv22, habv12)
    empc <- cbind(pos_cunif(n=pc*k, rayon  = re, centre = c(1, 1), beta=beta), p=1, g=1)
    empv1 <- cbind(pos_cunif(n=(1-pc)/2*k, rayon = re*sqrt((1-pc)/2/pc), centre = c(.5, .5), beta=beta), p=1, g=3)
    empv12 <- t(t(empv1) + c(-.25, -.25, 0, 0))
    empv2 <- cbind(pos_cunif(n=(1-pc)/2*k, rayon = re*sqrt((1-pc)/2/pc), centre = c(1.5, 1.5), beta=beta), p=1, g=2)
    empv22 <- t(t(empv2) + c(.25, .25, 0, 0))
    emp <- rbind(empc, empv2, empv1)
    emp2<- rbind(empc, empv22, empv12)
    
    make_tibs(emp, hab, binwidth)
  })
  
  output$s1map <- renderPlot({
    coords <- coord_equal(xlim=c(0,2), ylim=c(0,2))
    scn <- s1()
    (ggplot()+
           stat_binhex(data = scn$habs,
                       aes(x=x,y=y, fill=100*after_stat(density)), binwidth=binwidth)+
           scale_fill_distiller(palette="Greens", direction=1, name = "densité\nd'habitants",
                                breaks = c(1,2))+
           coords +
           geom_text(data = scn$hgroupes, aes(x=x, y=y, label = g_label), nudge_y = 0.3,  size = 2) +
           labs(title = "Habitants")+
           theme_ofce_void(base_size = 9, base_family = "sans")+ 
           theme(plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b=4)),
                 plot.margin = margin(6,6,6,6),
                 panel.background = element_rect(fill="grey97")))+
        (ggplot()+
           stat_binhex(data=scn$emps,
                       aes(x=x,y=y, fill=100*after_stat(density)), binwidth=binwidth)+
           scale_fill_distiller(palette = "Oranges", direction=1, name = "densité\nd'emplois", 
                                breaks = c(1,2))+
           coords +
           geom_text(data = scn$egroupes, aes(x=x, y=y, label = g_label), size = 2, nudge_y = 0.2) +
           labs(title = "Emplois")+
           theme_ofce_void(base_size = 9, base_family="sans")+ 
           theme(plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b=4)),
                 plot.margin = margin(6,6,6,6),
                 panel.background = element_rect(fill="grey97"))) + 
        plot_layout(guides = 'collect')
    
  })
  
}

shinyApp(ui = ui, server = server)