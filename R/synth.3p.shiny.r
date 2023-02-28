library(shiny)
library(rmeaps)
library(tidyverse)
library(matrixStats)
library(ofce)
library(progressr)
library(rmeaps)
library(patchwork)
library(scales)
library(promises)
library(future)
library(shinyWidgets)
library(shinydashboard)
library(shinybusy)
library(markdown)
library(conflicted)
source("radiation functions.r")
plan("multisession")
conflicts_prefer(shinydashboard::box)
binwidth <- 0.075
beta <- 1.6
rh <- 0.25
re <- rh*sqrt(9/10)
coords <- coord_equal(xlim=c(0,2), ylim=c(0,2))
nshuf <- 256
steps <- 32
nthreads <- 2
coords <- coord_equal(xlim=c(0,2), ylim=c(0,2))

safe_meaps <- purrr::safely(.f = rmeaps::meaps_alt)

options(ofce.base_size = 10, 
        ofce.base_family = "arial", 
        ofce.background_color = "grey99")
set.seed(1942) 

scn_agg <- function(flux, scn) {
  mmb <- meaps_summary(scn$emp, scn$hab, scn$dist, flux)
  mmb$meaps <- NULL
  return(append(mmb, list(egroupes = scn$egroupes, hgroupes = scn$hgroupes)))
}

# Define UI for app that draws a histogram ----
ui <- dashboardPage(
  skin = "black",
  
  # App title ----
  dashboardHeader(title = "MEAPS"),
  # Sidebar layout with input and output definitions ----
  dashboardSidebar(
    width = "250px",
    fluidRow(
      column(6, numericInput(inputId = "n_habitants",
                             label = "habitants",
                             min = 0,
                             max = 5000,
                             step = 100,
                             value = 400)),
      column(6, numericInput(inputId = "k_emplois",
                             label = "emplois",
                             min = 0,
                             max = 5000,
                             step = 100,
                             value = 350))),
    sliderInput(inputId = "pc_habitants",
                label = "part du centre",
                min = 0,
                max = 1,
                step = 0.1,
                value = .7),
    sliderInput(inputId = "distance_c2",
                label = "distance du pôle 2 au pôle 1", 
                min = 0.1,
                value = 0.7,
                step = 0.1,
                max = 1),
    sliderInput(inputId = "angle_c2",
                label = "angle du pôle 2", 
                min = -90,
                value = 45,
                step = 15,
                max = 90),
    sliderInput(inputId = "distance_c3",
                label = "distance du pôle 3 au pôle 1", 
                min = 0.1,
                value = 0.7,
                step = 0.1,
                max = 1),
    sliderInput(inputId = "angle_c3",
                label = "angle du pôle 3", 
                min = -90,
                value = -45,
                step = 15,
                max = 90)),
  ## body ----------
  dashboardBody(
    tabBox(
      width = 12,
      tabPanel("Cartes des habitants et des emplois", 
               fluidRow(
                 column(width=6, plotOutput(outputId = "s1map_h")),
                 column(width=6, plotOutput(outputId = "s1map_e")))),
      tabPanel("Cartes des distances", 
               fluidRow(
                 column(width=6, plotOutput(outputId = "mapdensite_h")),
                 column(width=6, plotOutput(outputId = "mapdensite_e")))),
      tabPanel("Distribution des distances", 
               fluidRow(
                 column(width=6, plotOutput(outputId = "distdensite_h")),
                 column(width=6, plotOutput(outputId = "distdensite_e")))),
      tabPanel("Tension", 
               plotOutput(outputId = "tension"))),
    box(textOutput("jobstatus")),
    box(includeMarkdown("signature.md"))
  ))

# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  
  calculation <- reactiveVal("done")
  step <- reactiveVal(1)
  s1 <- reactive({
    n <- input$n_habitants
    k <- input$k_emplois
    if(n==0) n <- 1
    if(k==0) k <- 1
    f <- if(n>k) (n-k)/n else 0
    pc <- input$pc_habitants
    pc <- min(0.99, max(0.01, pc))
    r2 <- input$distance_c2
    r3 <- input$distance_c3
    a2 <-  input$angle_c2
    a3 <- input$angle_c3
    n1 <- round(pc*n)
    n2 <- round((1-pc)/2*n)
    n3 <- n - n1 - n2
    k1 <- round(pc*k)
    k2 <- round((1-pc)/2*k)
    k3 <- k - k1 - k2
    rh23 <- rh*sqrt((1-pc)/2/pc)
    c1 <- c(1,1)
    c2 <- c1 + c(r2*cos(a2/90*pi/2), r2*sin(a2/90*pi/2))
    c3 <- c1 + c(-r3*cos(-a3/90*pi/2), -r3*sin(-a3/90*pi/2))
    habc <- cbind(pos_cunif(n=n1, rayon = rh, centre = c1, beta=beta), f=f, g = 1)
    habv1 <- cbind(pos_cunif(n=n2, rayon = rh23, centre = c2, beta=beta), f=f, g = 3)
    habv2 <- cbind(pos_cunif(n=n3, rayon = rh23, centre = c3, beta=beta), f=f, g = 2)
    hab <- rbind(habc, habv2, habv1)
    
    empc <- cbind(pos_cunif(n=k1, rayon  = re, centre = c1, beta=beta), p=1, g=1)
    empv1 <- cbind(pos_cunif(n=k2, rayon = re*sqrt((1-pc)/2/pc), centre = c2, beta=beta), p=1, g=3)
    empv2 <- cbind(pos_cunif(n=k3, rayon = re*sqrt((1-pc)/2/pc), centre = c3, beta=beta), p=1, g=2)
    emp <- rbind(empc, empv2, empv1)
    
    shufs <- do.call(rbind, purrr::map(1:nshuf, ~sample.int(n,n)))
    scn <- append(make_tibs(emp, hab, binwidth), list(shufs = shufs))
    
    calculation("pending")
    step(1)
    
    n <- nrow(scn$rk)
    k <- ncol(scn$rk)
    un_shuf <- (scn$shufs)[1,, drop=FALSE]
    
    raw <- meaps_alt(
      rkdist=scn$rk, 
      emplois = rep(1,k), 
      actifs = rep(1,n),
      modds = matrix(1, ncol=k, nrow=n),
      f = scn$f,
      shuf = un_shuf,
      nthreads = 1,
      progress = FALSE
    )
    acc_flux(raw)
    agg_flux(scn_agg(raw, scn))
    plan("multisession")
    
    return(scn)
  })
  
  acc_flux <- reactiveVal(NULL)
  agg_flux <- reactiveVal(NULL)
  fluxs1HD <- reactive({
    scn <- s1()
    n <- nrow(scn$rk)
    k <- ncol(scn$rk)
    le_step <- step()
    if(le_step > nshuf %/% steps)
      return()
    les_ranks <- scn$rk
    la_fuite <- scn$f
    les_shufs <- (scn$shufs)[((le_step-1)*steps+1):(le_step*steps),]
    future({
      rmeaps::meaps_alt(
        rkdist = les_ranks, 
        emplois = rep(1,k), 
        actifs = rep(1,n),
        modds = matrix(1, ncol=k, nrow=n),
        f = la_fuite,
        shuf = les_shufs,
        nthreads = 4,
        progress = FALSE)
    },
    seed=TRUE, 
    globals = c("les_ranks", "la_fuite", "les_shufs", "n", "k"),
    packages = c("rmeaps"))
  })
  
  observe({
    le_step <- step()
    if(resolved(fluxs1HD())) {
      if(le_step == nshuf %/% steps) 
        calculation("done")
      else {
        scn <- isolate(s1())
        if(le_step==1) {
          new_flux <- value(isolate(fluxs1HD()))
          acc_flux(new_flux)
        } else {
          past_flux <- isolate(acc_flux())
          new_flux <- value(isolate(fluxs1HD()))
          new_flux <- (le_step * past_flux + new_flux)/(le_step+1)
          acc_flux(new_flux)
        }
        agg_flux(scn_agg(new_flux, scn))
        step(le_step+1)
      }
    }
    if(calculation()!="done")
      invalidateLater(500, session)
  })
  ## job status --------------
  output$jobstatus <- renderText({
    le_step <- step()
    if(calculation()!="done")
      glue::glue("Montecarlo sur {steps * le_step}/{nshuf} tirages, en cours")
    else
      glue::glue("Montecarlo sur {nshuf} tirages, terminé")
  })
  # carte du scénario ----------
  output$s1map_h <- renderPlot({
    scn <- s1()
    ggplot()+
      stat_binhex(data = scn$habs,
                  aes(x=x,y=y, fill=100*after_stat(density)), binwidth=binwidth)+
      scale_fill_distiller(palette="Greens", direction=1, name = "densité\nd'habitants")+
      coords +
      geom_text(data = scn$hgroupes, aes(x=x, y=y, label = g_label), nudge_y = 0.3,  size = 4) +
      labs(title = "Habitants")+
      theme_ofce_void()+ 
      theme(plot.margin = margin(6,6,6,6),
            panel.background = element_rect(fill="grey99"))})
  
  output$s1map_e <- renderPlot({
    scn <- s1()
    ggplot()+
      stat_binhex(data=scn$emps,
                  aes(x=x,y=y, fill=100*after_stat(density)), binwidth=binwidth)+
      scale_fill_distiller(palette = "Oranges", direction=1, name = "densité\nd'emplois")+
      coords +
      geom_text(data = scn$egroupes, aes(x=x, y=y, label = g_label), size = 4, nudge_y = 0.2) +
      labs(title = "Emplois")+
      theme_ofce_void()+ 
      theme(plot.margin = margin(6,6,6,6)) })
  ## carte densité --------------
  output$mapdensite_h <- renderPlot({
    mmb <- agg_flux()
    ggplot()+
      stat_summary_hex(data = mmb$hab, aes(x=x, y=y, z=d), binwidth=binwidth)+
      guides(fill=guide_colourbar("Distance\npar habitant"))+
      scale_fill_distiller(palette="Greens", direction=1)+
      xlab(NULL)+ylab(NULL)+labs(title="Distance moyenne parcourue par habitant")+
      geom_text(data = mmb$hgroupes, aes(x=x, y=y, label = g_label), nudge_y = -0.2,  size = 4)+
      coords +
      theme_ofce_void() +
      theme(legend.position="right", 
            plot.margin = margin(6,6,6,6))
  })
  output$mapdensite_e <- renderPlot({
    mmb <- agg_flux()
    ggplot()+
      stat_summary_hex(data = mmb$emp, aes(x=x, y=y, z=d), binwidth=binwidth)+
      guides(fill=guide_colourbar("Distance\npour un emploi"))+
      scale_fill_distiller(palette="Oranges", direction=1)+
      xlab(NULL)+ylab(NULL)+labs(title="Distance moyenne parcourue pour un emploi") +
      geom_text(data = mmb$egroupes, aes(x=x, y=y, label = g_label), size = 4, nudge_y = -0.2) +
      coords +
      theme_ofce_void()+
      theme(legend.position="right",
            plot.margin = margin(6,6,6,6))
  })
  ## densite des distances ---------
  output$distdensite_h <- renderPlot({
    le_step <- step()
    if(le_step<nshuf %/% steps) {
      alpha <- 0.75+0.5*(le_step-nshuf %/% steps)/(nshuf %/% steps -1)
      line <- "dashed"} 
    else {
      alpha <- 0.75
      line <- "solid"
    }
    mmb <- agg_flux()
    ggplot(mmb$hab)+
      geom_density(aes(x=d, fill=factor(g), col=factor(g), 
                       y=after_stat(count)/nrow(mmb$hab)), 
                   position="stack", alpha = alpha, linetype = line)+
      scale_color_brewer(
        palette ="Set2",
        name="pôle d'habitation",
        aesthetics = c('color','fill'), 
        labels = c("h1", "h2", "h3"))+
      xlim(c(0,1))+
      scale_y_continuous(limits=c(0,10), oob = squish)+
      ylab(NULL)+xlab("distance")+
      labs(title="Distances parcourue par habitant") +
      theme_ofce(legend.position = "right")
    })
  
  output$distdensite_e <- renderPlot({
    le_step <- step()
    if(le_step<nshuf %/% steps) {
      alpha <- 0.75+0.5*(le_step-nshuf %/% steps)/(nshuf %/% steps -1)
      line <- "dashed"} 
    else {
      alpha <- 0.75
      line <- "solid"
    }
    mmb <- agg_flux()
    ggplot(mmb$emp)+
      geom_density(aes(x=d, fill=factor(g), col=factor(g), y=after_stat(count)/nrow(mmb$emp)), 
                   position="stack", alpha = alpha, linetype = line)+
      scale_color_brewer(
        palette ="Set3",
        name="pôle d'emploi",
        aesthetics = c('color','fill'), 
        labels = c("e1", "e2", "e3"))+
      xlim(c(0,1))+scale_y_continuous(limits=c(0,10), oob = squish)+
      ylab(NULL)+xlab("distance")+
      labs(title="Distances parcourues pour un emploi") +
      theme_ofce(legend.position = "right")
  })
  output$tension <- renderPlot({
    ggplot()
  })
}

shinyApp(ui = ui, server = server)
