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
library(fresh)
library(shinydashboard)
library(shinybusy)
library(markdown)
library(conflicted)
library(fontawesome)
source("radiation functions.r")
plan("multisession")
conflicts_prefer(shinydashboard::box)
conflicts_prefer(dplyr::filter)
binwidth <- 0.075
beta <- 1.6
rh <- 0.25
re <- rh*sqrt(9/10)
nshuf <- 256
steps <- 16
nthreads <- 4
coords <- coord_equal(xlim=c(-1,1), ylim=c(-1,1))

safe_meaps <- purrr::safely(.f = rmeaps::meaps_alt)
safe_resolved <- purrr::safely(.f = future::resolved)

options(ofce.base_size = 10, 
        ofce.base_family = "arial", 
        ofce.background_color = "grey99")
set.seed(1942) 

scn_agg <- function(flux, scn) {
  mmb <- meaps_summary(scn$emp, scn$hab, scn$dist, flux)
  mmb$meaps <- NULL
  return(append(mmb, list(egroupes = scn$egroupes, hgroupes = scn$hgroupes)))
}

# fresh theme -------------
mytheme <- create_theme(
  adminlte_color(
    light_blue = "#434C5E"
  ),
  adminlte_sidebar(
    width = "300px",
    dark_bg = "#efefef",
    dark_hover_bg = "#81A1C1",
    dark_color = "#2E3440",
    dark_submenu_color = "#111"
  ),
  adminlte_global(
    content_bg = "#FFF",
    box_bg = "#efefef", 
    info_box_bg = "#efefef"
  ),
  theme="flatly"
)

# Define UI for app that draws a histogram ----
ui <- dashboardPage(
  skin = "black",
  # App title ----
  dashboardHeader(title = "MEAPS"),
  # Sidebar layout with input and output definitions ----
  dashboardSidebar(
    fluidRow(
      column(6, numericInput(inputId = "n_habitants",
                             label = list(fa("people-group", fill="green" ), "actifs"),
                             min = 100,
                             max = 5000,
                             step = 100,
                             value = 1000)),
      column(6, numericInput(inputId = "k_emplois",
                             label = list(fa("person-digging", fill="orange"), "emplois"),
                             min = 100,
                             max = 5000,
                             step = 100,
                             value = 800))),
    fluidRow(
      column(6, numericInput(inputId = "fuite",
                             label = list(icon("person-running"), "fuite"),
                             min = .01,
                             max = .2,
                             step = .01,
                             value = .05)),
      column(6, numericInput(inputId = "beta",
                             label = list(icon("bezier-curve"), "beta"),
                             min = 0.5,
                             max = 2,
                             step = .05,
                             value = 1.5))),
    fluidRow(
      column(6, numericInput(inputId = "rayon",
                             label = list(icon("ruler-horizontal"),"rayon du centre"),
                             min = 0.1,
                             max = 1.5,
                             step = .1,
                             value = .5))),
    fluidRow(
      column(6, numericInput(inputId = "pc_habitants",
                             label = "% actifs au centre",
                             min = 0,
                             max = 100,
                             step = 5,
                             value = 70)),
      column(6, numericInput(inputId = "pc_emplois",
                             label = "% emplois au centre",
                             min = 0,
                             max = 100,
                             step = 5,
                             value = 70))),
    fluidRow(
      column(6, numericInput(inputId = "distance_c2",
                             label = list(icon("ruler-horizontal"),"dist. p2-centre"),
                             min = 0,
                             max = 1.5,
                             step = .05,
                             value = .75)),
      column(6, numericInput(inputId = "distance_c3",
                             label = list(icon("ruler-horizontal"),"dist. p3-centre"),
                             min = 0,
                             max = 1.5,
                             step = .05,
                             value = .75))),
    fluidRow(
      column(6, numericInput(inputId = "angle_c2",
                             label = list(icon("angle-down"), "angle du p2 (°)"),
                             min = -90,
                             value = 45,
                             step = 15,
                             max = 90)),
      column(6, numericInput(inputId = "angle_c3",
                             label = list(icon("angle-up"), "angle du p3 (°)"),
                             min = -90,
                             value = -45,
                             step = 15,
                             max = 90)))),
  # Body ----------
  dashboardBody(
    use_theme(mytheme),
    chooseSliderSkin(skin = "Modern"),
    tags$head(
      tags$style(HTML("
                      .sidebar {
                        color: #333333;
                        font-weight: normal;}
                      .irs-grid-text {
                        color: #555555;
                      }"))),
    tabBox(
      width = 12,
      tabPanel("Cartes des actifs et des emplois", 
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
      tabPanel("Tension localisée", 
               plotOutput(outputId = "tension"))),
    box(
      progressBar("job_pb", value=0, display_pct = TRUE),
      textOutput("jobstatus")),
    box(includeMarkdown("signature.md"))
  ))

# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  
  calculation <- reactiveVal("done")
  step <- reactiveVal(1)
  s1 <- reactive({
    scn <- genere_3p(
      n = input$n_habitants,
      k = input$k_emplois,
      f = input$fuite,
      part_h = input$pc_habitants/100,
      part_e = input$pc_emplois/100,
      d_cp2 = input$distance_c2,
      d_cp3 = input$distance_c3,
      theta2 = input$angle_c2,
      theta3 = input$angle_c3,
      beta = input$beta,
      rayon = input$rayon,
      nshuf = nshuf,
      binwidth = binwidth) |> 
      add_dist()
    step(0)
    return(scn)
  }) |> debounce(millis = 1000)
  
  step1 <- observe({
    scn <- s1()
    un_shuf <- (scn$shufs)[1:steps,]
    raw <- rmeaps::meaps_tension_alt(
      rkdist=scn$rk, 
      emplois = rep(1,scn$k), 
      actifs = rep(1,scn$n),
      modds = matrix(1, ncol=scn$k, nrow=scn$n),
      f = scn$f,
      shuf = un_shuf,
      nthreads = 4,
      progress = FALSE
    )
    acc_flux(raw$flux)
    acc_tension(list(tension = raw$tension,
                     n = scn$n, k = scn$k,
                     emps = scn$emps, egroupes = scn$egroupes))
    agg_flux(scn_agg(raw$flux, scn))
    calculation("pending")
    step(1)
    updateProgressBar(id="job_pb", value = 1, total = nshuf %/% steps)
    plan("multisession")
  })
  
  acc_flux <- reactiveVal(NULL)
  acc_tension <- reactiveVal(NULL)
  agg_flux <- reactiveVal(NULL)
  
  fluxs1HD <- reactive(
    {
      le_step <- step()
      if(le_step==0|le_step > nshuf %/% steps)
        return()
      scn <- s1()
      n <- nrow(scn$rk)
      k <- ncol(scn$rk)
      les_ranks <- scn$rk
      la_fuite <- scn$f
      les_shufs <- (scn$shufs)[((le_step-1)*steps+1):(le_step*steps),]
      future({
        rmeaps::meaps_tension_alt(
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
  
  observe(
    {
      le_step <- step()
      if(le_step==0)
        return()
      state <- safe_resolved(fluxs1HD())
      if(!is.null(state$error))
        return()
      if(state$result==TRUE) {
        if(le_step == nshuf %/% steps) {
          updateProgressBar(id="job_pb", value = le_step, total = nshuf %/% steps)
          calculation("done")
        }
        else {
          scn <- isolate(s1())
          past_flux <- isolate(acc_flux())
          past_tension <- isolate(acc_tension())$tension
          new_data <- value(isolate(fluxs1HD()))
          new_flux <- (le_step * past_flux + new_data$flux)/(le_step+1)
          new_tension <- (le_step * past_tension + new_data$tension)/(le_step+1)
          acc_flux(new_flux)
          acc_tension(list(
            tension = new_tension,
            n = scn$n, k = scn$k,
            emps = scn$emps, egroupes = scn$egroupes))
          agg_flux(scn_agg(new_flux, scn))
          updateProgressBar(id="job_pb", value = le_step, total = nshuf %/% steps)
          step(le_step+1)
        }
      }
      if(calculation()!="done")
        invalidateLater(250, session)
    })
  
  ## job status --------------
  output$jobstatus <- renderText({
    le_step <- step()
    if(calculation()!="done")
      glue::glue("Monte-Carlo sur {steps * (le_step-1)}/{nshuf} tirages, en cours")
    else
      glue::glue("Monte-Carlo sur {nshuf} tirages, terminé")
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
    at <- acc_tension()
    rangns_max <- at$emps |> 
      mutate(tension = at$tension) |> 
      mutate(tension  = (at$n - tension)/(at$n))
    ggplot(rangns_max)+
      stat_summary_hex(aes(x=x,y=y, z=tension), binwidth=binwidth)+
      scale_fill_distiller(palette="Spectral", direction=-1, name = "Tension", 
                           limits = c(0, .5), breaks = c(0, .1, .2, .3, .4, .5),
                           oob =squish)+
      coords +
      geom_text(data = at$egroupes, aes(x=x, y=y, label = g_label), nudge_y = 0.3,  size = 2) +
      labs(title = "Tension")+
      theme_ofce_void(base_size = 8)+ 
      theme(plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b=4)),
            plot.margin = margin(6,6,6,6),
            panel.background = element_rect(fill="grey97"))
  })
}

shinyApp(ui = ui, server = server)
