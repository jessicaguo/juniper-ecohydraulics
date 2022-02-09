#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(tidyverse)
library(shiny)
library(ggplot2)
load("met_daily.Rdata")
load("psy_all.Rdata")
load("maint.Rdata")


# Define UI for application with two tabs
ui <- navbarPage("Daily psychrometer data",
  tabPanel("VPD",
           titlePanel("Psychrometer vs. VPD"),
           # Sidebar with drop down input for House number
           sidebarLayout(
             sidebarPanel(
               # Select Tree
               selectInput(inputId = "Tree1",
                           label = "Select Tree",
                           choices = levels(psy_all$Tree)),
               # Select Logger
               selectInput(inputId = "Logger1",
                           label = "Select Logger",
                           choices = levels(psy_all$Logger)),
               # Select range of dates
               sliderInput("slider1", label = h3("Date Range"), 
                           min = min(psy_all$date), 
                           max = max(psy_all$date), 
                           value = range(psy_all$date)
               ),
             ),
             # Show a size plot for selected species
             mainPanel(
               fluidRow(plotOutput("psyVPD_ts", width = "100%", height = "250px")),
               fluidRow(plotOutput("psyVPD_scatter", width = "800px", height = "400px"))
             )
           )),
  tabPanel("VWC",
           titlePanel("Psychrometer vs. VWC"),
           sidebarLayout(
             sidebarPanel(
               # Select Tree
               selectInput(inputId = "Tree2",
                           label = "Select Tree",
                           choices = unique(psy_all$Tree)),
               # Select Logger
               selectInput(inputId = "Logger2",
                           label = "Select Logger",
                           choices = levels(psy_all$Logger)),
               # Select range of dates
               sliderInput("slider2", label = h3("Date Range"), 
                           min = min(psy_all$date), 
                           max = max(psy_all$date), 
                           value = range(psy_all$date)
               ),
             ),
             
             # Show a size plot for selected species
             mainPanel(
               fluidRow(plotOutput("psyVWC_ts", width = "100%", height = "250px")),
               fluidRow(plotOutput("psyVWC_scatter", width = "800px", height = "800px"))
             )
           ))
  
)


server <- function(input, output) {
  # Function to plot psy and VPD time series
  output$psyVPD_ts <- renderPlot({
    
    maint_temp <- maint %>%
      filter(st_date >= input$slider1[1],
             en_date <= input$slider1[2])
    
    psy_temp <- psy_all %>%
      filter(Tree == input$Tree1,
             Logger == input$Logger1,
             date >= input$slider1[1],
             date <= input$slider1[2])
    
    met_temp <- met_daily %>%
      filter(date >= input$slider1[1],
             date <= input$slider1[2])
    
    ggplot() +
      geom_rect(data = maint_temp, 
                aes(xmin = st_date, xmax = en_date, 
                    ymin= -Inf, ymax = Inf),
                alpha = 0.5) +
      geom_point(data = psy_temp,
                 aes(x = date,y = MD, col = "MD (MPa)"),
                 size = 2) +
      geom_point(data = psy_temp,
                 aes(x = date,y = PD, col = "PD (MPa)"),
                 size = 2) +
      geom_line(data = met_temp,
                aes(x = date,y = VPD_max, col = "VPD (kPa)"),
                size = 1.5) +
      geom_bar(data = met_temp,
               aes(x = date,y = Precip/10, fill = "Precip"),
               stat = "identity") +
      scale_x_date(date_labels = "%b",
                   date_breaks = "1 month") +
      scale_y_continuous(expression(paste(Psi, " | VPD | Precip (cm)"))) +
      theme_bw(base_size = 16) +
      theme(strip.background = element_blank(),
            axis.title.x = element_blank(),
            legend.title = element_blank()) +
      scale_color_manual(values = c("darkkhaki",
                                    "aquamarine3",
                                    "coral")) +
      scale_fill_manual(values = "cornflowerblue") +
      guides(fill = "none",
             color = guide_legend(
               override.aes = list(linetype = c(0, 0, 1),
                                   shape = c(16, 16, NA),
                                   size = c(3, 3, 1.5))))
  })
  
  # Function to plot VPD vs. MD or PD relationship
  output$psyVPD_scatter <- renderPlot({
    psy_all %>%
      pivot_longer(cols = MD:PD,
                   names_to = "type",
                   values_to = "WP") %>%
      filter(Tree == input$Tree1,
             Logger == input$Logger1,
             date >= input$slider1[1],
             date <= input$slider1[2]) %>%
      ggplot(mapping = aes(x = VPD_max, y = WP, color = type)) +
      geom_point(size = 2) +
      geom_smooth(method = lm, 
                  se = FALSE,
                  fullrange = TRUE,
                  lty = 2,
                  lwd = 0.75) +
      facet_wrap(~type, ncol = 2) +
      scale_x_continuous(expression(paste(VPD[max], " (kPa)"))) +
      scale_y_continuous(expression(paste(Psi, " (MPa)"))) +
      theme_bw(base_size = 16) +
      scale_color_manual(values = c("darkkhaki",
                                      "aquamarine3")) +
      theme(strip.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      guides(color = "none")
  })

  # Function to plot psy and VWC time series
  output$psyVWC_ts <- renderPlot({
    
    maint_temp <- maint %>%
      filter(st_date >= input$slider2[1],
             en_date <= input$slider2[2])
    
    psy_temp <- psy_all %>%
      filter(Tree == input$Tree2,
             Logger == input$Logger2,
             date >= input$slider2[1],
             date <= input$slider2[2])
    
    met_temp <- met_daily %>%
      filter(date >= input$slider2[1],
             date <= input$slider2[2])
    
    ggplot() +
      geom_rect(data = maint_temp, 
                aes(xmin = st_date, xmax = en_date, 
                    ymin= -Inf, ymax = Inf),
                alpha = 0.5) +
      geom_point(data = psy_temp,
                 aes(x = date,y = MD, col = "MD"),
                 size = 2) +
      geom_point(data = psy_temp,
                 aes(x = date,y = PD, col = "PD"),
                 size = 2) +
      geom_bar(data = met_temp,
               aes(x = date,y = Precip/10, fill = "Precip"),
               stat = "identity") +
      geom_line(data = met_temp,
                aes(x = date,y = VWC_5cm/2, col = "VWC (5 cm)"),
                size = 1.5) +
      geom_line(data = met_temp,
                aes(x = date,y = VWC_10cm/2, col = "VWC (10 cm)"),
                size = 1.5) +
      geom_line(data = met_temp,
                aes(x = date,y = VWC_20cm/2, col = "VWC (20 cm)"),
                size = 1.5) +
      scale_x_date(date_labels = "%b",
                   date_breaks = "1 month") +
      scale_y_continuous(expression(paste(Psi, " (MPa) | Precip (cm)")),
                         sec.axis = sec_axis(~.*2, 
                                             name = "VWC (%)")) +
      theme_bw(base_size = 16) +
      theme(axis.title.x = element_blank(),
            legend.title = element_blank()) +
      scale_color_manual(values = c("aquamarine3",
                                    "darkkhaki",
                                    "mistyrose4",
                                    "darkgoldenrod",
                                    "chocolate4"),
                         breaks = c("PD", "MD",
                                    "VWC (5 cm)",
                                    "VWC (10 cm)",
                                    "VWC (20 cm)")) +
      scale_fill_manual(values = "cornflowerblue") +
      guides(fill = "none",
             color = guide_legend(
               override.aes = list(linetype = c(0, 0, 1, 1, 1),
                                   shape = c(16, 16, NA, NA, NA),
                                   size = c(3, 3, 1.5, 1.5, 1.5))))
  })
  
  # Function to plot VWC vs. MD or PD relationship
  output$psyVWC_scatter <- renderPlot({
    psy_all %>%
      pivot_longer(cols = MD:PD,
                   names_to = "type",
                   values_to = "WP") %>%
      pivot_longer(cols = VWC_5cm:VWC_100cm,
                   names_to = "depth",
                   names_prefix = "VWC_",
                   values_to = "VWC") %>%
      filter(depth %in% c("5cm", "10cm", "20cm")) %>%
      mutate(depth = factor(depth, levels = c("5cm", "10cm", "20cm"))) %>%
      filter(Tree == input$Tree2,
             Logger == input$Logger2,
             date >= input$slider2[1],
             date <= input$slider2[2]) %>%
      ggplot(mapping = aes(x = VWC, y = WP, color = type)) +
      geom_point(size = 2) +
      geom_smooth(method = lm, 
                  se = FALSE,
                  fullrange = TRUE,
                  lty = 2,
                  lwd = 0.75) +
      facet_grid(rows = vars(depth),
                 cols = vars(type),
                 scales = "free") +
      scale_x_continuous(expression(paste(VWC, " (%)"))) +
      scale_y_continuous(expression(paste(Psi, " (MPa)"))) +
      theme_bw(base_size = 16) +
      scale_color_manual(values = c("darkkhaki",
                                    "aquamarine3")) +
      theme(strip.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      guides(color = "none")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
