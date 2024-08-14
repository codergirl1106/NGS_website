library(shiny)
library(DT)
library(ggplot2)
library(data.table)
library(bslib)

options(shiny.maxRequestSize=100*1024^2)

shinyUI(fluidPage(
  sidebarLayout(
    sidebarPanel(
      fileInput("r1_file", "upload R1 fastq file"),
      fileInput("r2_file", "upload R2 fastq file"),
      fileInput("metadata_file", "upload PCA primers + barcodes file")
    ),
    mainPanel(
      plotOutput("plot1"),
      headerPanel(""),
      DT::dataTableOutput("sequencetable"),
      headerPanel(""),
      downloadButton("savebutton1", "save DNA sequences table to .csv"),
      headerPanel(""),
      column(width = 12, verbatimTextOutput("savetext1")),
      headerPanel(""),
      plotOutput("plot2"),
      headerPanel(""),
      DT::dataTableOutput("aminoacidtable"),
      headerPanel(""),
      downloadButton("savebutton2", "save amino acid sequences table to .csv"),
      headerPanel(""),
      column(width = 12, verbatimTextOutput("savetext2"))
    )
  )
))
