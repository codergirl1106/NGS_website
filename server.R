library(shiny)
library(Biostrings)
library(Rfastp)
library(DT)
library(ggplot2)
library(data.table)
library(stringr)
library(tibble)
library(dplyr)
library(bslib)

merging <- function(file1, file2, out) {
  o <- rfastp(read1 = file1, read2 = file2,
              merge = TRUE, outputFastq = "", mergeOut = out)
}

convert_StringSet_to_df = function(stringset, combine_duplicates = T) {
  sequence_df = as_tibble(as.character(stringset))
  names(sequence_df) = "sequence"
  if (combine_duplicates) {
    sequence_df = sequence_df %>%
      group_by(sequence) %>%
      summarise(read_count = n()) %>%
      ungroup() %>%
      arrange(desc(read_count), sequence)
  } else {
    sequence_df = sequence_df %>%
      mutate(read_count = 1)
  }
  return(as.data.frame(sequence_df))
}

translate_reads = function(reads) {
  translations = translate(reads, if.fuzzy.codon = "solve")
  return(as.character(translations))
}

shinyServer(function(input, output) {
  set.seed(1234)
  merged <- reactive({
    req(input$r1_file)
    req(input$r2_file)
    req(input$metadata_file)
    
    r1_file <- file.path(getwd(), "r1_file.fastq")
    
    writeLines(readLines(input$r1_file$datapath), con = r1_file)
    
    r2_file <- file.path(getwd(), "r2_file.fastq")
    
    writeLines(readLines(input$r2_file$datapath), con = r2_file)
    
    metadata_file <- file.path(getwd(), "metadata_file.csv")
    
    writeLines(readLines(input$metadata_file$datapath), con = metadata_file)
    
    temp_output <- file.path(getwd(), "temp_fastq_output.fastq.gz")
    
    merging(r1_file, r2_file, temp_output)
    
    fasta <- readLines(temp_output)
    
    file.remove(temp_output)
    file.remove(r1_file)
    file.remove(r2_file)
    
    sequences <- fasta[seq(2, length(fasta), by = 4)]
    sequences <- sequences[sequences != ""]
    
    primers <- fread(metadata_file, header=TRUE)
    
    file.remove(metadata_file)
    
    out <- NULL
    for (row in 1:nrow(primers)) {
      F_seq <- as.character(primers[row,"Fwd Seq"])
      R_seq <- as.character(primers[row,"Rev Seq"])
      
      n1 <- unlist(str_extract_all(sequences, paste0("^", F_seq, "(.*)", as.vector(reverseComplement(DNAStringSet(R_seq)))[1], "$")))
      n2 <- unlist(str_extract_all(sequences, paste0("^", R_seq, "(.*)", as.vector(reverseComplement(DNAStringSet(F_seq)))[1], "$")))
      
      n1 <- substring(n1, nchar(F_seq) + 1, nchar(n1) - nchar(R_seq))
      n2 <- substring(n2, nchar(R_seq) + 1, nchar(n2) - nchar(F_seq))
      
      n1 <- n1[n1 != ""] %>% convert_StringSet_to_df()
      n2 <- n2[n2 != ""] %>% convert_StringSet_to_df()
      
      colnames(n1)[colnames(n1) == "read_count"] <- primers[row,]$Primer
      colnames(n2)[colnames(n2) == "read_count"] <- primers[row,]$Primer
      
      comb <- rbind(n1, n2)
      
      comb <- comb %>% group_by(sequence) %>% summarise(across(primers[row,]$Primer, ~sum(.x, na.rm = TRUE)))
      
      if (is.null(out)) {
        out <- comb
      } else {
        out <- merge(out, comb, all = TRUE, by = "sequence")
        out[is.na(out)] <- 0
      }
    }
    
    return(out)
  })
  
  amino <- reactive({
    req(input$r1_file)
    req(input$r2_file)
    req(input$metadata_file)
    
    dna <- merged()
    
    dna$sequence <- translate_reads(DNAStringSet(dna$sequence))
    
    dna <- dna %>% group_by(sequence) %>% summarise(across(where(~ !is.character(.) & !is.factor(.)), ~sum(.x, na.rm = TRUE)))
    
    dna <- dna[dna$sequence != "", ]
    
    dna
  })
  
  output$sequencetable <- DT::renderDataTable({
    sequences <- merged()
    sequences$sequence <- sapply(sequences$sequence, function(x) {
      if (nchar(x) > 78) {paste0(substr(x, 1, 78), "...")}
      else {x}
    })
    
    if (nrow(sequences) == 0) {
      return(data.frame(message = "no data"))
    } else {
      return(sequences)
    }
  })
  
  output$aminoacidtable <- DT::renderDataTable({
    sequences <- amino()
    
    if (nrow(sequences) == 0) {
      return(data.frame(message = "no data"))
    } else {
      return(sequences)
    }
  })
  
  output$plot1 <- renderPlot({
    req(input$r1_file)
    req(input$r2_file)
    req(input$metadata_file)
    
    ggplot(data = data.frame(x = nchar(merged()$sequence)), aes(x = x)) +
      geom_density(alpha = 0.5) +
      labs(x = "length of DNA sequences", title = "distribution of DNA sequence lengths")
  }, res = 96)
  
  output$plot2 <- renderPlot({
    req(input$r1_file)
    req(input$r2_file)
    req(input$metadata_file)
    
    ggplot(data = data.frame(x = nchar(amino()$sequence)), aes(x = x)) +
      geom_density(alpha = 0.5) +
      labs(x = "length of amino acid sequences", title = "distribution of amino acid sequence lengths")
  }, res = 96)
  
  output$savebutton1 <- downloadHandler(
    filename = function() {
      paste0("dna_sequences.csv")
    },
    content = function(file) {
      # Write the dataset to the `file` that will be downloaded
      write.csv(merged(), file = file, row.names = FALSE)
      output$savetext1 <- renderText({
        paste("Data saved to:", file)
      })
    }
  )
  
  output$savebutton2 <- downloadHandler(
    filename = function() {
      paste0("amino_acid_sequences.csv")
    },
    content = function(file) {
      write.csv(amino(), file = file, row.names = FALSE)
      
      output$savetext2 <- renderText({
        paste("Data saved to:", file)
      })
    }
  )
})
