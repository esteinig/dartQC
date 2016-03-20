library(shiny)
library(shinydashboard)


shinyServer(function(input, output, session) {
  
  base.path = c('~')
  names(base.path) = 'Home'
  
  shinyFileChoose(input, 'data_file', root=base.path, session=session)
  shinyFileChoose(input, 'pop_file', root=base.path, session=session)
  
  # Function to construct Configuration File
  
  get.config <- function(data.path, pop.path){
    
    config.list = list(
      
      "project" = input$project,
      "data_file" = data.path$datapath,
      "data_format" = input$data_format,
      "data_type" = input$data_type,
      "pop_file" = pop.path$datapath,
      "output_format" = input$output_format,
      "maf" = input$maf,
      "call" = input$call,
      "rep" = input$rep,
      "seq_identity" = input$identity_threshold,
      "identity_selector" = input$identity_selector,
      "clone_selector" = input$clone_selector,
      "one_ratio_ref" = -1,
      "one_ratio_snp" = -1,
      "freq_homozygous_ref" = -1,
      "freq_homozygous_snp" = -1,
      "pic_ref" = -1,
      "pic_snp" = -1,
      "average_pic" = -1,
      "average_read_count_ref" = -1,
      "average_read_count_snp" = -1,
      "data_row" = input$data_row,
      "sample_row" = input$id_row,
      "id_col" = input$id_col,
      "clone_col" = input$clone_col,
      "seq_col" = input$seq_col,
      "snp_col" = input$snp_col,
      "snp_pos_col" = input$snp_pos_col,
      "call_rate_dart_col" = input$call_rate_dart,
      "one_ratio_ref_col" = input$one_ratio_ref_col,
      "one_ratio_snp_col" = input$one_ratio_snp_col,
      "freq_homozygous_ref_col" = input$freq_homozygous_ref_col,
      "freq_homozygous_snp_col" = input$freq_homozygous_snp_col,
      "freq_heterozygous_col" = input$freq_heterozygous_col,
      "pic_ref_col" = input$pic_ref_col,
      "pic_snp_col" = input$pic_snp_col,
      "average_pic_col" = input$average_pic_col,
      "average_read_count_ref_col" = input$average_read_count_ref_col,
      "average_read_count_snp_col" = input$average_read_count_snp_col,
      "rep_col" = input$rep_col,
      "call_col" = input$call_col,
      "homozygous_major" = input$homozygous_major,
      "homozygous_minor" = input$homozygous_minor,
      "heterozygous" = input$heterozygous,
      "missing" = input$missing,
      "verbose" = input$verbose,
      "keep" = input$keep
    )
    
    print(names(config.list))
    print(unlist(config.list))
    
  }
  
  observeEvent(input$run_qc, {
    
    data.path = parseFilePaths(base.path, input$data_file)
    pop.path = parseFilePaths(base.path, input$pop_file)
    
    print(data.path)
    print(pop.path)
    
    withProgress({
      
      incProgress(1)
      
      config.list <- get.config(data.path, pop.path)
      
      incProgress(1)
      
      out <- system2(input$python_path, "dartQC.py", stdout = TRUE, stderr = TRUE)
      
      incProgress(1)
      
    }, max = 3, value=0, message="Running DartQC script in Python...")
    
    output$log <- renderUI({
      
      text <- paste(out, collapse = "<br/>")
      HTML(text)
      
    })
    
  })
  
})