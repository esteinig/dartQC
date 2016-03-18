library(shiny)

library(shinydashboard)

shinyUI(
  dashboardPage(
    
    skin = 'black',
    
    dashboardHeader(title = "DartQC", titleWidth = 250),
    
    dashboardSidebar(width = 250, 
                     
                     sidebarMenu(
                       br(),
                       menuItem(" Files", tabName = "input", icon = icon("file-o")),
                       menuItem(" DartQC", tabName = "quality_control", icon = icon("cogs")),
                       menuItem(" Results", icon=icon("bar-chart"),
                                menuSubItem(" Report", tabName = "reports", icon = icon("newspaper-o")),
                                menuSubItem(" Graphs", tabName="graphs", icon = icon("area-chart"))),
                       menuItem(" About", tabName = "about", icon=icon("heart", lib = "glyphicon"))
                     )
    ),
    
    dashboardBody(
      tabItems(
       
        tabItem(tabName="input",
                fluidRow(
                  
                          column(4, 
                                 box(align='center', width='100%', br(),
                                        h1('Files'),
                                        hr(),
                                        br(),
                                        helpText("Input format is a comma-delimited version of the Excel spreadsheet from DArT. Currently only supports double-row format and SNPs. For appropriate output formatting in STRUCTURE and PLINK formats, you can also upload a comma-delimited population file, containing two columns with header: ID and Population."),
                                        br(),
                                        br(),
                                        br(),
                                        hr(),
                                        h2('Data'),
                                        fileInput("data_file", NULL, multiple=F),
                                        h2('Populations'),
                                        fileInput("pop_file", NULL, multiple=F),
                                        hr(),
                                        br(),
                                        div(textInput("project", "Project", "Dart"),  style="text-align:center"),
                                        br(),
                                        selectInput('data_type','Data Type', choices=list("SNP" = "SNP",
                                                                                          "DART" = "DART"), selected='SNP'),
                                        selectInput('data_format','Data Format', choices=list("Double Row" = "double",
                                                                                       "Single Row" = "single"), selected='double'),
                                        selectInput('output_format','Output Format', choices=list("PLINK" = "plink",
                                                                                           "STRUCTURE" = "structure"), selected='structure'),
                                        br(),
                                        checkboxInput('keep', 'Keep Temporary Files'),
                                        checkboxInput('verbose', 'Verbose')
                                        
                                        
                                  )
                          ),
                          
                          column(4, 
                                 box(align='center', width='100%', br(),
                                     h1('Formatting'),
                                     hr(),
                                     br(),
                                     helpText("Details of the data spreadsheet formatting for proper parsing of data from DArT. Make sure you have the right rows and column selected,
                                              otherwise the pipeline will not return the correct calculations or return Errors. Hover over the input boxes to get details about the parameters."),
                                     br(),
                                     br(),
                                     br(),
                                     hr(),
                                     h2('Rows'),
                                     numericInput('id_row', "Sample IDs", value = 7, min = 1),
                                     numericInput('data_row', "Start of Allele Data", value = 8, min = 1),
                                     h2('Columns'),
                                     numericInput('id_col', "Allele IDs", value = 1, min = 1),
                                     numericInput('clone_col', "Clone IDs", value = 2, min = 1),
                                     numericInput('seq_col', 'Allele Sequences', value = 3, min = 1),
                                     numericInput('snp_col', 'SNP', value = 4, min = 1),
                                     numericInput('snp_pos_col', 'SNP Position', value = 5, min = 1),
                                     numericInput('call_rate_dart', 'Call Rate', value = 6, min = 1),
                                     numericInput('one_ratio_ref_col', 'OneRatio REF', value = 7, min = 1),
                                     numericInput('one_ratio_snp_col', 'OneRatio SNP', value = 8, min = 1),
                                     numericInput('freq_homozygous_ref_col', 'Frequency Homozygous REF', value = 9, min = 1),
                                     numericInput('freq_homozygous_snp_col', 'Frequency Homozygous SNP', value = 10, min = 1),
                                     numericInput('freq_heterozygous_col','Frequency Heterozygous', value = 11, min = 1),
                                     numericInput('pic_ref_col', 'PIC REF', value = 12, min = 1),
                                     numericInput('pic_snp_col', 'PIC SNP', value = 13, min = 1),
                                     numericInput('average_pic_col', 'Average PIC', value = 14, min = 1),
                                     numericInput('average_read_count_ref_col', 'Average Read Count REF', value = 15, min = 1),
                                     numericInput('average_read_count_snp_col', 'Average Read Count SNP', value = 16, min = 1),
                                     numericInput('rep_col', 'Average Replication',value = 17, min = 1),
                                     numericInput('call_col', 'Start of Allele Calls', value = 18, min = 1)
                                 )
                          ),
                          
                          column(4, 
                                 box(align='center', width='100%', br(),
                                     h1('Alleles'),
                                     hr(),
                                     br(),
                                     helpText("Here you can specify the encoding for homozygous major, homozygous minor, heterozygous and missing alleles."),
                                     br(),
                                     br(),
                                     br(),
                                     br(),
                                     br(),
                                     hr(),
                                     h2('Encoding'),
                                     textInput('homozygous_major', 'Homozygous Major', '10'),
                                     textInput('homozygous_minor', 'Homozygous Minor', '01'),
                                     textInput('heterozygous', 'Heterozygous', '11'),
                                     textInput('missing', 'Missing', '--')
                                 )
                          )
                          
                )
        ),
      
        tabItem(tabName="quality_control",
                fluidRow(
                  column(6, 
                         box(align='center', width='100%', br(),
                             h1('Quality Control'),
                             hr(),
                             br(),
                             helpText("Setting and filter values for running DartQC. Reference allele sequences can be clustered by identity using CD-HIT, and selectors for retaining best-picks when removing identitcal sequences and duplicate clones can be selected under Identity. Filter thresholds are currently implemented only for SNPs."),
                             br(),
                             hr(),
                             h3("Identity Clustering"),
                             br(),
                             checkboxInput('identity', 'Remove Identical Sequences', TRUE),
                             numericInput("identity_threshold", "Identity Threshold", value=0.95, min=0, max=1),
                             selectInput('identity_selector','Retain best identical sequence by:', choices=list("Minor Allele Frequency" = "maf",
                                                                                       "Call Rate" = "call_rate",
                                                                                       "Replication Average" = "rep"), selected='maf'),
                             selectInput('clone_selector','Retain best duplicate clone by:', choices=list("Minor Allele Frequency" = "maf",
                                                                                                      "Call Rate" = "call_rate",
                                                                                                      "Replication Average" = "rep"), selected='maf'),
                             br(),
                             hr(),
                             h3("Filters"),
                             br(),
                             numericInput("maf", "Minor Allele Frequency (<=)", value=0.02, min=0, max=1),
                             
                             numericInput("call", "Call Rate (<=)", value=0.70, min=0, max=1),
                             
                             numericInput("rep", "Replication Average (<=)", value=0.95, min=0, max=1),
                             br(),
                             hr(),
                             br()
                             
                             
                         )
                  ),
                  column(6, box(align='center', width='100%', br(),
                            h1("Pipeline"),
                            hr(),
                            actionButton('run_qc', 'Start DartQC', icon=icon("circle-arrow-right", lib = "glyphicon")),
                            br(),
                            hr(),
                            htmlOutput("log")
                            
                         )
                  )
                  
                )
        )
      
      )
    )
  )
)