library(shiny)

shinyUI(pageWithSidebar(

    headerPanel("Zectran the explorer"),

    sidebarPanel(

        conditionalPanel(condition="input.tabs == 'Dataset Summary'",
                         h5("Select Dataset(s)"),
                         wellPanel(
                             uiOutput("analysis_select_ds"),
                             uiOutput("value_select_ds"),
                             uiOutput("study_select_ds")
                             )
                         ),
        
        conditionalPanel(condition="input.tabs == 'Cohort Builder'",
                     tags$div(class = "row-fluid",
                              tags$div(class = "span6",
                                       h5("(sub)cohort A"),
                                       wellPanel(
                                           uiOutput("cohort_a_select_cb"),
                                           textInput("samplecolumn_a_query_cb", "sample filter column",
                                                     value="sample_type_id"),
                                           uiOutput("samplefilter_a_select_cb"),
                                           textInput("patientcolumn_a_query_cb", "patient filter column",
                                                     value="gleason_score"),
                                           uiOutput("patientfilter_a_select_cb"),
                                           uiOutput("patientselect_a_select_cb")
                                           )
                                       ),
                              tags$div(class = "span6",
                                       h5("(sub)cohort B"),
                                       wellPanel(
                                           uiOutput("cohort_b_select_cb"),
                                           textInput("samplecolumn_b_query_cb", "Sample Column",
                                                     value="sample_type_id"),
                                           uiOutput("samplefilter_b_select_cb"),
                                           textInput("patientcolumn_b_query_cb", "patient filter column",
                                                     value="gleason_score"),
                                           uiOutput("patientfilter_b_select_cb"),
                                           uiOutput("patientselect_b_select_cb")
                                           )
                                       )
                              ),
                         actionButton("cohorts_ab_action_cb", "update cohorts")
                         ),
        
        conditionalPanel(condition="input.tabs == 'Gene Query'",
                         h5("Select Gene"),
                         wellPanel(
                             textInput("hgnc_query_gq", "HUGO (HGNC) gene name"),
                             uiOutput("select_ucsc_gq")
                             #actionButton("plot_ensg_gq", "plot gene model")
                             )
                         ##
                         ## h5("Refine region"),
                         ## wellPanel(
                         ##     uiOutput("select_enst_gq"),
                         ##     selectInput("region_select_gq", "select region",
                         ##                 choices=c(
                         ##                     "TSS slim (500, 200)",
                         ##                     "TSS wide (2000, 500)",
                         ##                     "gene regulatory region")),
                         ##     uiOutput("query_region_limit_gq"),
                         ##     actionButton("plot_enst_gq", "plot selected region")
                         ##     )
                         )
    ),
    
    mainPanel(
        
        tabsetPanel(

            tabPanel("Dataset Summary",
                     h4("Dataset Summary"),
                     tableOutput("dataset_table_ds")
                     ),

            tabPanel("Cohort Builder",
                     h4("Patient Inspector"),
                     dataTableOutput("patient_clinical_cb"),
                     h4("Cohort Summary"),
                     tableOutput("cohort_table_cb")
                     ),

            tabPanel("Gene Query",
                     h4("Ensembl Gene Model"),                     
                     ##plotOutput("ensts_tracks_gq", width=800),

                     h4("Selected Transcript")
                     ##plotOutput("ensts_select_tracks_gq", width=800, height=300)
                     ),
            
            id="tabs")
    )
))
print("ui")
