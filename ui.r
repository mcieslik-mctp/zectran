library(shiny)

shinyUI(pageWithSidebar(

    headerPanel("Zectran: the methyl explorer"),

    sidebarPanel(

        conditionalPanel(condition="input.tabs == 'Gene Query'",
                         h5("Select gene"),
                         wellPanel(
                             textInput("query_hgnc_gq", "HUGO (HGNC) gene name"),
                             uiOutput("select_ensg_gq"),
                             actionButton("plot_ensg_gq", "plot gene model")
                             ),
                         ##
                         h5("Refine region"),
                         wellPanel(
                             uiOutput("select_enst_gq"),
                             selectInput("region_select_gq", "select region",
                                         choices=c(
                                             "TSS slim (500, 200)",
                                             "TSS wide (2000, 500)",
                                             "gene regulatory region")),
                             uiOutput("query_region_limit_gq"),
                             actionButton("plot_enst_gq", "plot selected region")
                             )
                         ),

        conditionalPanel(condition="input.tabs == 'Probe Analysis'",
                         h5("Probe Analysis"),
                         wellPanel(
                             selectInput("analysis_select_gq", "select analysis",
                                         choices=c("", "probe summary", "probe correlation",
                                             "probe boxplot", "methylation track")),
                             uiOutput("select_cg_gq"),
                             conditionalPanel(condition=
                                              "input.analysis_select_gq == 'probe boxplot'",
                                              selectInput("comparison_select_gq", "comparison",
                                                          choices=c("tumor", "normal",
                                                              "tumor-normal",
                                                              "tumor-normal-stage"),
                                                          selected="tumor-normal")
                                              ),
                             uiOutput("cohort_select_gq"))
                         ),

        conditionalPanel(condition="input.tabs == 'Sample Analysis'",
                         h5("Sample Analysis"),
                         wellPanel(
                             selectInput("analysis_select_sa", "select analysis",
                                         choices=c("PCA", "RPMM clusters")),

                             selectInput("probetypes_select_sa", "select probe location",
                                         choices=c("gene regulatory regions", "enhancers")),

                             selectInput("probetypes_select_sa", "select probe types",
                                         choices=c("CG Island", "CG shore")),
                             
                             sliderInput("cg_variance_sa", "probe variance percentile:",
                                         min=0, max=100, value=c(75,100), step=5),
                             uiOutput("cohort_select_sa"),
                             selectInput("samples_select_sa", "select samples",
                                         choices=c("all", "tumor", "normal")),
                             selectInput("comparison_select_sa", "select comparison",
                                         choices=c("tumor vs normal", "stage")),
                             
                             conditionalPanel(condition="input.analysis_select_sa == 'PCA'",
                                              actionButton("update_plot_pca_sa", "update PCA")),
                             
                             conditionalPanel(condition="input.analysis_select_sa == 'RPMM clusters'",
                                          actionButton("update_plot_rpmm_sa", "update RPMM"))
                             )
                         ),

        conditionalPanel(condition="input.tabs == 'Transcript Expression'",
                         h5("Select Table"),
                         wellPanel(
                             uiOutput("analysis_select_te"),
                             uiOutput("value_select_te"),
                             uiOutput("study_select_te"),
                             uiOutput("cohort_select_te")
                             ),
                         h5("Select Data"),
                         wellPanel(
                             uiOutput("sample_select_te"),
                             textInput("id_query_te", "id")
                             )
                         ),
        
        conditionalPanel(condition="input.tabs == 'Dataset Summary'",
                         h5("Select cohorts"),
                         wellPanel(
                             uiOutput("cohort_select_ds")
                             )
                         )
    ),
    
    mainPanel(
        
        tabsetPanel(
            
            tabPanel("Gene Query",
                     h4("Ensembl Gene Model"),                     
                     plotOutput("ensts_tracks_gq", width=800),

                     h4("Selected Transcript"),                     
                     plotOutput("ensts_select_tracks_gq", width=800, height=300)
                     ),

            tabPanel("Probe Analysis",
                     conditionalPanel(condition="input.analysis_select_gq == ''",
                                      h5("select analysis, probes, and cohorts")),
                     conditionalPanel(condition="input.analysis_select_gq == 'probe summary'",
                                      h4("CpG info"),
                                      tableOutput("cpg_info"),
                                      br(),
                                      h4("CpG islands info"),
                                      tableOutput("cpg_islands_info")),                
                     conditionalPanel(condition="input.analysis_select_gq == 'probe boxplot'",
                                      h4("Methylation Boxplot"),
                                      plotOutput("cg_boxplot_gq", width=800, height=600)),
                     conditionalPanel(condition="input.analysis_select_gq == 'probe correlation'",
                                      h4("Probe Correlation"),
                                      plotOutput("cg_corrplot_gq", width=500)),
                     conditionalPanel(condition="input.analysis_select_gq == 'methylation track'",
                                      h4("Methylation Track"),
                                      plotOutput("cg_track_gq", width=800, height=400))
                     ),

            

            tabPanel("Sample Analysis",
                     conditionalPanel(condition="input.analysis_select_sa == ''",
                                      h5("select analysis, samples")),

                     conditionalPanel(condition="input.analysis_select_sa == 'PCA'",
                                      plotOutput("plot_pca_sa", width=800)
                                      ), 
                     
                     conditionalPanel(condition="input.analysis_select_sa == 'RPMM clusters'",
                                      plotOutput("plot_rpmm_sa", width=800)
                                      )

                     ),

            tabPanel("Transcript Expression",
                     h5("KLK3 - uc010eof.1"),
                     textOutput("data_te")
                     ),
 
            
            tabPanel("Dataset Summary",
                     h4("Cohort Summary"),
                     verbatimTextOutput("cohort_description"),
                     tableOutput("cohort_summary"),
                     h4("Cohort Samples"),
                     dataTableOutput("cohort_samples"),
                     h4("Query Data"),
                     tableOutput("debug_data"),
                     h4("ENST Debug probes"),
                     textOutput("enst_debug_probes")
                     ),
            
            id="tabs")
    )
))
print("ui")
source("thetbl.r")
