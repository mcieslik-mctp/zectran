library(shiny)

shinyUI(pageWithSidebar(

    headerPanel("Zectran the explorer"),

    sidebarPanel(

        conditionalPanel(condition="input.tabs == 'Select Dataset'",
                         h5("Select Dataset(s)"),
                         wellPanel(
                             uiOutput("analysis_select_ds"),
                             uiOutput("value_select_ds"),
                             uiOutput("study_select_ds")
                             )
                         ),
        
        conditionalPanel(condition="input.tabs == 'Cohort Builder'",
                         h5("(sub)cohort A"),
                         wellPanel(
                             uiOutput("cohort_a_select_cb"),
                             textInput("samplecolumn_a_query_cb", "sample filter column",
                                       value="sample_type_id"),
                             uiOutput("samplefilter_a_select_cb"),
                             textInput("patientcolumn_a_query_cb", "patient filter column",
                                       value="gleason_score"),
                             uiOutput("patientfilter_a_select_cb")
                             ),
                         h5("(sub)cohort B"),
                         wellPanel(
                             uiOutput("cohort_b_select_cb"),
                             textInput("samplecolumn_b_query_cb", "Sample Column",
                                       value="sample_type_id"),
                             uiOutput("samplefilter_b_select_cb"),
                             textInput("patientcolumn_b_query_cb", "patient filter column",
                                       value="gleason_score"),
                             uiOutput("patientfilter_b_select_cb")
                             )
                         ),
        
        conditionalPanel(condition="input.tabs == 'Patient Summary'",
                         h5("Select Patient"),
                         wellPanel(
                             uiOutput("patient_a_select_ps"),
                             uiOutput("patient_b_select_ps")
                             )
                         ),
                
        conditionalPanel(condition="input.tabs == 'Gene Query'",
                         h5("Select Gene"),
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
                         )
        
    ),
    
    mainPanel(
        
        tabsetPanel(

            tabPanel("Select Dataset",
                     h4("Dataset Summary")
                     ),

            tabPanel("Cohort Builder",
                     h4("Cohort Summary")
                     ),

            tabPanel("Patient Summary",
                     h4("Patient Summary"),
                     dataTableOutput("patient_clinical_ps")
                     ),            
            
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
            
            id="tabs")
    )
))
print("ui")
