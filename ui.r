library(shiny)

dataset_tab = tabPanel("Dataset Summary", 
    h4("Dataset Summary"),
    tableOutput("dataset_table_ds"))

cohort_tab = tabPanel("Cohort Builder",
    h4("Patient Inspector"),
    dataTableOutput("patient_clinical_cb"),
    h4("Cohort Summary"),
    tableOutput("cohort_table_cb"))

meth_tab = tabPanel("Analysis")

exp_tab = tabPanel("Analysisx")

methexp_tab = tabPanel("Analysisz")

dataset_well = wellPanel(
    uiOutput("analysis_select_ds"),
    uiOutput("value_select_ds"),
    uiOutput("study_select_ds")
    )


shinyUI(
    navbarPage("Zectran the explorer",
               
               header=p("Marcin Cieslik (mcieslik@med.umich.edu)", align="right"),
               
               tabPanel("Dataset Selection",
                        sidebarLayout(
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
                                        id="tabs")
                                )
                            )),
                        
               tabPanel("Methylation Analysis",
                        sidebarLayout(
                            sidebarPanel(
                                h5("Select Gene"),
                                wellPanel(
                                    tags$div(class = "row-fluid",
                                             tags$div(class = "span6", textInput("hgnc_query_me", "HUGO (HGNC) gene name")),
                                             tags$div(class = "span6", textOutput("name_text_me"))
                                             ),
                                    actionButton("ucscmodel_action_me", "plot gene model")
                                    ),
                                                 h5("Refine Region"),
                                wellPanel(
                                    uiOutput("ucsc_select_me"),
                                    selectInput("region_select_me", "select region",
                                                choices=c(
                                                    "TSS slim (500, 200)",
                                                    "TSS wide (2000, 500)",
                                                    "gene regulatory region")),
                                    uiOutput("refinedregion_query_me"),
                                    actionButton("ucscregion_plot_me", "plot refined region")
                                    )
                                ),
                            
                            mainPanel(
                                h4("UCSC Gene Model"),
                                plotOutput("ucscmodel_plot_me", width=800),
                                h4("Selected UCSC Transcript"),
                                plotOutput("ucsctranscript_plot_me", width=800, height=300)
                                )
                            )
                        ),
               tabPanel("Expression Analysis"
                        ##
                   )
               )
    )
print("ui")
