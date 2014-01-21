library(shiny)

shinyUI(
    navbarPage("Zectran the explorer",
               header=p("contact: Marcin Cieslik (mcieslik@med.umich.edu)", align="right"),

               ## IDS
               tabPanel("Individual Dataset Selection",
                        sidebarLayout(
                            sidebarPanel(
                                conditionalPanel(condition="input.tabs_ids == 'Dataset Summary'",
                                                 h5("Select Dataset(s)"),
                                                 wellPanel(
                                                     uiOutput("analysis_select_ids"),
                                                     uiOutput("value_select_ids"),
                                                     uiOutput("study_select_ids")
                                                     )
                                                 ),
                                conditionalPanel(condition="input.tabs_ids == 'Cohort Builder'",
                                                 fluidRow(
                                                     column(width=6,
                                                            h5("(sub)cohort A"),
                                                            wellPanel(
                                                                uiOutput("cohort_a_select_ids"),
                                                                textInput("samplecolumn_a_query_ids", "sample filter column",
                                                                          value="sample_type_id"),
                                                                uiOutput("samplefilter_a_select_ids"),
                                                                textInput("patientcolumn_a_query_ids", "patient filter column",
                                                                          value="gleason_score"),
                                                                uiOutput("patientfilter_a_select_ids"),
                                                                uiOutput("patientselect_a_select_ids")
                                                                       )
                                                            ),
                                                     column(width=6,
                                                                   h5("(sub)cohort B"),
                                                                   wellPanel(
                                                                       uiOutput("cohort_b_select_ids"),
                                                                       textInput("samplecolumn_b_query_ids", "Sample Column",
                                                                                 value="sample_type_id"),
                                                                       uiOutput("samplefilter_b_select_ids"),
                                                                       textInput("patientcolumn_b_query_ids", "patient filter column",
                                                                                     value="gleason_score"),
                                                                       uiOutput("patientfilter_b_select_ids"),
                                                                       uiOutput("patientselect_b_select_ids")
                                                                       )
                                                            )
                                                     ),
                                                 actionButton("cohorts_ab_action_ids", "update cohorts")
                                                 )
                                    ),
                            
                            mainPanel(
                                    tabsetPanel(
                                        tabPanel("Dataset Summary",
                                                 h4("Dataset Summary"),
                                                 tableOutput("dataset_table_ids")
                                                 ),
                                        tabPanel("Cohort Builder",
                                                 h4("Patient Inspector"),
                                                 dataTableOutput("patient_clinical_ids"),
                                                 h4("Cohort Summary"),
                                                 tableOutput("cohort_table_ids")
                                                 ),
                                        id="tabs_ids")
                                )
                            )),

               ## PDS
               tabPanel("Paired Dataset Selection",
                        sidebarLayout(
                            sidebarPanel(
                                conditionalPanel(condition="input.tabs_pds == 'Dataset Summary'",
                                                 h5("Build Paired Dataset"),
                                                 wellPanel(
                                                     fluidRow(
                                                         column(width=6,
                                                                h5("dataset A"),
                                                                wellPanel(
                                                                    uiOutput("analysis_a_select_pds"),
                                                                    uiOutput("value_a_select_pds"),
                                                                    uiOutput("study_a_select_pds")
                                                                    )
                                                                ),
                                                         column(width=6,
                                                                h5("dataset B"),
                                                                wellPanel(
                                                                    uiOutput("analysis_b_select_pds"),
                                                                    uiOutput("value_b_select_pds"),
                                                                    uiOutput("study_b_select_pds")
                                                                    )
                                                                )
                                                     ))
                                                 ),
                                conditionalPanel(condition="input.tabs_pds == 'Cohort Builder'",
                                                 h5("Build Cohort"),
                                                 wellPanel(
                                                     uiOutput("cohort_select_pds"),
                                                     textInput("samplecolumn_query_pds", "sample filter column",
                                                               value="sample_type_id"),
                                                     uiOutput("samplefilter_select_pds"),
                                                     textInput("patientcolumn_query_pds", "patient filter column",
                                                               value="gleason_score"),
                                                     uiOutput("patientfilter_select_pds"),
                                                     uiOutput("patientselect_select_pds")
                                                     )
                                                 )
                                    ),
                            
                            mainPanel(
                                    tabsetPanel(
                                        tabPanel("Dataset Summary",
                                                 h4("Dataset Summary")
                                                 ),
                                        tabPanel("Cohort Builder",
                                                 h4("Patient Inspector"),
                                                 h4("Cohort Summary")
                                                 ),
                                        id="tabs_pds")
                                )
                        )),
                        

               ## MA
               tabPanel("Methylation Analysis",
                        sidebarLayout(
                            sidebarPanel(
                                conditionalPanel(condition="input.tabs_ma == 'Gene Query'",
                                                 h5("Select Gene"),
                                                 wellPanel(
                                                     fluidRow(
                                                         column(width=6, textInput("hgnc_query_ma", "HUGO (HGNC) gene name")),
                                                         column(width=6, textOutput("name_text_ma"))
                                                         ),
                                                     actionButton("ucscmodel_action_ma", "plot gene model")
                                                     ),
                                                 h5("Refine Region"),
                                                 wellPanel(
                                                     uiOutput("ucsc_select_ma"),
                                                     selectInput("region_select_ma", "select region",
                                                                 choices=c(
                                                                     "TSS slim (500, 200)",
                                                                     "TSS wide (2000, 500)",
                                                                     "gene regulatory region")),
                                                     uiOutput("refinedregion_query_ma")
                                                     ),
                                                 h5("Select Probes"),
                                                 wellPanel(
                                                     uiOutput("cpg_select_ma"),
                                                     actionButton("probes_action_ma", "plot transcript & probes")
                                                     )
                                                 ),
                                conditionalPanel(condition="input.tabs_ma == 'Individual Dataset Analysis'",
                                                 h5("Individual")
                                                 ),
                                conditionalPanel(condition="input.tabs_ma == 'Paired Dataset Analysis'",
                                                 h5("Paired")
                                                 )
                                             ),
                            mainPanel(
                                tabsetPanel(
                                    tabPanel("Gene Query",
                                             h4("UCSC Gene Model"),
                                             plotOutput("ucscmodel_plot_ma", width=800),
                                             h4("Selected UCSC Transcript"),
                                             plotOutput("probes_plot_ma", width=800, height=300)
                                             ),
                                    tabPanel("Individual Dataset Analysis",
                                             h4("Individual Dataset Analysis"),
                                             htmlOutput("dataset_text_ma")
                                             ),
                                    tabPanel("Paired Dataset Analysis",
                                             h4("Paired Dataset Analysis")
                                             ),
                                    id="tabs_ma"
                                    )
                                )
                            )
                        ),
               
               ## EA
               tabPanel("Expression Analysis",
                        sidebarLayout(
                            sidebarPanel(
                                conditionalPanel(condition="input.tabs_ea == 'Gene Query'",
                                                 h5("Select Gene"),
                                                 wellPanel(
                                                     fluidRow(
                                                         column(width=6, textInput("hgnc_query_ea", "HUGO (HGNC) gene name")),
                                                         column(width=6, textOutput("name_text_ea"))
                                                         ),
                                                     uiOutput("ucsc_select_ea"),
                                                     actionButton("ucscmodel_action_ea", "plot gene model")
                                                     )
                                                 ),
                                conditionalPanel(condition="input.tabs_ea == 'Individual Dataset Analysis'",
                                                 h5("Dataset"),
                                                 htmlOutput("dataset_text_ea"),
                                                 h5("Cohorts"),
                                                 htmlOutput("cohorts_text_ea")
                                                 ),
                                conditionalPanel(condition="input.tabs_ea == 'Paired Dataset Analysis'",
                                                 h5("Paired")
                                                 )
                                ),
                            mainPanel(
                                tabsetPanel(
                                    tabPanel("Gene Query",
                                             h4("UCSC Gene Model"),
                                             plotOutput("ucscmodel_plot_ea", width=800)
                                             ),
                                    tabPanel("Individual Dataset Analysis",
                                             h4("Individual Dataset Analysis")
                                             ),
                                    tabPanel("Paired Dataset Analysis",
                                             h4("Paired Dataset Analysis")
                                             ),
                                    id="tabs_ea"
                                    )
                                )
                            )
                        )
               )
    )
print("ui")




## conditionalPanel(condition="input.tabs == 'Gene Query'",
##                  h5("Select Transcript"),
##                  wellPanel(
##                      fluidRow(
##                          column(width=6, textInput("hgnc_query_gq", "HUGO (HGNC) gene name")),
##                          column(width=6, textOutput("name_text_gq"))
##                          ),
##                      actionButton("ucscmodel_action_gq", "plot gene model"),
##                      hr(),
##                      uiOutput("ucsc_select_gq")
##                      )
##                  )
## conditionalPanel(condition="input.paired_tabs == 'Gene Query'",
##                  h5("Select Transcript"),
##                  wellPanel(
##                      ## fluidRow(
##                      ##     column(width=6, textInput("hgnc_query_gq", "HUGO (HGNC) gene name")),
##                      ##     column(width=6, textOutput("name_text_gq"))
##                      ##     ),
##                      ## actionButton("ucscmodel_action_gq", "plot gene model"),
##                      ## hr(),
##                      ## uiOutput("ucsc_select_gq")
##                      )
##                  )
## fluidRow(
##     column(width=6,
##            h5("(sub)cohort A"),
##            wellPanel(
##                       )
##            ),
##     column(width=6,
##                   h5("(sub)cohort B"),
##                   wellPanel(
##                       uiOutput("cohort_b_select_ids"),
##                       textInput("samplecolumn_b_query_ids", "Sample Column",
##                                 value="sample_type_id"),
##                       uiOutput("samplefilter_b_select_ids"),
##                       textInput("patientcolumn_b_query_ids", "patient filter column",
##                                     value="gleason_score"),
##                       uiOutput("patientfilter_b_select_ids"),
##                       uiOutput("patientselect_b_select_ids")
##                       )
##            )
##     ),
## actionButton("cohorts_ab_action_ids", "update cohorts")
## tabPanel("Gene Query",
##      h4("UCSC Gene Model"),
##      plotOutput("ucscmodel_plot_gq", width=800),
##      textOutput("ucscselect_text_gq")
##      ),

