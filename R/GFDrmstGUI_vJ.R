library("plyr")
library("GFDmcv")
library("lpSolve")
library("shinyMatrix")
library("shinyWidgets")
library("GFDrmst")

GFDrmtlGUI <- function(){
  requireNamespace("shiny", quietly = FALSE)
  if (!("package:shiny" %in% search())) {
    attachNamespace("shiny")
  }
  requireNamespace("shinyMatrix", quietly = FALSE)
  if (!("package:shinyMatrix" %in% search())) {
    attachNamespace("shinyMatrix")
  }
  requireNamespace("tippy", quietly = FALSE)
  if (!("package:tippy" %in% search())) {
    attachNamespace("tippy")
  }
  ui <- fluidPage(theme = shinythemes::shinytheme("cerulean"),
                  shinyjs::useShinyjs(), titlePanel("Tests for GFDrmtl"),
                  sidebarLayout(sidebarPanel(splitLayout(cellWidths = c("40%","15%","20%","20%","0%"), fileInput("infile",
                                                                                                                 "Choose CSV File", accept = c("text/csv", "text/comma-separated-values,text/plain",".csv")),
                                                         actionButton("infoButton2", "", icon = icon("info-circle")),
                                                         checkboxInput("header", "Header", TRUE),
                                                         selectInput("sep", "Separator in csv", c(",", ";", ".", "|")),
                                                         tippy::tippy_this("infoButton2", "The csv file should be structured in the following way: different rows must contain the data for the independent individuals, different columns contain the different values (at least the event time, censoring status, and one factor variable).   ")
                  ),
                  tags$head(tags$style(HTML("\n .shiny-split-layout > div {\n overflow: visible;\n }\n"))),
                  tags$style(HTML("\n input[type=number] {\n -moz-appearance:textfield;\n }\n input[type=number]::{\n -moz-appearance:textfield;\n }\n
                  input[type=number]::-webkit-outer-spin-button,\n input[type=number]::-webkit-inner-spin-button {\n -webkit-appearance: none;\n margin: 0;\n }\n")),
                  h3(id = "titleLoadData", "Load dataset first!", style = "color:red"),
                  splitLayout(cellWidths = c("80%"), shinyjs::hidden(selectInput("Method", "Select Testing Method:", c(`Asymptotic` = "asymptotic", `Permutation` = "permutation"), selected = "permutation"))),
                  splitLayout(cellWidths = c("50%", "50%"), uiOutput(outputId = "dynamicInput"), uiOutput(outputId = "dynamicInput2")),
                  splitLayout(cellWidths = c("50%"), shinyjs::hidden(checkboxInput("stepwise", "Stepwise extension", TRUE))),
                  splitLayout(cellWidths = c("50%"), shinyjs::hidden(checkboxInput("plots", "Plot the confidence intervals", FALSE))),
                  splitLayout(cellWidths = c("50%", "5%", "45%"), uiOutput(outputId = "dynamicInput3"),
                              #shinyjs::hidden(textInput("formula", "Formula ", "timeFactor ~ FactorA * FactorB")),
                              shinyjs::hidden(actionButton("infoButton", "", icon = icon("info-circle"))),
                              tippy::tippy_this("infoButton", "Example:<br><br>\n - 1 Factor: <br>\n  time ~ factorA <br><br>\n
                                    - 2 Factors:\n <br> time ~ factorA * factorB \n <br> or\n <br> time ~ factorA + factorB  ",
                                                placement = "right")),
                  shinyjs::hidden(h5(id = "titleWeights", strong("Further Arguments"), style = "color:grey")),
                  splitLayout(cellWidths = c("30%"), shinyjs::hidden(selectInput("hyp_mat", "Select contrast matrix:", c(`Dunnett` = "Dunnett", `Tukey` = "Tukey", `Center` = "center", `Other` = "Other"), selected = "Tukey"))),
                  # Other
                  splitLayout(cellWidths = c("70%","30%"), uiOutput(outputId = "out_mat"), uiOutput(outputId = "out_vec")),
                  conditionalPanel(condition = "input.hyp_mat == 'Other'",
                                   splitLayout(cellWidths = c("50%","50%"),
                                               actionButton(inputId = "add", label = "Add matrix"),
                                               conditionalPanel(condition = "input.add > input.del",actionButton(inputId = "del", label = "Delete matrix"))
                                   )),

                  splitLayout(cellWidths = c("15%"), shinyjs::hidden(numericInput("tau", "Endpoint tau of the relevant time window [0,tau]", 1))),
                  splitLayout(cellWidths = c("40%", "40%"), shinyjs::hidden(autonumericInput("alpha", "Level of significance alpha", value = 0.05)),
                              shinyjs::hidden(numericInput("nres", "Number of resampling repetitions", value = 4999))),
                  shinyjs::hidden(actionButton("process", "Calculate", class = "btn-primary")), width = 6),
                  mainPanel(verbatimTextOutput("group"), verbatimTextOutput("result"), plotOutput("result_plot"), width = 6)))

  server <- function(input, output, session){

    ### Part 1 - Data Loading ###
    # Purpose: This reactive function reads a CSV file uploaded by the user in a Shiny app.

    datasetInput <- reactive({
      req(input$infile)                                                         # Check if a file is uploaded
      if (is.null(input$infile))                                                # Return NULL if no file is uploaded
        return(NULL)
      read.csv(input$infile$datapath, header = input$header,                    # Read the uploaded CSV file with user-specified options
               sep = as.character(input$sep))
    })

    # This allows the Shiny app to dynamically read and process a CSV file based on user inputs.

    ### Part 2 - Visibility ###

    # Purpose: Controls the visibility of various UI elements based on file upload status and a checkbox.

    observeEvent(input$infile, {
      if (is.null(input$infile)) {                                              # Hide UI elements if no file is uploaded
        shinyjs::hide(id = "Method")
        shinyjs::hide(id = "formula")
        shinyjs::hide(id = "infoButton")
        shinyjs::hide(id = "hyp_mat")
        shinyjs::hide(id = "add")
        shinyjs::hide(id = "del")
        shinyjs::hide(id = "alpha")
        shinyjs::hide(id = "nres")
        shinyjs::hide(id = "tau")
        shinyjs::hide(id = "process")
        shinyjs::hide(id = "titleWeights")
        shinyjs::hide(id = "stepwise")
        shinyjs::hide(id = "plots")
      }
      else {                                                                    # Show UI elements if a file is uploaded
        shinyjs::show(id = "Method")
        shinyjs::show(id = "formula")
        shinyjs::show(id = "infoButton")
        shinyjs::show(id = "stepwise")
        shinyjs::show(id = "process")
        shinyjs::show(id = "titleWeights")
        shinyjs::show(id = "hyp_mat")
        shinyjs::show(id = "tau")
        shinyjs::show(id = "alpha")
        shinyjs::show(id = "nres")
        shinyjs::hide(id = "titleLoadData")
        observeEvent(input$stepwise, {                                          # Observe changes in the stepwise checkbox
          if (input$stepwise == TRUE){                                          # Hide the plot option if stepwise is TRUE
            shinyjs::hide(id = "plots")
            updateCheckboxInput(inputId = "plots", value = FALSE)
          }
          else{                                                                 # Show the plot option if stepwise is FALSE
            shinyjs::show(id = "plots")
          }
        })
      }
    })

    # Ensures that only relevant UI elements are displayed based on the file upload status and user selections.


    ### Part 3 - Dynamic UI Elements ###

    ## Part 3 a) - Dropdown for Presets ##

    # Purpose: Creates a dynamic dropdown and stores the selected value based on a condition.

    values <- reactiveValues()                                                              # Create a reactive object to store dynamic values

    # Dynamic Dropdown #
    output$dynamicInput <- renderUI({                                                       # Dynamically generate a dropdown menu
      selectInput(inputId = "dynamic", label = "Name of censoring status variable",
                  choices = colnames(datasetInput()))                                       # Choices are the column names of the dataset
    })

    # Reactive Observation #
    observe({ if (input$Method == "asymptotic" || input$Method == "permutation") {
      values$dyn <- input$dynamic                                                                                     # Store the selected value in values$dyn
    }
      else {                                                                                # Reset values$dyn if the method is not one of the specified values
        values$dyn <- NULL
      }
    })

    # Stores a dynamically selected value from the dropdown based on the chosen method (asymptotic, permutation).

    ## Part 3 b) - Setting (In)dependent Variables ##

    # Purpose: Dynamically creates UI elements and stores their values based on user inputs and dataset properties.

    values2 <- reactiveValues()                                                             # Create another reactive object to store dynamic values

    # Dropdown for Censored Variable #
    output$dynamicInput2 <- renderUI({
      selectInput(inputId = "dynamic2", label = "Label of censored variable",               # Create a dropdown for the censored variable
                  choices = sort(unique(datasetInput()[, values$dyn])),
                  selected = suppressWarnings(min(unique(datasetInput()[, values$dyn]))))
    })

    # Text Input for Formula #
    output$dynamicInput3 <- renderUI({                                                      # Dynamically generate a text input field for the formula
      if(length(grep("time", colnames(datasetInput())[colnames(datasetInput()) != input$dynamic])) > 0){
        textInput(inputId = "dynamic3", label = "Formula",
                  value = suppressWarnings(paste(colnames(datasetInput())[colnames(datasetInput()) != input$dynamic][(grep("time", colnames(datasetInput())[colnames(datasetInput()) != input$dynamic])[1])],
                                                 "~", paste(colnames(datasetInput())[colnames(datasetInput()) != input$dynamic][-(grep("time", colnames(datasetInput())[colnames(datasetInput()) != input$dynamic])[1])], collapse = " * ", sep = " ")) ))
      }
      else{
        textInput(inputId = "dynamic3", label = "Formula",
                  value = suppressWarnings(paste(colnames(datasetInput())[colnames(datasetInput()) != input$dynamic][1],
                                                 "~", paste(colnames(datasetInput())[colnames(datasetInput()) != input$dynamic][-1], collapse = " * ", sep = " ")) ))
      }
    })

    # Store Input Values #
    observe({                                                                               # Observe changes in input$dynamic2 and input$dynamic3
      values2$dyn <- input$dynamic2
      values2$dyn3 <- input$dynamic3
    })

    observe({
      # Ensure that dynamic3 (the formula) is not empty before proceeding
      req(input$dynamic3)

      # Load data and process dynamic inputs
      data <- as.data.frame(datasetInput())
      event <- input$dynamic
      formula <- isolate(input$dynamic3)

      # Use tryCatch to handle potential errors in num_factor_levels
      tryCatch({
        # Try to get the number of factor levels
        factor_levels <- num_factor_levels(formula, event, data)

        # Check if the condition is met (2 levels for 2 factors)
        if (length(factor_levels) == 2 && all.equal(as.integer(c(2, 2)), factor_levels)) {

          # If condition is met, add the "2by2" and "2by2 cause-wisely" options
          updateSelectInput(session, "hyp_mat",
                            choices = c(`Dunnett` = "Dunnett",
                                        `Tukey` = "Tukey",
                                        `Center` = "center",
                                        `Other` = "Other",
                                        `2by2` = "2by2",
                                        `2by2 cause-wisely` = "2by2 cause-wisely"),
                            selected = "Tukey")
        } else {
          # If condition is not met, show the default options
          updateSelectInput(session, "hyp_mat",
                            choices = c(`Dunnett` = "Dunnett",
                                        `Tukey` = "Tukey",
                                        `Center` = "center",
                                        `Other` = "Other"),
                            selected = "Tukey")
        }
      }, error = function(e) {
        # If num_factor_levels throws an error, show the default options
        updateSelectInput(session, "hyp_mat",
                          choices = c(`Dunnett` = "Dunnett",
                                      `Tukey` = "Tukey",
                                      `Center` = "center",
                                      `Other` = "Other"),
                          selected = "Tukey")
      })
    })

    # Allows a dynamic and interactive UI that responds to user inputs and dataset structure.


    ## Part 3 c) - 100 Hidden Matrices ##

    # Purpose: Generates dynamic UI elements (matrix inputs) based on user inputs and dataset properties.

    values3 <- reactiveValues()                                                             # Create another reactive object to store dynamic values

    # Define a Dynamic UI Element #
    output$out_mat <- renderUI({                                                            # Dynamically generate UI elements
      if (input$hyp_mat == "Other") {

        # Load and Prepare Dataset #
        data <- as.data.frame(datasetInput())                                               # Load and convert dataset to a data frame
        event <- data[, input$dynamic]                                                      # Extract specific column based on input$dynamic
        censored <- (event == input$dynamic2)
        data[!censored, input$dynamic] <- as.numeric(droplevels(factor(event[!(censored)])))
        data[censored, input$dynamic]  <- 0

        # Calculate Number of Unique Groups #
        k    <- length(unique(formula2input(input$dynamic3, input$dynamic, data)$group))    # Calculate unique groups based on the input formula

        # Generate Matrix Input Fields #
        lapply(1:100, function(i){                                                          # Create 100 hidden matrix input fields
          shinyjs::hidden( matrixInput(inputId = paste0("new_mat", i), label = paste0("H_",i," ="),
                                       value = matrix("", 1, k), rows = list(extend = TRUE, names = FALSE),
                                       cols = list(n = k, names = FALSE, extend = TRUE, delta = 0), class = "numeric"))
        })
      } else {
        return(NULL)
      }
    })

    # Generates dynamic matrix input fields based on user inputs and dataset characteristics, only if a specific condition is met.

    ### Part 4 - Matrix Input Fields ###

    # Purpose: Creates a UI with a list of matrix input fields, shown only if input$hyp_mat == "Other". These fields are hidden when created.

    output$out_vec <- renderUI({
      vecinputs <- lapply(1:100, function(i){
        if (input$hyp_mat == "Other") {
          return(matrixInput(inputId = paste0("new_vec", i), label = paste0("c_", i, " ="),
                             value = matrix(rep("", 1)),
                             rows = list(extend = TRUE, names = FALSE),
                             cols = list(n = 1, names = FALSE, extend = TRUE,
                                         delta = 0), class = "numeric"))
        } else {
          return(NULL)
        }
      })
      if (input$hyp_mat == "Other") shinyjs::hidden(vecinputs[!sapply(vecinputs, is.null)])
    })

    # Summary: Generates up to 100 hidden matrix input fields based on the condition input$hyp_mat == "Other".


    ### Part 5 - Add and Delete ###

    # Purpose: Dynamically shows or hides matrix and vector input fields based on user interactions with input$add and input$del.

    observeEvent(eventExpr = input$add, handlerExpr = {
      shinyjs::show(id = paste0("new_mat", input$add - input$del))
      shinyjs::show(id = paste0("new_vec", input$add  - input$del))
    })
    observeEvent(eventExpr = input$del, handlerExpr = {
      shinyjs::hide(id = paste0("new_mat", input$add - input$del + 1))
      shinyjs::hide(id = paste0("new_vec", input$add  - input$del + 1))
    })


    # Summary: Controls the visibility of input fields based on add/delete actions.


    observe({
      if (input$hyp_mat == "Other") {
        for (i in 1:(input$add - input$del)) {
          values3[[paste0("new_mat", i)]] <- input[[paste0("new_mat", i)]]
          values3[[paste0("new_vec", i)]] <- input[[paste0("new_vec", i)]]
        }
      } else {
        for (i in 1:(input$add - input$del)) {
          values3[[paste0("new_mat", i)]] <- NULL
          values3[[paste0("new_vec", i)]] <- NULL
        }
      }
    })

    # Summary: Ensures values3 holds the current input values or NULL based on input$hyp_mat.


    ### Part 7 - Error or Execution ###

    # Purpose: Monitors user interaction with input$process, performs validations, and executes hypothesis tests, displaying results.

    ## Part 7 a) - Error Handling ##

    observeEvent(input$process, {  # Triggered when `input$process` is activated.

      # Check for exactly two or three unique censoring types #
      if ((length(unique(as.data.frame(datasetInput())[, input$dynamic])) < 2)) {
        output$result <- renderPrint({
          "ERROR: Less than two censoring types"
        })
      }else{

        # Check if at least one contrast matrix is specified #
        if(input$hyp_mat == "Other" && input$add - input$del == 0){
          output$result <- renderPrint({
            "ERROR: Specify at least one contrast matrix."
          })}
        else {

          # Ensure all matrix rows are fully filled #
          if(input$hyp_mat == "Other" && any(unlist(sapply(1:(input$add - input$del),
                                                           function(i){
                                                             (apply(input[[paste0("new_mat",i)]], 1, function(row){
                                                               if(length(row) > 1){
                                                                 return(anyNA(row) && !all(is.na(row)))
                                                               }else{return(FALSE)}
                                                             }) )
                                                           })))){
            output$result <- renderPrint({
              "ERROR: Matrices are not completely filled in."
            })}else{

              # Check for empty matrices #
              if(input$hyp_mat == "Other" && any((sapply(1:(input$add - input$del), function(i) (nrow(na.omit(input[[paste0("new_mat",i)]])) == 0))))){
                output$result <- renderPrint({
                  "ERROR: Delete empty matrices."
                })}
              else {

                # Ensure vector lengths match matrix rows #
                if(input$hyp_mat == "Other" && any(sapply(1:(input$add - input$del),
                                                          function(i){
                                                            (nrow(na.omit(input[[paste0("new_mat",i)]])) !=
                                                             length(na.omit(input[[paste0("new_vec",i)]])))}))  ){
                  output$result <- renderPrint({
                    "ERROR: The vectors c must have the same length as number of rows in H."
                  })}else{

                    # Create and validate contrast matrices and vectors
                    if(input$hyp_mat == "Other"){
                      my_hyp_mat <- lapply(1:(input$add - input$del), function(i) na.omit(input[[paste0("new_mat",i)]]))
                      my_hyp_vec <- lapply(1:(input$add - input$del), function(i) na.omit(input[[paste0("new_vec",i)]]))
                      solution <- logical(length(my_hyp_mat))
                      for(i in 1:length(my_hyp_mat)){
                        solution[i] <- existence(my_hyp_mat[[i]], my_hyp_vec[[i]], input$tau)
                      }
                    }
                    else{
                      my_hyp_mat <- input$hyp_mat
                      my_hyp_vec <- NULL
                      solution <- TRUE
                    }
                    if(any(!solution)){
                      output$result <- renderText({
                        paste0("ERROR: Hypothesis", which(!solution),
                               " does not have a possible solution in [0,tau]^k.")
                      })

                      # Summary: If no solution exists for the hypothesis, an error is returned.


                      ## Part 7 b) - Data transformation and result display ##

                      # Purpose: Transforms data and performs the selected hypothesis test (asymptotic, permutation).
                      #   Displays a modal dialog during computation and shows results and plots.

                    } else{
                      data <- as.data.frame(datasetInput())
                      event <- data[, input$dynamic]
                      censored <- (event == input$dynamic2)
                      data[!censored, input$dynamic] <- as.numeric(droplevels(factor(event[!(censored)])))
                      data[censored, input$dynamic]  <- 0

                      output$group <- renderPrint({
                        cat("Group assignment: \n")
                        print(formula2groups(formula = isolate(input$dynamic3),
                                             event = input$dynamic,
                                             data = isolate(data),
                                             cens = input$dynamic2),
                              row.names = FALSE)
                      })



                      showModal(modalDialog("Calculating!"))

                      if (input$Method == "asymptotic") {
                        output_asy <- RMTL.asymptotic.test(formula = isolate(input$dynamic3), event = input$dynamic, Nres = input$nres,
                                                           data = isolate(data), hyp_mat = my_hyp_mat, hyp_vec = my_hyp_vec,
                                                           tau = input$tau, stepwise = input$stepwise, alpha = input$alpha)
                        removeModal()
                        output$result <- renderPrint({
                          GFDrmst:::summary.GFDrmst(output_asy)
                        })
                        if (input$plots) {
                          output$result_plot <- renderPlot({
                            GFDrmst:::plot.GFDrmst(output_asy)
                          })
                        }
                      }
                      # if (input$Method == "groupwise") {
                      #   output_grp <- RMTL.groupwise.test(formula = isolate(input$dynamic3),
                      #                                     event = input$dynamic,
                      #                                     data = isolate(data), hyp_mat = my_hyp_mat, hyp_vec = my_hyp_vec,
                      #                                     tau = input$tau, stepwise = input$stepwise, alpha = input$alpha)
                      #
                      #   removeModal()
                      #   output$result <- renderPrint({
                      #     GFDrmst:::summary.GFDrmst(output_grp)
                      #   })
                      #   if (input$plots) {
                      #     output$result_plot <- renderPlot({
                      #       GFDrmst:::plot.GFDrmst(output_grp)
                      #     })
                      #   }
                      # }
                      if (input$Method == "permutation") {
                        output_perm <- RMTL.permutation.test(formula = isolate(input$dynamic3),
                                                             event = input$dynamic, Nres = input$nres,
                                                             data = isolate(data), hyp_mat = my_hyp_mat, hyp_vec = my_hyp_vec,
                                                             tau = input$tau, stepwise = input$stepwise, alpha = input$alpha)

                        removeModal()
                        output$result <- renderPrint({
                          GFDrmst:::summary.GFDrmst(output_perm)
                        })
                        if (input$plots) {
                          output$result_plot <- renderPlot({
                            GFDrmst:::plot.GFDrmst(output_perm)
                          })
                        }
                      }
                    }
                  }}}}}
      #}
    })
  }

  # Summary: Monitors input$process, performs validations, and executes the appropriate hypothesis test, displaying the results.


  ### Part 8 - Shiny App ###

  # Purpose: Execute the function

  shinyApp(ui = ui, server = server)
}


