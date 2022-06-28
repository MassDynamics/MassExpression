#' Generate QC report
#'
#' @param listResults list. List of results from `runGenericDiscovery` 
#' @param output_folder str. Output folder path 
#' @param format str. 'pdf' or 'html'
#'
#' @return NULL
#' @export generate_qc_report
#'

generate_qc_report <- function(listResults, output_folder, format="html"){
  
  if(format == 'html'){
    # Render and save QC report
    dir.create(output_folder, showWarnings = FALSE)
    
    qc_report <- system.file("rmd","QC_report.Rmd", package = "MassExpression")
    rmarkdown::render(qc_report,
                      params = list(listInt = listResults,
                                    experiment = "Mass Dynamics QC report",
                                    output_figure = file.path(output_folder, "figure_html/"),
                                    format = "html"),
                      output_file = file.path(output_folder, "QC_Report.html"),
                      output_format=rmarkdown::html_document(
                        self_contained=FALSE,
                        lib_dir=file.path(output_folder,"qc_report_files"),
                        code_folding= "hide",
                        theme="united",
                        toc = TRUE,
                        toc_float = TRUE,
                        fig_caption= TRUE,
                        df_print="paged"))
  }else{
    # Render PDF
    output_folder_pdf <- file.path(output_folder, "pdf")
    dir.create(output_folder_pdf, showWarnings = FALSE)
    
    rmarkdown::render(qc_report,
                      params = list(listInt = listResults,
                                    experiment = "Mass Dynamics QC report",
                                    output_figure = file.path(output_folder_pdf, "figure_pdf/"),
                                    format = "pdf"),
                      output_file = file.path(output_folder_pdf, "QC_Report.pdf"),
                      output_format=rmarkdown::pdf_document(
                        toc = TRUE,
                        fig_caption= TRUE))
    }
  
}


#' Generate separate QC reports
#'
#' @param listResults list. List of results from `runGenericDiscovery` 
#' @param output_folder str. Output folder path 
#'
#' @return NULL
#' @export generate_separate_qc_reports
#'

generate_separate_qc_reports <- function(listResults, output_folder){
  
  # Render separate QCs
  qc_names <- get_names_qc()
  
  for(name in qc_names){
    print(paste0("Creating QC for: ",name))
    qc_report_name <- paste0("QC_", name, ".Rmd")
    qc_report_output <- paste0("QC_", name, ".html")
    qc_report <- system.file("rmd", qc_report_name, package = "MassExpression")
    
    rmarkdown::render(qc_report,
                      params = list(listInt = listResults,
                                    output_figure = file.path(output_folder, "figure_html_separate/")),
                      output_file = file.path(output_folder, qc_report_output),
                      output_format = rmarkdown::html_document(
                        self_contained=FALSE,
                        lib_dir = file.path(output_folder,"qc_report_files"),
                        theme = "united",
                        fig_caption = TRUE,
                        df_print = "paged"))
  }
}