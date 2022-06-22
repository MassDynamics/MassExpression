#' Sanitize strings in dataframe
#' 
#' @param df data frame
#' 
#' @description Sanitize all character entries in dataframe. sanitize from invisible and unicode characters.
#' @export sanitize_strings_in_dataframe

sanitize_strings_in_dataframe <- function(df){
  characters_columns <- colnames(df)[sapply(df, class) == "character"]
  df_subset <- df[, characters_columns, drop=FALSE]
  
  if(ncol(df_subset) > 0){
    sanitize_vector <- Vectorize(sanitize_strings, vectorize.args = "any_string")
    
    df_subset_cleaned <- apply(df_subset, 2, sanitize_vector)
   
    df[, characters_columns] <- df_subset_cleaned
  }
  
  return(df)
} 


#' Sanitize strings for export
#' 
#' @param any_string string
#' 
#' @description Remove invisible characters from string and substitute Unicode characters with empty space
#' @export sanitize_strings

sanitize_strings <- function(any_string){
  # Substitute unicode with empty space
  # \u00A0 represents no break?
  any_string <- gsub("\u00A0", " ", any_string)
  # \xa0 represents a hard break space
  any_string <- gsub("\xa0", " ", any_string)
  
  # Remove invisible characters
  chars <- strsplit(any_string, "")[[1]]
  chars <- chars[nchar(chars, type = "width") > 0]
  good_string <- paste(chars, collapse = "")
  
  good_string <- trimws(good_string)
  
  return(good_string)
}




