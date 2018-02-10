
dropdownButton <- function(label = "", status = c("default", "primary", "success", "info", "warning", "danger"), ..., width = NULL) {

	  status <- match.arg(status)
	  # dropdown button content
	  html_ul <- list(
		class = "dropdown-menu",
		style = if (!is.null(width)) 
		  paste0("width: ", validateCssUnit(width), ";"),
		#lapply(X = list(...), FUN = tags$li, style = "margin-left: 10px; margin-right: 10px;")
		lapply(X = list(...), FUN = tags$li, style = "margin: auto;")
	  )
	  # dropdown button apparence
	  html_button <- list(
		class = paste0("btn btn-", status," dropdown-toggle"),
		type = "button", 
		`data-toggle` = "dropdown"
	  )
	  html_button <- c(html_button, list(label))
	  html_button <- c(html_button, list(tags$span(class = "caret")))
	  # final result
	  tags$div(
		class = "dropdown", 
		do.call(tags$button, html_button),
		do.call(tags$ul, html_ul),
		tags$script(
		  "$('.dropdown-menu').click(function(e) {
		  e.stopPropagation();
			});"
		)
	  )
}