Sys.setenv(R_HISTFILE = file.path(getwd(), ".Rhistory"))

# Ensure history is saved on session exit
.Last <- function() {
	if (interactive()) try(savehistory(), silent = TRUE)
}

# Prevent restoring workspace on startup
options(save.workspace = FALSE)
