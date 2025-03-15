#############################
### 00 CREATE DIRECTORIES ###
#############################

#renv::init()
#renv::update()
#renv::snapshot()

# Create data directory
data_dir <- dir.create("data")

# create data subdirectories
data_subdirectories <- c("raw",
                         "processed",
                         "metadata")
for(i in 1:length(data_subdirectories)) {
  subdirectories <- dir.create(paste0("data/", data_subdirectories[i]))
}

# Create script directory
scripts_dir <- dir.create("scripts")

# Create results directory
results_dir <- dir.create("results")
results_subdirectories <- c("figures",
                           "tables")

for(i in 1:length(results_subdirectories)){
  subdirectories <- dir.create(paste0("results/", results_subdirectories[i]))
}

# functions directory
function_dir <- dir.create("functions")

# references directory
references_dir <- dir.create("references")


## if you want to delete a directory
#dir_path <- "data"
#if (dir.exists(dir_path)) {
#  unlink(dir_path, recursive = TRUE)
#  message("Directory deleted successfully.")
#} else {
#  message("Directory does not exist.")
#}
