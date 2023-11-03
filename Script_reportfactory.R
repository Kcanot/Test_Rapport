#Création de rapports
library(reportfactory)

new_factory(
  factory = "qPCR_TP_TAD",
  path = ".",
  report_sources = "report_sources",
  outputs = "outputs",
  move_in = TRUE,
  create_README = TRUE,
  create_example_report = TRUE,
  create_data_folders = TRUE,
  create_scripts_folder = TRUE,
  use_here = TRUE,
  use_rproj = TRUE,
  create_gitignore = TRUE
)

#Compiler = exécuter les rapports

compile_reports("example_report")

  #A partir de sous-dossiers
  compile_reports(
    reports = "summary_for_partners.Rmd",
    subfolder = "for_partners")
  
  #Avec paramétrisation
  compile_reports(
    reports = "daily_sitrep.Rmd",
    params = list(most_recent_data = TRUE,
                  region = "NORTHERN",
                  rates_denominator = 10000),
    subfolder = "regional"
  )

  #Visualisation
  factory_overview()