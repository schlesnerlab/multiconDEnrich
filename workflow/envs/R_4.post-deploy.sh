#!/bin/bash
# Install packages not available or problematic in conda
Rscript -e 'install.packages(c("DendSer", "mcprogress","glmmSeq"), upgrade = "never", dependencies = F, repos = "https://cloud.r-project.org")'
Rscript -e 'BiocManager::install(c("OmnipathR"), update = F, ask = F)'
Rscript -e 'devtools::install("workflow/scripts/RNAscripts", upgrade = "never", dependencies = F)'
