

survData <- readRDS(file.path("F:\\Zhaolab2020\\virome\\TGScompare\\02genome_bins", "survData.rda"))

survfitted <- survfit(Surv(time = stim_at_thresh, event = does_respond) ~ 1, 
                      data = survData, 
                      conf.int = FALSE)

ggsurvplot(fit = survfitted,
           fun = "events",
           conf.int = FALSE)