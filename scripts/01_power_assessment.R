(result_file_path <- rev(dir("results/power/"))[1])
pg <- readRDS(paste0("results/power/",result_file_path))

pg |> glimpse()

pg |> unnest(out) %>%
 group_by(N) %>%
 summarise(
   power = mean(detected, na.rm = TRUE),
   .groups = "drop"
 )
