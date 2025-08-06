library(dplyr)
seed <-read.table("~/Desktop/CLY25_SEED_COUNT.TAB", sep= "\t", header=TRUE)
count <- with(seed, table(DESC,PACKED)) %>% as.data.frame()
count %>%
  mutate(PACKED = as.character(PACKED) %>% as.numeric(),
         total = PACKED*Freq)

     