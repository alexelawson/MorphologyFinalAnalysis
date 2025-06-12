
library(dplyr)
library(stringr)

data_frame_final <- read.csv("/Users/alexlawson/Desktop/Masters-Work/data-frames/processed-dataframe.csv")

# Step 1: Get unique combinations (PVN as example, but can be changed)
combinations <- data_frame_final %>%
  filter(BrainRegion == "HYPO") %>%
  distinct(MouseID, SlideNumber, Treatment, Sex, BrainRegion)

# Step 2: Loop through each row
for (i in 1:nrow(combinations)) {
  mouse_id <- combinations$MouseID[i]
  slide_number <- combinations$SlideNumber[i]
  treatment <- combinations$Treatment[i]
  sex <- combinations$Sex[i]
  region <- combinations$BrainRegion[i]
  
  # Step 3: Filter original dataframe for this combo
  colorbycluster <- data_frame_final %>%
    filter(BrainRegion == region,
           MouseID == mouse_id,
           SlideNumber == slide_number) %>%
    dplyr::select(Cluster, ID)
  
  # Step 4: Construct file name
  filename <- paste0(region, "_", mouse_id, "_", treatment, "_", sex, "_", slide_number, ".csv")
  filepath <- file.path("/Users/alexlawson/Masters-Data-Final/Representative-images/Ciernia Lab/Hypothalamus/dapionly/v6/csv", filename)
  
  # Step 5: Write to CSV
  write.csv(colorbycluster, filepath, row.names = FALSE)
}
